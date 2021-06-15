/*
  Fuse given HDF *coord*.h5 files into one monolithic HDF
  Type 'Converter' (w/o arguments) for usage.
  last change: 06may14 JK
 */

#include "Converter.H"

Attribute::Attribute (string attr_name) {
	name = attr_name;
}

void Attribute::fill_s (int attr_size, string value) {
	type = 1;
	size = 1;
	svalue = new string [1];
	svalue[0] = value;
}

template<class T>
void Attribute::fill_n (int attr_size, T value[], int attr_type) {
	type = attr_type;
	size = attr_size;
	if (type == 2) {
		ivalue = new int [size];
		for (int pos = 0; pos < size; pos++)
			ivalue[pos] = int(value[pos]);
	} else if (type == 3) {
		fvalue = new float [size];
		for (int pos = 0; pos < size; pos++)
			fvalue[pos] = float(value[pos]);
	} else if (type == 4) {
		dvalue = new double [size];
		for (int pos = 0; pos < size; pos++)
			dvalue[pos] = value[pos];
	}
}


// bool file_exists (string fname)
// {
// 	FILE *file;
// 	file = fopen (fname.c_str(), "r");
// 	if (file != NULL) {
// 		fclose(file);
// 		return true;
// 	}
// 	return false;
// }


bool Converter::findfname (string &filename, int step, int nx, int ny, int nz)
{
	// build filename string and see if it exists
	
	string cstep_str = "error";
	if (step != -1) {
		std::stringstream ss;
		ss << step;
		cstep_str = ss.str();
	}
	string cco_str = "";
	if (nx > -1 and ny > -1 and nz > -1) {  // {-1,-1,-1} => omit
		char ccoord[255];                    //  coord part of filename
		sprintf (ccoord, "_coord%i-%i-%i", nx,ny,nz);
		cco_str = string(ccoord);
	}
	filename = poub + "/" + pname + "_float" + "/" +
		pname + "_flt_step" + cstep_str + cco_str + ".h5";
	
	return file_exists (filename);
}


Converter::Converter (string req_poub, string req_pname)
{
	this->poub  = req_poub;
	this->pname = req_pname;
	if (file_exists (poub +"/"+ pname +".cat")) {
		setenv ("poub" ,   poub.c_str(), 1);
		setenv ("pname" , pname.c_str(), 1);
		EdgeGridding = (value_exists((char*)"Nx") and
		                value_exists((char*)"Ny") and
		                value_exists((char*)"Nz"));
	} else {
		cout << "No cat file found. Assume EdgeGridding?\n";
		string choice = "";
		while (true) {
			cout << " (y)es / (n)o / (a)bort : ";
			std::cin >> choice;
			if      (choice == "y") { EdgeGridding = true;  break; }
			else if (choice == "n") { EdgeGridding = false; break; }
			else if (choice == "a") { exit(-18); }
		}
	}
	
	for (short dir = 0; dir < 3; ++dir)  // to be filled on
		{ nproc[dir] = 0; }                 //  first conversion call
}


void Converter::convert (int itime)
{
	if (nproc[0]*nproc[1]*nproc[2] == 0) { // first call?
		int np[] = {0,0,0};                 //  -> get nproc{x|y|z}
		string fdum; fdum = fdum;
		for (np[0] = 0; findfname (fdum, itime, np[0], 0, 0);) { np[0]++; }
		for (np[1] = 0; findfname (fdum, itime, 0, np[1], 0);) { np[1]++; }
		for (np[2] = 0; findfname (fdum, itime, 0, 0, np[2]);) { np[2]++; }
		for (short dir = 0; dir < 3; ++dir)
			{ nproc[dir] = np[dir]; }
	}
	
	if (verbose)
		cout << "\n"
		     << "---------------------------------"
		     << " Conversion for timestep: " << itime
		     << "--------------------------------- \n";	
	
	string filename;    // Get filename of master rank
	bool test = findfname (filename, itime, 0,0,0);
	numFields = GetNumFields (filename);  // Get number of fields
	
	// Make as many fields as needed
	om = new NumMatrix<float,3>[numFields];
	
	for (int nx = 0; nx < nproc[0]; nx++)
		for (int ny = 0; ny < nproc[1]; ny++)
			for (int nz = 0; nz < nproc[2]; nz++) {
				test = findfname (filename, itime, nx,ny,nz);
				if (test)  ReadHDF (filename);
				else {
					cerr << "ERROR: Failure to open file '"
					     << filename << "'\n";
					exit(-5);
				}
			}
	WriteHDF (itime);
}


void Converter::get_itimeList (std::vector<int> &it_list)
{
	// --> build list of available files and store their itime parts
	DIR *dir;
	struct dirent *dirpointer;
	struct stat buf;
	const string dirname = poub + "/" + pname + "_float";
	
	if ((dir = opendir(dirname.c_str())) == NULL) {
		cerr << " Error while opening dir... '" << dirname << "'\n";
		exit(22);
	}
	if (stat(dirname.c_str(), &buf) == -1) {
		cerr << " Error for stat '" << dirname << "'\n";
		exit(22);
	}
	
	while ((dirpointer = readdir(dir)) != NULL) {
		string fname = (*dirpointer).d_name;
		if (fname.find("coord0-0-0") < fname.length()) {
			string fullname = dirname + "/" + fname;
			int sPos =  dirname.size()+pname.size()+10;
			int sLen = fullname.size()-sPos        -14;
			int itime = atoi (fullname.substr(sPos,sLen).c_str());
			it_list.push_back (itime);
		}
	}
}

int Converter::GetNumFields (string filename)
{
	Hdf5iStream h5file (filename, -1);
	int entries;
	h5file.ReadGlobalAttr ("Entries", entries);
	return entries;
}

void Converter::ReadHDF (string filename)
{
	int mxloc[3], shift[3], mxdata[3], coords[3];
	int rank, numCellsLoc[3];
	double xbloc[3];
	
	if (verbose)
		cout << "Opening file '" << filename << "'\n";
	
	if (attr_list.size() == 0) { // --> read and store attributes
		herr_t status;
		H5O_info_t object_info; 
		hid_t hdf5file = H5Fopen (filename.c_str(),
		                          H5F_ACC_RDONLY, H5P_DEFAULT);
		hid_t group = H5Gopen2 (hdf5file, "Data", H5P_DEFAULT);
		if (H5Oget_info (group, &object_info) < 0) {
			cerr << "H5 group access error! \n";
			exit(1);
		}
		status = H5Oget_info (group, &object_info);
		int N_attr = object_info.num_attrs; 
		for (int ia = 0; ia < N_attr; ia++) {
			hid_t attr_id = H5Aopen_by_idx (group, ".", H5_INDEX_CRT_ORDER,
			                                H5_ITER_INC, ia,
			                                H5P_DEFAULT, H5P_DEFAULT);  	 	
			// --> get attribute name
			const size_t maxlen = 40;
			char namebuf[maxlen];  // for name of current attribute
			ssize_t len = H5Aget_name (attr_id, maxlen, namebuf);
			string attr_name = string (namebuf);
			
			// --> get attribute datatype
			hid_t   attr_type = H5Tget_class (H5Aget_type (attr_id));
			size_t  attr_size = H5Tget_size  (H5Aget_type (attr_id));
			hsize_t attr_strs = H5Aget_storage_size (attr_id);
			Attribute attr (attr_name);
			int type = 0;
			if      (attr_type == H5T_STRING)  { type = 1; } // string
			else if (attr_type == H5T_INTEGER) { type = 2; } // int
			else if (attr_type == H5T_FLOAT)   { type = 3; } // float
			if (type == 3 and attr_size == 8)  { type = 4; } // double
			
			if (type != 0) {
				int arrsize = int (attr_strs / attr_size);
				if (type == 1) {
					char c_value [int(maxlen)];
					hid_t ntype = H5Tget_native_type (H5Aget_type (attr_id),
					                                  H5T_DIR_ASCEND);
					status = H5Aread (attr_id, ntype, c_value);
					attr.fill_s (1, string(c_value));
				} else if (type == 2) {
					int * value = new int[arrsize];
					status = H5Aread (attr_id, H5T_NATIVE_INT, value);
					attr.fill_n (arrsize, value, type);
				} else if (type == 3) {
					float * value = new float[arrsize];
					status = H5Aread (attr_id, H5T_NATIVE_FLOAT, value);
					attr.fill_n (arrsize, value, type);
				} else if (type == 4) {
					double * value = new double[arrsize];
					status = H5Aread (attr_id, H5T_NATIVE_DOUBLE, value);
					attr.fill_n (arrsize, value, type);
					if (attr_name == "time" and verbose)
						cout << "R: t=" << *value << std::endl;
				}
				attr_list.push_back (attr);
			}
		}
	}
	
	Hdf5iStream h5file (filename, -1);
	int entries;
	
	h5file.ReadGlobalAttr ("rim"    ,  rim    );
	h5file.ReadGlobalAttr ("rank"   ,  rank   );
	h5file.ReadGlobalAttr ("coords" , *coords );
	h5file.ReadGlobalAttr ("Entries",  entries);
	if (entries != numFields) {
		cerr << " ERROR: Number of entries differs: "
		     << entries << " - " << numFields << "\n";
		exit(2);
	}
	
	omNames = new string[numFields];
	for (int q = 0; q < numFields; q++) {
		
		string omName = h5file.GetDatasetName(q);
		omNames[q] = omName;
		
		h5file.getSize(omName, mxloc, 3);	// get size of data (w/o ghosts)
		for (short dir = 0; dir < 3; ++dir) {
			mxloc[dir] -= 2*rim+1;
			if (EdgeGridding) {
				numCellsLoc[dir] = mxloc[dir] + 1;
				mx[dir] = numCellsLoc[dir] * nproc[dir] - 1;
			} else {
				mx[dir] = mxloc[dir] * nproc[dir];
			}
		}
		
		NumMatrix<float,3> data;
		if (!h5file.Read3DMatrix(omName, data, xbloc, dx)) {
			cerr << " Reading not successful for om: " << omName << "\n";
		}
		for (short dir = 0; dir < 3; ++dir) {
			mxdata[dir] = data.getHigh(dir) - data.getLow(dir);
			shift[dir]  =     ( mxdata[dir] -       mxloc[dir] )/2;
		}
		if (shift[0] != rim or shift[1] != rim or shift[2] != rim) {
			cerr << "ERROR: rim mismatch detected!" << "\n";
			exit(-3);
		}
		
		if (verbose)
			cout << " Reading field q=" << q
			     << " : " << omName << " ...";
		if (rank == 0) {
			om[q].resize(Index::set(     -rim,     -rim,     -rim),
			             Index::set(mx[0]+rim,mx[1]+rim,mx[2]+rim));
			for (short dir = 0; dir < 3; ++dir)
				{ xb[dir] = xbloc[dir]; }   // save xmin from coord0-0-0
		}
		
		if (q == 0) 	// New computation of coords acc. to info in file
			for (short dir = 0; dir < 3; ++dir) {
				if (EdgeGridding) coords[dir] *= numCellsLoc[dir];
				else 					coords[dir] *= mxloc[dir];
			}
		
		for (int k = -rim; k <= mxloc[2]+rim; k++)
			for (int j = -rim; j <= mxloc[1]+rim; j++)
				for (int i = -rim; i <= mxloc[0]+rim; i++)
					om[q] (i+coords[0],
					       j+coords[1],
					       k+coords[2]) = data (i+shift[0],
					                            j+shift[1],
					                            k+shift[2]);
		if (verbose)  { cout << " done. \n"; }
	}
	if (verbose) { cout << " finished. \n"; }
}

void Converter::WriteHDF (int itime)
{
	string filename;
	bool fdum = findfname (filename, itime, -1,-1,-1);
	Hdf5Stream h5out (filename, numFields);
	cout << "Setting up output HDF5 \"" << filename << "\"\n";
	
	for (unsigned ia = 0; ia < attr_list.size(); ia++) {
		Attribute c_attr = attr_list.at(ia);
		string name = c_attr.name;
		if (name != "Entries" and name.substr(0,7) != "Name_om" and
		    name != "using_cbase")	{
			short type = c_attr.type;
			int   size = c_attr.size;
			if (name == "time" and verbose)
				cout << "W: t=" << *(c_attr.dvalue) << std::endl;
			if (size == 1) {
				if (type == 1)  h5out.AddGlobalAttr (name, *(c_attr.svalue));
				if (type == 2)  h5out.AddGlobalAttr (name, *(c_attr.ivalue));
				if (type == 3)  h5out.AddGlobalAttr (name, *(c_attr.fvalue));
				if (type == 4)  h5out.AddGlobalAttr (name, *(c_attr.dvalue));
			} else {
				if (type == 2)  h5out.AddGlobalAttr (name, c_attr.ivalue, size);
				if (type == 3)  h5out.AddGlobalAttr (name, c_attr.fvalue, size);
				if (type == 4)  h5out.AddGlobalAttr (name, c_attr.dvalue, size);
				}
			}
	}

	//	--> write 'heavy' data
	for (int q = 0; q < numFields; q++) {
		string omName = omNames[q];
		std::stringstream q_ss;
		q_ss << q;
		string outp_txt =
			" Writing field q=" + q_ss.str() + " : " + omName + " ... ";
		bool res_ok = h5out.Write3DMatrix(omName, om[q], xb, dx);
		if (res_ok) {
			if (verbose)
				cout << outp_txt + "done. \n";
		} else {
			cout << outp_txt + "FAIL! \n";
			exit (-15);
		}
	}
	
	for (unsigned ia = 0; ia < attr_list.size(); ia++) {
		Attribute c_attr = attr_list.at(ia);
		if (c_attr.type == 0)
			cout << "Warning: attribute '" << c_attr.name
			     <<  "' has unknown datatype, skipped. \n";
	}
	attr_list.clear();
}


int main (int argc, char *argv[])
{
	const string usage_text[] = {
		" ",
		"Possible options for converter usage:  ",
		" 1) Converter [<poub> <pname>] [-q] <itime> ",
		" 2) Converter [<poub> <pname>] [-q] --all   ",
		" 3) Converter [<poub> <pname>] [-q] --new   ",
		" 4) Converter [<poub> <pname>] [-q] --ask (-> prompt for selection) ",
		" 5) Converter [<poub> <pname>] [-q]       (-> just output this text) ",
		" (Use -q option to suppress non-critical output.) "
	};
	
	std::vector<string> arglist;
	
	bool quiet = false; // -> check for -q (quiet) option
	for (int ia = 1; ia < argc; ia++) {
		string arg = string(argv[ia]);
		if (arg.compare("-q") == 0)  { quiet = true; }
		else               { arglist.push_back(arg); }
	}
	
	int itime = -1;
	string poub, pname;
	bool env_ok = (getenv("poub") != NULL and getenv("pname") != NULL);
	
	if (arglist.size() == 1) {
		if (env_ok) {
			poub  = string (getenv("poub"));
			pname = string (getenv("pname"));
		} else {
			cerr << "ERROR: poub and/or pname not set! \n";
			exit(3);
		}
	} else if (arglist.size() == 3) {
		poub  = arglist[0];
		pname = arglist[1];
	} else {    // no valid args => print usage
		for (int line = 0; line < 8; line++)
			cerr << usage_text[line] << "\n";
		exit(-17);
	}
	
	const string arg = arglist.back();
	bool do_all = (arg == "--all"); // convert all existing coord files
	bool do_new = (arg == "--new"); // convert only those that lack joined *.h5
	bool do_ask = (arg == "--ask"); // ask user for <itime>
	if (!do_all and !do_new and !do_ask) {
		itime = atoi(arg.c_str());
		cout << "parsing result: itime = " << itime << "\n";
	}
	
	Converter conv (poub, pname);
	conv.verbose = !quiet;
	std::vector<int>::iterator iter;
	std::vector<int> it_list;   // list of itimes to convert
	
	if (itime == -1) {
		conv.get_itimeList (it_list);
		if (it_list.size() == 0)
			cerr << "ERROR: No coord files matched. \n";
		else {
			if (!do_all and !do_new) {     // -> show menu
				string filename;
				cout << " Please choose from available output times: \n";
				for (iter = it_list.begin(); iter != it_list.end(); ++iter) {
					double time_loc;
					bool tdum = conv.findfname (filename, *iter, 0,0,0);
					Hdf5iStream h5file (filename, -1);
					h5file.ReadGlobalAttr ("time", time_loc);
					cout << *iter << " at time " << time_loc << "\n";
				}
				cout << "\n";
				while (conv.findfname (filename, itime, 0,0,0) == false) {
					cout << " Your choice: ";
					std::cin >> itime;
				}
				it_list.clear();
				it_list.push_back (itime);
			}
		}
	} else {  // -> place <itime> from arg as only element in list
		string filename;
		if (conv.findfname (filename, itime, 0,0,0))
			it_list.push_back (itime);
		else
			cerr << "ERROR: '" + filename + "' not found. \n";
	}
	
	if (it_list.size() == 0)
		cerr << "Warning: nothing to do. \n";
	else {
		string filename;
		std::sort (it_list.begin(), it_list.end());
		for (iter = it_list.begin(); iter != it_list.end(); ++iter) {
			bool skip = false;
			if (do_new)
				{ skip = conv.findfname (filename, *iter, -1,-1,-1); }
			if (skip) { cout << "skipping file '" + filename + "'\n"; }
			else      { conv.convert (*iter); }
		}
	}
	return 0;
}
