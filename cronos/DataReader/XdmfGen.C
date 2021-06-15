/*
 *  PROGRAM:  XdmfGen V0.94-r
 *
 *  PURPOSE:  generate Xdfm "light data" to describe "heavy" HDF5 data
 *            as written by Cronos (mandatory for ParaView visualization)
 *
 *  REQUIRES:  text file "DataReader/xdmf_body.partxmf"
 *             (as defined via FN_XBODY in Makefile)
 *
 *  CHANGELOG:  Feb 2011  (initial)
 *              Jul 2011  output to single *.xmf
 *              Dec 2012  adaption to string size bug in HFA_read (minor)
 *              Jul 2014  extension to multi-fluid
 *              Jul 2018  == revision work pending ==
 *              Dec 2019  bugfix related to fixed/variable names in float.h5's
 *              Mar 2020  bugfix, abort due to faulty variable consistency error
 *  AUTHOR:  Jens Kleimann <jk@tp4.rub.de>
 *
 */

#define VERSION  "XdmfGen V0.95"
#define USAGE    "XdmfGen <poub> <pname> <itime_1> [<itime_2> [...] ]"
// #define FN_XBODY "xdmf_body.partxmf" // definition moved to Makefile

#include <iostream>
#include <vector>
#include <cstdlib>
#include "hdf5.h"
#include <fstream>
using namespace std;

class quant {  // name class for phys. quantities in HDF
public:
	string name, // name as in HDF ("v_z")
		para,     // as used for PV ("F0__v_rad") (multi-fluid, else w/o F part)
		path;     // path in HDF ("/fluid0/v_z")
};

bool operator==( quant qu1, quant qu2 ) {
	return ((qu1.name == qu2.name) and
	        (qu1.para == qu2.para) and
	        (qu1.path == qu2.path));
}

bool file_missing( string filename ) {
	FILE *file;
	if ((file = fopen( filename.c_str(), "r"))) {
		fclose( file );
		return false;
	} else {
		cerr << "Error, file \"" + filename + "\" not found!" << endl;
		return true;
	}
}

void abort( string message = "" ) {
	if (message != "") {
		cerr << "Error: " + message + ",";
	}
	cerr << " aborting." << endl;
	exit(1);
}

string remove_last2( string input ) {
	int len = input.size();
	if (len < 2) {
		abort( "string too short for removal." );
	}
	return input.substr( 0, len-2 );
}

int listpos( const string & str, const vector<quant> & vect ) {
	// --> return position of <str>'s first occurrence in <vect>.para
	for (int pos = 0; pos < int(vect.size()); ++pos)
		if (vect[pos].para == str)
			return pos;
	return -1;  // report item_not_found
}

void paste_line( FILE * xmf, string line ) {
	fprintf( xmf, (line + "\n").c_str() );
}

void paste_block( FILE * xmf, string tag,
						string ins1 = "", string ins2 = "", string ins3 = "" ) {
	// --> copy tagged <text> from stencil file, optionally replacing
	//     "%{1|2|3}%" by strings ins1, ins2, ins3 in the process 
	//     pattern: [tag]:<text>
	//              012345678
	string token[] = {"%1%", "%2%", "%3%"}; // patterns to locate in line
	string insrt[] = { ins1,  ins2,  ins3}; // what to replace them with
	ifstream xbody ( FN_XBODY );
	bool tag_found = false;
	if (xbody.is_open()) {
		while (! (xbody.eof()) ) {
			string line;
			getline( xbody, line );
			if (line.compare (0, 6, "["+ tag +"]:" ) == 0) { // tag found?
				tag_found = true;
				for (int is = 0; is < 3; is++) {
					size_t pos = line.find( token[is] );
					if (pos != string::npos)     // perform replacements
						line = line.substr(0, pos) + insrt[is] + line.substr(pos+3);
				}			
				fprintf( xmf, (line.substr(6)+"\n").c_str() );
			}
		}
		if (tag_found == false)
			cerr << "Warning, tag '" << tag << "' not found." << endl;
		xbody.close();
	} else
		abort( "missing XDMF stencil file '" + string(FN_XBODY) + "'" );
}


int main (int argc, char *argv[]) {
	cout << " ** " << VERSION << " ** " << endl;
	const int Nf = argc-3;   // arg: total number of frames
	if (Nf < 1) {
		abort( "Usage: " + string(USAGE) );
	}
	const string poub  = argv[1];
	const string pname = argv[2];
	const string gridfile = pname + "_grid/" + pname + "_grid.h5";
	const string xdmffile = pname            +         ".xmf"    ;
	const string x_dirfile = poub + "/" + xdmffile;
	const string g_dirfile = poub + "/" + gridfile;
	int idummy;   // all-purpose dummy integer return value
	
	// --> check for grid file & get name of coordinate system
	if (file_missing (g_dirfile)) {
		abort("mandatory grid file not found");
	}
	char *attr_value;        // buffer for names of fields
	hid_t id_file, id_group, id_rootgr, id_dset, id_space, id_attr;
	herr_t status;
	
	id_file   = H5Fopen ( g_dirfile.c_str(), H5P_DEFAULT, H5F_ACC_RDONLY );
	id_group  = H5Gopen2( id_file, "/Data" , H5P_DEFAULT );
	id_attr   = H5Aopen_name( id_group, "Geom_Type" );
	status    = H5Aread ( id_attr, H5Aget_type(id_attr), &attr_value );
	status    = H5Aclose( id_attr );	
	status    = H5Gclose( id_group ); 
	status    = H5Fclose( id_file );
	string geom_name = string(attr_value);
	
	// (Car)tesian label subscripts (may) get replaced by specific ones
	const string subs_Car[] = {  "x",   "y",   "z"};
	const string subs_Cyl[] = {  "r", "phi",   "z"};
	const string subs_Sph[] = {"rad", "tet", "phi"};
	const string subs_gen[] = {  "1",   "2",   "3"};
	vector<string> subs;
	if      (geom_name == "Cartesian"  ) { subs.assign( subs_Car, subs_Car+3 ); }
	else if (geom_name == "Cylindrical") { subs.assign( subs_Cyl, subs_Cyl+3 ); }
	else if (geom_name == "Spherical"  ) { subs.assign( subs_Sph, subs_Sph+3 ); }
	else                                 { subs.assign( subs_gen, subs_gen+3 ); }
	
	string * itime    = new string[Nf];  // itime strings (from argv[])
	string * time     = new string[Nf];  // simulation timestamps
	string * hdf5file = new string[Nf];  // HDF5 file names 
	int    * rim      = new int   [Nf];  // number of ghost cells
	
	hsize_t fulldimsZYX[3];  // (full) array dimensions
	vector<quant> vstruct;   // list of vector quantities (v, B, ...)
	vector<quant> qstruct;   // list of scalar quantities (rho, v_x, ...)
	vector<quant> qstruct0;  //  copy thereof (for consistency check)
	int Nfgrid[3] = {-1};    // grid extension (w/ ghosts)
	bool same_grid = true;
	bool same_vars = true;
	bool same_rim  = true;
	
	for (int it = 0; it < Nf; it++) {  // loop over HDF files
		itime   [it] = argv[3+it];
		hdf5file[it] = (pname + "_float/"   +
		                pname + "_flt_step" + itime[it] + ".h5");
		const string h_dirfile = poub + "/" + hdf5file[it];
		
		cout << " Opening data file '" << h_dirfile << "'...";
		if (file_missing (h_dirfile)) { abort( "Data file not found" ); }
		else                          { cout << "OK"; }
		
		id_file = H5Fopen( h_dirfile.c_str(), H5P_DEFAULT, H5F_ACC_RDONLY );
		
		// --> get { rim, time } from datafile
		double r_time;
		id_group = H5Gopen2     ( id_file, "/Data", H5P_DEFAULT );
		id_attr  = H5Aopen_name ( id_group, "rim" );
		status   = H5Aread      ( id_attr, H5Aget_type(id_attr), &(rim[it]) );
		status   = H5Aclose     ( id_attr );
		id_attr  = H5Aopen_name ( id_group, "time" );
		status   = H5Aread      ( id_attr, H5Aget_type(id_attr), &r_time );
		status   = H5Aclose     ( id_attr );
		status   = H5Gclose     ( id_group );
		
		const size_t maxlen = 20;
		char time_cbuf[maxlen]; // buffer for timestamp conv. (dbl -> string)
		idummy = sprintf( time_cbuf, "%f", r_time );
		time[it] = string( time_cbuf );
		cout << " (time = " << time[it] << ")" << endl;
		
		// --> generate list of all om fields in /Data/ or /fluid0, ...
		//     (both scalars and vector comp.s)
		H5G_info_t group_info;
		H5O_info_t obj_info;
		id_rootgr = H5Gopen2( id_file, "/", H5P_DEFAULT );
		status    = H5Gget_info( id_rootgr, &group_info );  // get # of groups
		int n_link = group_info.nlinks;
		qstruct.clear();
		
		vector<string> list_groups;  // -> build list of group names
		list_groups.push_back( "/Data" );  // always present
		for (int ilink = 0; ilink < n_link; ilink++) {  // iterate over links
			status = H5Oget_info_by_idx( id_file, "/", H5_INDEX_NAME,
			                             H5_ITER_INC, ilink, &obj_info, NULL );
			id_group =   H5Oopen_by_idx( id_file, "/",  H5_INDEX_NAME,
												  H5_ITER_INC, ilink, NULL ); 
			
			if (obj_info.type == H5O_TYPE_GROUP) {  // group found!
				char group_name[maxlen];   // -> get group name
				ssize_t len = H5Iget_name( id_group, group_name, maxlen );
				string group_string = string (group_name); 
				if (group_string.compare( 0, 5, "/flui" ) == 0) // fluid group
					list_groups.push_back( group_string );  // ("/Data" has 5 letters!)
			}
		} // -> list_groups now complete!
		
		if (list_groups.back() != "/Data") {          // -> remove /Data if
			list_groups.erase( list_groups.begin() );  //    not a fluid group
		}
		int n_flu = list_groups.size();
		for (int iflu = 0; iflu < n_flu; iflu++) { // iterate over fluid groups
			string group_string = list_groups[iflu];
			string flu = "0";   // default "/Data" case
			if (group_string != "/Data") {    // --> get fluid number (defined as
 				flu = group_string.substr(6);  //     whatever comes after "/fluid")
				if (n_flu <= 10) { flu = flu.substr(1); }  // ...w/o leading '0's
			}
			id_group = H5Gopen2( id_file, group_string.c_str(), H5P_DEFAULT );
			const int n_attr = H5Aget_num_attrs( id_group );
			char attr_name[maxlen];  // for name of current attribute
			for (int iattr = 0; iattr < n_attr; ++iattr) {
				id_attr = H5Aopen_idx( id_group, iattr );
				idummy  = int( H5Aget_name( id_attr, maxlen, attr_name ) );
				const string att_string = string(attr_name) + "      ";
				if (att_string.compare( 0, 7, string("Name_om") ) == 0) {
					hid_t id_memtype = H5Aget_type( id_attr );
					string field_string;   // name of field to add
					if (H5Tis_variable_str( id_memtype) ) {   // -> variable size
						status = H5Aread( id_attr, id_memtype, &attr_value );
						field_string = string(attr_value);
						cout << "VAR " << field_string << endl;
					} else {    // -> fixed size
						hsize_t size = H5Aget_storage_size( id_attr );
						char attr_value[size+1];
						status = H5Aread( id_attr, id_memtype, &attr_value );
						field_string = string(attr_value).substr(0, size);  // cut trailing rubbish
					}
					quant newquant;  // --> build quantity object
				
					string prefix = "";
					bool is_AorB = false; 
					const string AB_comp[] = {"A_x", "A_y", "A_z", "B_x", "B_y", "B_z"};
					for (int ic = 0; ic < 6; ic++) {
						is_AorB = is_AorB or (field_string.substr(0, 2) == AB_comp[ic]);
					}
					if ((n_flu > 1) and (is_AorB == false)) {
						prefix = "F" + flu + "__";
					}
					newquant.name = field_string;
					newquant.path = string(group_string) +"/"+ field_string;
					newquant.para = prefix + field_string;
					qstruct.push_back( newquant );     // add quantity to list
				}
				status = H5Aclose( id_attr );
			}
			status = H5Gclose( id_group );
		}
		status = H5Gclose( id_rootgr );
		
		if (it == 0) {    // output this only once
			if (n_flu == 1) { cout << " single-fluid mode is on." << endl; }
			else            { cout << " " << n_flu << " fluids found." << endl; }
		}
		
		// --> get grid size, save as 'fulldimsZYX'
		id_dset  = H5Dopen1( id_file, (qstruct[0].path).c_str() );
		id_space = H5Dget_space( id_dset );
		idummy   = H5Sget_simple_extent_dims( id_space, fulldimsZYX, NULL );
		status   = H5Dclose( id_dset );
		status   = H5Fclose( id_file );
		
		int Nfc[3];   // current (full) grid size
		for (short comp = 0; comp < 3; comp++) {
			Nfc[comp] = int(fulldimsZYX[comp]);
			if (it == 0)  { Nfgrid[comp] = Nfc[comp]; }
			else { same_grid = same_grid and (Nfc[comp] == Nfgrid[comp]); }
		}
		
		printf( " fulldimsXYZ = [ %d %d %d ], rim = %d, vars = { ",
		        Nfc[2], Nfc[1], Nfc[0], rim[0] );
		for (int ient = 0; ient < int(qstruct.size()); ++ient)
			cout << qstruct[ient].name << " ";
		cout << "}" << endl;
		
	   // --> sort out vector components into separate list
		vstruct.clear();
		for (int ient = 0; ient < int(qstruct.size()); ++ient) {
			string cu_para = (qstruct[ient]).para;  // current PV name
			string last_2c = cu_para.substr(cu_para.size()-2, 2 );
			if (last_2c == "_x") {     // --> look for (y,z) components
				string vecname = remove_last2( cu_para );
				if (listpos( vecname+"_y", qstruct ) > -1 and  // all comp.s
					 listpos( vecname+"_z", qstruct ) > -1) {  //   present?
					quant newvec;
					newvec.name = remove_last2( qstruct[ient].name );
					newvec.path = remove_last2( qstruct[ient].path );
					newvec.para = vecname;
					vstruct.push_back( newvec );   // => add to vector list
					// cout << newvec.name << endl;
					// --> replace subscripts {x,y,z} for non-Cartesian coords.
					for (unsigned int iq = 0; iq < qstruct.size(); iq++)
						for (int id = 0; id < 3; id++)
							if (qstruct[iq].name == vecname + "_" + subs_Car[id]) {
								int len = (qstruct[iq].para).size();
								qstruct[iq].para = (qstruct[iq].para).substr( 0, len-1 ) + subs[id];
							}
				}
			}
			qstruct0 = qstruct; // save frame #0 as reference
		}
		same_vars = same_vars and (qstruct == qstruct0);
	}
	cout << " Reading done." << endl;
	
	// --> report consistency check
	string res[2];
	res[true] = "--OK--"; res[false] = "-FAIL-";
	for (int it = 1; it < Nf; it++)
		same_rim = same_rim and (rim[it] == rim[0]);
	cout << " [consistency] grid size: " << res[same_grid] << endl;
	cout << "                boundary: " << res[same_rim ] << endl;
	cout << "               variables: " << res[same_vars] << endl;	
	if (same_grid and same_rim and same_vars == false) {
		abort( "inconsistent input data" );
	}
	int n_vec = vstruct.size();
	if (n_vec == 0)
		cout << " No vectors found." << endl;
	else {
		cout << " Set of vectors found: {";
		for (int iv = 0; iv < n_vec; ++iv)
			cout << " " << vstruct[iv].name << " ";
		cout << "}" << endl;
	}
	
	// --> staggered B needs at least one ghost layer!
	bool flag_mag = false;       // flag if B present
	int   pos_mag = -1;          // position of B in list
	for (int iv = 0; iv < n_vec; ++iv)
		if (vstruct[iv].name == "B") {
			flag_mag = true;
			pos_mag = iv;
		}
	if (flag_mag and rim == 0) {   // -> remove B from vector list
		cout << "Warning: datafile lacks boundary,"
			  << " B vector field access disabled." << endl;
		vstruct.erase( vstruct.begin() + pos_mag ); 
		--n_vec;
	}
	
	const string tss = string( "        <!-- ====== time slice ")
		+               string( "separator ====== -->\n        " );
	FILE *xmf = 0;   // --> create and open Xdmf output file
	xmf = fopen( x_dirfile.c_str(), "w" );
	
	paste_block( xmf, "beg" );  // write Xdmf head(er) part
	paste_line( xmf, "<!ENTITY grid_file \"" + gridfile + "\">");
	for (int it = 0; it < Nf; it++)
		paste_line( xmf, "<!ENTITY data_file" + itime[it] + " \""
		            + hdf5file[it] + "\">" );
	fprintf( xmf, "<!ENTITY fullDimsZYX \"%d %d %d\">\n",
	         Nfgrid[0], Nfgrid[1], Nfgrid[2] );
	fprintf (xmf, "<!ENTITY cellDimsZYX \"%d %d %d\">\n",
	         Nfgrid[0]-2*rim[0], Nfgrid[1]-2*rim[0], Nfgrid[2]-2*rim[0] );
	fprintf( xmf, "<!ENTITY nodeDimsZYX \"%d %d %d\">\n",
	         Nfgrid[0]-2*rim[0]+1, Nfgrid[1]-2*rim[0]+1, Nfgrid[2]-2*rim[0]+1 );
	fprintf( xmf, "<!ENTITY ghost_cells \"%d\">\n", rim[0]   );
	fprintf( xmf, "<!ENTITY ghostminus1 \"%d\">\n", rim[0]-1 );
	paste_line( xmf, "<!ENTITY coordName_1 \"" + subs[0] + "\">" );
	paste_line( xmf, "<!ENTITY coordName_2 \"" + subs[1] + "\">" );
	paste_line( xmf, "<!ENTITY coordName_3 \"" + subs[2] + "\">" );
	paste_line( xmf, "]>");
	
	paste_block (xmf, "dTA");  // def. of transfo matrix + auxilliary terms 
	
	// --> write definition for vector components
	for (int it = 0; it < Nf; it++) {
		string s_sep = "{define vector components, (" + itime[it] + ")=itime}";
		paste_line( xmf, "  <!-- BEGIN " + s_sep + " -->\n  " );
		for (int iv = 0; iv < n_vec; ++iv)
			for (int id = 0; id < 3; ++id) {
				string compname = vstruct[iv].para + "_" + subs_Car[id];
				string pathname = vstruct[iv].path + "_" + subs_Car[id];
				string label = "gVe";              // generic vector
				if (vstruct[iv].name == "B")
					label = "dB" + subs_Car[id];        // B vector
				paste_block( xmf, label, compname, itime[it], pathname );
			}
		paste_line( xmf, ("  <!-- END " + s_sep + " -->\n  ") );
	}
	paste_block( xmf, "do1" ); // start domain declaration, begin Cgrid
	
	// --> C grid: print vector attributes, looping over time
	for (int it = 0; it < Nf; it++) {
		if (it > 0)  { paste_line( xmf, tss ); }
		paste_block( xmf, "Cg1", time[it] );
		for (int iv = 0; iv < n_vec; ++iv)
			paste_block( xmf, "AtV", vstruct[iv].name, itime[it] );
		paste_block( xmf, "Atu" );  // attributes for unit vectors
	}
	paste_block( xmf, "do2" );     // end Cgrid, begin Ngrid
	
	// --> N grid: print scalar attributes, looping over time
	for (int it = 0; it < Nf; it++) {
		if (it > 0)  { paste_line( xmf, tss ); }
		paste_block( xmf, "Ng1", time[it] );
		for (int is = 0; is < int(qstruct.size()); ++is) {
			paste_block( xmf, "AtS", qstruct[is].para,
							 itime[it],  qstruct[is].path );
		}
		paste_line( xmf, "        </Grid>\n        " );
	}
	paste_block( xmf, "end" );   // closing tags
	fclose( xmf );
	cout << " Light data written to '" << x_dirfile << "'." << endl;
	
	delete[] itime;
	delete[] time;
	delete[] hdf5file;
	delete[] rim;
	return 0;
}
