#include "data.H"

#include <iomanip>
#include <sstream>

#include "timewrapper.H"
#include "Hdf5File_cbase.H"
#include "CException.H"

using namespace std;

Pot::Pot() {
	//FieldType = "undetermined";
}

Pot::Pot(const int *mx)
{
	// FieldType = "undetermined";
	resize(mx);
}

Pot::Pot(const NumMatrix<DATA_TYPE, DIM>& _matr)
	: NumMatrix<DATA_TYPE, DIM>(_matr) {
	// FieldType = "undetermined";
}

void Pot::resize(const int *mx)
{
	int l[DIM], h[DIM];
	for (int d = 0; d < DIM; d++) {
		l[d] = -B;
		h[d] = mx[d]+B;
	}
	NumMatrix<DATA_TYPE,DIM>::resize(l,h);
}


void Pot::resize(const int *mx, int rim)
{
	int l[DIM], h[DIM];
	for (int d = 0; d < DIM; d++) {
		l[d] = -rim;
		h[d] = mx[d]+rim;
	}
	NumMatrix<DATA_TYPE,DIM>::resize(l,h);
}


void Pot::resize_ll(const int *ib, const int *ie)
{
	const int* ib_old = NumMatrix<DATA_TYPE,DIM>::getLow();
	const int* ie_old = NumMatrix<DATA_TYPE,DIM>::getHigh();

	int ib_loc[DIM], ie_loc[DIM];

	for(int i=0; i<DIM; ++i) {
		ib_loc[i]=std::max(ib_old[i], ib[i]);
		ie_loc[i]=std::min(ie_old[i], ie[i]);
	}

	NumMatrix<DATA_TYPE, DIM> matr_back(*this);

	for (int i = 0; i < size; i++)
		matr_back[i] = matr[i];

	NumMatrix<DATA_TYPE,DIM>::resize(ib,ie);

	if(DIM == 3) {
		for (int k=ib_loc[2]; k<=ie_loc[2]; ++k) {
			for (int j=ib_loc[1]; j<=ie_loc[1]; ++j) {
				for (int i=ib_loc[0]; i<=ie_loc[0]; ++i) {
					NumMatrix<DATA_TYPE,DIM>::operator()(i,j,k) = matr_back(i,j,k);
				}
			}
		}
	} else if (DIM == 2) {
		for (int j=ib_loc[1]; j<=ie_loc[1]; ++j) {
			for (int i=ib_loc[0]; i<=ie_loc[0]; ++i) {
				NumMatrix<DATA_TYPE,DIM>::operator()(i,j) = matr_back(i,j);
			}
		}
	} else if (DIM == 1) {
		for (int i=ib_loc[0]; i<=ie_loc[0]; ++i) {
			NumMatrix<DATA_TYPE,DIM>::operator()(i) = matr_back(i);
		}
	} else {
		cerr << " resize_ll not implemented for DIM > " << DIM << endl;
		exit(-3);
	}


}


void Pot::set_max(const double &val)
{
	for(int i=0; i < this->size; i++) {
		if(this->matr[i] > val) {
			this->matr[i] = val;
		}
	}
	return;
}


double Pot::get_max()
{
	double maximum = static_cast<double>(this->matr[0]);
	for(int i=0; i<this->size; i++) {
		maximum = std::max(maximum, static_cast<double>(this->matr[i]));
	}
	return maximum;
}


Pot& Pot::operator=(const double& val)
{
	for (int i = 0; i < size; i++)
		matr[i] = val;
	return *this;
}

void Pot::set_min(const double &val)
{
	for(int i=0; i < this->size; i++) {
		if(this->matr[i] < val) {
			this->matr[i] = val;
		}
	}
	return;
}


double Pot::get_min()
{
	double minimum = static_cast<double>(this->matr[0]);
	for(int i=0; i<this->size; i++) {
		minimum = std::min(minimum, static_cast<double>(this->matr[i]));
	}
	return minimum;
}

Data::Data() : cflSYCL (celerity::buffer<double,1>(celerity::range{1})), nomSYCL (CelerityBuffer<nom_t, 3>(celerity::range<3>(mx[0]+6 +1, mx[1]+6+1, mx[2]+6+1)))
{

	// Start runtime timer
	gettimeofday(&tick_start, 0);

#if(CRSWITCH_DUAL_ENERGY == CRONOS_ON)
	fluid.set_dualEnergy(AUX_ENERGY);
#else
	fluid.unset_dualEnergy();
#endif
	fluid.setup(FLUID_TYPE, ENERGETICS, N_ADD, N_SUBS, N_OMINT_USER);
	int N_OM = fluid.get_N_OM();
	int n_omInt = fluid.get_N_OMINT();
	int n_omUser = fluid.get_N_OM_USER();
//	N_OM = 8;

	this->rim = B;
	constexpr int numElements = 12;
	om = new Pot[numElements];
	//om = new Pot[N_OM+N_P]; //necessary for MHD

	for (int i = 0; i < numElements; i++) {
		omSYCL.push_back(CelerityBuffer<double, 3>(Range<3>(mx[0]+6 +1, mx[1]+6+1, mx[2]+6+1)));
		omSYCL_out.push_back(CelerityBuffer<double, 3>(Range<3>(mx[0]+2 +1, mx[1]+2+1, mx[2]+2+1)));
		omSYCL_out_flt.push_back(CelerityBuffer<float, 3>(Range<3>(mx[0]+1, mx[1]+1, mx[2]+1)));
	}

	for (int i = 0; i < 1; i++) {
		pThermSYCL.push_back(CelerityBuffer<double, 3>(Range<3>(mx[0]+6+1, mx[1]+6+1, mx[2]+6+1)));
		carbuncleFlagSYCL.push_back(CelerityBuffer<int, 3>(Range<3>(mx[0]+6+1, mx[1]+6+1, mx[2]+6+1)));
	}

	nom = new NumMatrix<double,3> [n_omInt];

	for (int q = 0; q < n_omInt; ++q) {
		nom[q].resize(Index::set(0,0,0), Index::set(mx[0],mx[1],mx[2]));
	}

//	exit(3);
//	om.clear();
//	om.reserve(N_OM+N_P);
	for(int q=0; q<N_OM+N_P; ++q) {
//		cout << " kuh " << q << endl;
		// Make om fields
//		om.push_back(Pot(mx));
//		cout << " next " << endl;
		om[q].clear();
		// om[q].set_VariableName("om");
	}

#if (OMS_USER == TRUE)
	int n_omIntUser = fluid.get_N_OMINT_USER();

	om_user = new Pot[n_omUser];
//	nom_user = new NumMatrix<double,3> [n_omUser];
	nom_user = new NumMatrix<double,3> [n_omIntUser];

//	for (int q = 0; q < n_omUser; ++q) {
	for (int q = 0; q < n_omIntUser; ++q) {
		nom_user[q].resize(Index::set(0,0,0), Index::set(mx[0],mx[1],mx[2]));
	}
	for(int q=0; q<n_omUser; ++q) {
		om_user[q].clear();
		// om_user[q].set_VariableName("om_user");
	}
#else
	nom_user = NULL;
#endif

	
	t_end = value((char*)"t_end");
	dt = value((char*)"dt");
	cfl = 0.;
	tstep = 0; // Timestep

	for (int ll=0; ll<N_OM+N_P; ++ll) {
		string name = "om";
		char cll[255];
		sprintf(cll,"%2.2d",ll);
		name += cll;
    
//		om.push_back(new Pot(mx));
		om[ll].resize(mx);
		om[ll].rename(name);
		om[ll].clear();
	}

#if (OMS_USER == TRUE)
    for (int q=0; q<n_omUser; ++q) {
		string name = "om_user";
		char cll[255];
		sprintf(cll,"%2.2d",q);
		name += cll;
    
		om_user[q].resize(mx);
		om_user[q].rename(name);
		om_user[q].clear();
    }
#endif

#if (GEOM == SPHERICAL)
    massFlux.resize(Index::set(0), Index::set(mx[0]+1));

    // Fill with corresponding area
    radialArea = 0.;
    for(int iPhi=0; iPhi<=mx[2]; ++iPhi) {
    	for(int iTheta=0; iTheta<=mx[1]; ++iTheta) {
#if (NON_LINEAR_GRID == CRONOS_OFF)
    		double f_geom_x = h2(0.5,iTheta,iPhi)/h1(0.5,iTheta,iPhi);
#else
    		double f_geom_x = h2(0,iTheta,iPhi,1,0,0)/h1(0,iTheta,iPhi,1,0,0);
#endif
    		radialArea += f_geom_x*getCen_dx(1,iTheta)*getCen_dx(2,iPhi);
    	}
    }

#endif

	// Check whether thermal pressure is needed
    if(value_exists("use_carbuncleFlag")) {
    	use_carbuncleFlag = true;
    } else {
    	use_carbuncleFlag = false;
    }
	if(use_carbuncleFlag) {
		storePressure = true;
		// Resize pressure matrix
		pTherm.resize(Index::set(-B,-B,-B), Index::set(mx[0]+B,mx[1]+B,mx[2]+B));
		carbuncleFlag.resize(Index::set(-B,-B,-B), Index::set(mx[0]+B,mx[1]+B,mx[2]+B));

	} else {
		storePressure = false;
	}

    mag = value((char*)"mag");
    time = 0.;

}

Data::~Data() {
	delete [] om;
	delete [] nom;
	omSYCL.clear();
	omSYCL_out.clear();
	omSYCL_out_flt.clear();
	pThermSYCL.clear();
	carbuncleFlagSYCL.clear();

#if (OMS_USER == TRUE)
	delete [] om_user;
	delete [] nom_user;
#endif

	// Timer for code finalisation
	gettimeofday(&tock_end, 0);

	double full_time_sec = ((tock_end.tv_sec + tock_end.tv_usec/1.e6) -
			(tick_start.tv_sec + tick_start.tv_usec/1.e6));
	double time_per_step = full_time_sec/(1.*tstep);
	if(rank==0) {
		cout << "======================================================" << endl;
		cout << " Total run time               " << git_humanReadable(full_time_sec) << endl;
		cout << " Average time step duration   " << git_humanReadable(time_per_step) << endl;
		cout << "======================================================" << endl;
	}
}



std::string Data::git_humanReadable(double t_sec) {
	int t_int_sec = static_cast<int>(t_sec);
	int full_time_min = t_int_sec/60;
	int full_time_hours = full_time_min/60;
	int full_time_days = full_time_hours/24;
	double full_time_sec = t_sec - 60*full_time_min;
	full_time_min -= 60*full_time_hours;
	full_time_hours -= 24*full_time_days;

	std::ostringstream t_string_d;
	t_string_d << full_time_days << "d:" << full_time_hours << "h:" << full_time_min << "m:" << full_time_sec << "s";
	std::ostringstream t_string_h;
	t_string_h << full_time_hours << "h:" << full_time_min << "m:" << full_time_sec << "s";
	std::ostringstream t_string_m;
	t_string_m << full_time_min << "m:" << full_time_sec << "s";
	std::ostringstream t_string_s;
	t_string_s << full_time_sec << "s";
	std::string t_string = t_string_s.str();
	if(full_time_days > 0) {
		t_string = t_string_d.str();
	} else if (full_time_hours > 0) {
		t_string = t_string_h.str();
	} else if (full_time_min > 0) {
		t_string = t_string_m.str();
	}
	return t_string;
}


double Data::computeInt(int q)
{
	double sum = 0.;
	double dV(0.);
#if (NON_LINEAR_GRID == CRONOS_OFF)
	dV = get_CellVolume(0,0,0);
#endif
	for(int k=0; k<numCellsEff[2]; ++k) {
		for(int j=0; j<numCellsEff[1]; ++j) {
			for(int i=0; i<numCellsEff[0]; ++i) {
#if (GEOM != CARTESIAN || NON_LINEAR_GRID == CRONOS_ON)
				dV = get_CellVolume(i,j,k);
#endif
				sum += om[q](i,j,k)*dV;
			}
		}
	}

	return sum;
}


double Data::computeRMS(int q)
{
	double sum(0.);
	double dV(dx[0]*dx[1]*dx[2]);
	double Vol(0.);
	for(int k=0; k<numCellsEff[2]; ++k) {
		for(int j=0; j<numCellsEff[1]; ++j) {
			for(int i=0; i<numCellsEff[0]; ++i) {
				sum += sqr(om[q](i,j,k))*dV;
				Vol += dV;
			}
		}
	}

	sum /= Vol;
	return sqrt(sum);
}


NumMatrix<double,3> Data::computeAbs(int qmin, int qmax) const
{
	int ibeg[3] = {om[qmin].getLow(0), om[qmin].getLow(1), om[qmin].getLow(2)};
	int iend[3] = {om[qmin].getHigh(0), om[qmin].getHigh(1), om[qmin].getHigh(2)};
  
	NumMatrix<double,3> Abs;
	Abs.resize(Index::set(ibeg[0],ibeg[1],ibeg[2]),
	           Index::set(iend[0],iend[1],iend[2]));
	Abs.clear();

	for(int q=qmin; q<=qmax; ++q) {
		for(int k=ibeg[2]; k<=iend[2]; ++k) {
			for(int j=ibeg[1]; j<=iend[1]; ++j) {
				for(int i=ibeg[0]; i<=iend[0]; ++i) {
					Abs(i,j,k) += sqr(om[q](i,j,k));
				}
			}
		}
	}

	for(int k=ibeg[2]; k<=iend[2]; ++k) {
		for(int j=ibeg[1]; j<=iend[1]; ++j) {
			for(int i=ibeg[0]; i<=iend[0]; ++i) {
				Abs(i,j,k) = sqrt(Abs(i,j,k));
			}
		}
	}

	return Abs;
}


double Data::getMin(int q)
{
	int ibeg[3] = {om[q].getLow(0), om[q].getLow(1), om[q].getLow(2)};
	int iend[3] = {om[q].getHigh(0), om[q].getHigh(1), om[q].getHigh(2)};
	double minimum(abs(om[q](ibeg[0],ibeg[1],ibeg[2])));

	for(int k=ibeg[2]; k<=iend[2]; ++k) {
		for(int j=ibeg[1]; j<=iend[1]; ++j) {
			for(int i=ibeg[0]; i<=iend[0]; ++i) {

				minimum = std::min(minimum, om[q](i,j,k));

			}
		}
	}

	return minimum;
}


double Data::getMin(int pos[3], int q)
{
	int ibeg[3] = {om[q].getLow(0), om[q].getLow(1), om[q].getLow(2)};
	int iend[3] = {om[q].getHigh(0), om[q].getHigh(1), om[q].getHigh(2)};
	pos[0] = ibeg[0];
	pos[1] = ibeg[1];
	pos[2] = ibeg[2];
	double minimum(abs(om[q](ibeg[0],ibeg[1],ibeg[2])));

	for(int k=ibeg[2]; k<=iend[2]; ++k) {
		for(int j=ibeg[1]; j<=iend[1]; ++j) {
			for(int i=ibeg[0]; i<=iend[0]; ++i) {

				if(om[q](i,j,k) < minimum) {
					minimum = om[q](i,j,k);
					pos[0] = i;
					pos[1] = j;
					pos[2] = k;
				}
			}
		}
	}

	return minimum;
}


double Data::getMax(int q)
{
	int ibeg[3] = {om[q].getLow(0), om[q].getLow(1), om[q].getLow(2)};
	int iend[3] = {om[q].getHigh(0), om[q].getHigh(1), om[q].getHigh(2)};
	double maximum(abs(om[q](ibeg[0],ibeg[1],ibeg[2])));

	for(int k=ibeg[2]; k<=iend[2]; ++k) {
		for(int j=ibeg[1]; j<=iend[1]; ++j) {
			for(int i=ibeg[0]; i<=iend[0]; ++i) {

				maximum = std::max(maximum, om[q](i,j,k));

			}
		}
	}

	return maximum;
}



double Data::getMax(int pos[3], int q)
{
	int ibeg[3] = {om[q].getLow(0), om[q].getLow(1), om[q].getLow(2)};
	int iend[3] = {om[q].getHigh(0), om[q].getHigh(1), om[q].getHigh(2)};
	pos[0] = ibeg[0];
	pos[1] = ibeg[1];
	pos[2] = ibeg[2];
	double maximum(abs(om[q](ibeg[0],ibeg[1],ibeg[2])));

	for(int k=ibeg[2]; k<=iend[2]; ++k) {
		for(int j=ibeg[1]; j<=iend[1]; ++j) {
			for(int i=ibeg[0]; i<=iend[0]; ++i) {

				if(om[q](i,j,k) > maximum) {
					maximum = om[q](i,j,k);
					pos[0] = i;
					pos[1] = j;
					pos[2] = k;
				}
			}
		}
	}

	return maximum;
}






NumMatrix<float,3> Data::float_data(int q, int R, bool generic)
{

	NumMatrix<float,3> datafloat(Index::set(-R,-R,-R),
	                             Index::set(mx[0]+R,mx[1]+R,mx[2]+R));
  
#if (OMS_USER == TRUE)
	if(generic) {
#endif
		for (int k = -R; k <= mx[2]+R; k++) {
			for (int j = -R; j <= mx[1]+R; j++) {
				for (int i = -R; i <= mx[0]+R; i++) {
					datafloat(i,j,k) = (q < 0) ? 0. : static_cast<float>(om[q](i,j,k));
				}
			}
		}
#if (OMS_USER == TRUE)
	} else {
		for (int k = -R; k <= mx[2]+R; k++) {
			for (int j = -R; j <= mx[1]+R; j++) {
				for (int i = -R; i <= mx[0]+R; i++) {
					datafloat(i,j,k) = (q < 0) ? 0. : static_cast<float>(om_user[q](i,j,k));
				}
			}
		}
	}
#endif
	return datafloat;

}



void Data::floatom_out(int R, Hdf5Stream & h5file, int om_lo, int om_hi)
{
	double xbvar[3];
	xbvar[0] = xb[0]-dx[0]*R;
	xbvar[1] = xb[1]-dx[1]*R;
	xbvar[2] = xb[2]-dx[2]*R;
	//   string omName;
	NumMatrix<float,3> data(Index::set(-R,-R,-R),
	                        Index::set(mx[0]+R,mx[1]+R,mx[2]+R));

	for (int q = om_lo; q <= om_hi; q++) {
		for (int k = -R; k <= mx[2]+R; k++) {
			for (int j = -R; j <= mx[1]+R; j++) {
				for (int i = -R; i <= mx[0]+R; i++) {
					data(i,j,k) = (q < 0) ? 0. : (float)om[q](i,j,k);
				}
			}
		}
		if(!h5file.Write3DMatrix(om[q].getName(), data, xbvar, dx)){
			cerr << " Writing not successful for om: " << om[q].getName() << endl;
			exit(-15);
		}
	}
}



void Data::CheckNan(int pos)
{
	char message[255];
	sprintf(message,"Check: %2.2d om is NAN",pos);

	int found_NAN = 0;
	int qNAN(-10), iNAN(-10), jNAN(-10), kNAN(-10);

	for (int q = 0; q < N_OMINT; ++q){

		int lo[3], up[3];
		lo[0] = om[q].getLow(0);
		lo[1] = om[q].getLow(1);
		lo[2] = om[q].getLow(2);
    
		up[0] = om[q].getHigh(0);
		up[1] = om[q].getHigh(1);
		up[2] = om[q].getHigh(2);
      
		for (int k = lo[2]; k <= up[2]; k++){
			for (int j = lo[1]; j <= up[1]; j++){
				for (int i = lo[0]; i <= up[0]; i++){
					if(std::isnan(om[q](i,j,k)) && !found_NAN) {
						found_NAN = 1;
						qNAN = q;
						iNAN = i;
						jNAN = j;
						kNAN = k;
					}
				}
			}
		}
	}

	if(found_NAN) {
		throw CException(message,pos,qNAN,iNAN,jNAN,kNAN);
	}
}



void Data::set_fieldIds() {
	// std::map<NumMatrix<double,3> *,int> id_numbers;

	for(int q=0; q<N_OM+N_P; ++q) {
		string name = om[q].getName();
		om_ids.enter_fieldId(name);
	}

#if (OMS_USER == TRUE)
	for(int q=0; q<N_OM_USER; ++q) {

		om_ids.enter_fieldId(om_user[q].getName());

	}
#endif

	// Print map to screen
	if(rank == 0) {
		cout << "------- Field ids-------------------------------------";
		cout << endl;
		om_ids.print_map();
		// cout << "------------------------------------------------------" << endl;
		// cout << endl;
	}

}

bool Data::is_userField(NumMatrix<double,3> &omField) {
	
	if(om_ids.get_fieldId(omField) >= om_ids.max_generic && 
	   om_ids.get_fieldId(omField) < om_ids.fieldIds.size()) {
		return true;
	} else {
		return false;
	}

}

void Data::fetch_cfl(Queue &queue) {

	double cfl_lin(0.);

	// queue.slow_full_sync();
	
	queue.submit(celerity::allow_by_ref, [=, &cfl_lin](celerity::handler& cgh) {
		celerity::accessor cflSYCL_acc{this->cflSYCL, cgh, celerity::access::all{}, celerity::read_only_host_task};
		cgh.host_task(celerity::experimental::collective, 
				[=, &cfl_lin](celerity::experimental::collective_partition /*part*/){	
			std::this_thread::sleep_for(std::chrono::milliseconds(100));
			// MPI_Comm comm = part.get_collective_mpi_comm();
        	// MPI_Barrier(comm);
			if (cflSYCL_acc[0] > 0.) {
				cfl_lin = cflSYCL_acc[0];
			}
		});
	});

	queue.slow_full_sync();
	
	this->cfl = std::max(cfl_lin, this->cfl);

}



id_handler::id_handler() {
	max_generic = N_OM+N_P;
}


void id_handler::enter_fieldId(string FieldName) {
	// Add entry to fieldIds map - only if not already used:
	size_t used = fieldIds.count(FieldName);

	if(used > 0) {
		
		typedef map<std::string, unsigned int> stringType;
		stringType::iterator stringIter;
		stringIter = fieldIds.find(FieldName);

		cerr << FieldName << " already in uses as ";
		cerr << stringIter->second;
		cerr << endl;
		exit(-3);

	} else {

		unsigned int id_number = (unsigned int)fieldIds.size();
		fieldIds[FieldName] = id_number;

	}
	

}



bool id_handler::is_generic(NumMatrix<double,3> &omField) {
	
	if(get_fieldId(omField) < max_generic) {
		return true;
	} else {
		return false;
	}

}



unsigned int id_handler::get_userFieldId(NumMatrix<double,3> &omField) {
	return get_fieldId(omField)-(max_generic);
}

unsigned int id_handler::get_fieldId(NumMatrix<double,3> &omField) {

	// get unique id number for field
	typedef map<std::string, unsigned int> id_map;
	id_map::iterator idIter;
	idIter = fieldIds.find(omField.getName());

	return idIter->second;

}


void id_handler::print_map() {

	typedef map<std::string, unsigned int>::iterator idIter;

	for(idIter id=fieldIds.begin(); id!=fieldIds.end(); ++id) {
		// cout << id->second << " -> " << id->first << endl;
		cout << setiosflags(ios::left);
		cout << setw(9) << id->first << " -> ";
		cout << resetiosflags(ios::left);
		cout << setiosflags(ios::right )<< setw(2) << id->second << endl;
		cout << resetiosflags(ios::right);
	}
	// cout << resetiosflags(ios::left);


}