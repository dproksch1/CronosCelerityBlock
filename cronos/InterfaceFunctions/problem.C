#include "problem.H"
#include "string.h"
#include <stdlib.h>

ProblemType::ProblemType(const Data &gdata) {

	//	if (strcmp(svalue((char*)"Fluid_Type"),"MHD") == 0) {
#if (FLUID_TYPE == CRONOS_MHD)
		mag = true;
#else
		//	} else {
		mag = false;
		//	}
#endif
	//	this->mag = value((char*)"mag");
	this->gamma = value((char*)"Adiabatic_exponent");
	this->cs2 = sqr(value((char*)"Isothermal_Soundspeed"));
	this->rho0 = value((char*)"Initial_density");
	this->info = true;
	this->asciiOut = true;
	this->type = static_cast<int>(value((char*)"type"));

	for(int i=0; i<9; ++i) {
		infodata[i] = false;
	}

#ifdef PHYSDISS
	this->nu0 = value((char*)"Viscosity");
	REAL Pm = value((char*)"PrandtlNumber");
	this->eta0  = nu0/Pm;
	REAL Re = sqrt(cs2)/nu0;
	REAL ReM = sqrt(cs2)/eta0;

	cout << endl << " kinetic Reynolds number:   " << Re  << endl;
	cout << " magnetic Reynolds number:  " << ReM << endl;
#endif


	// Get q-Values from fluid class:
#if(FLUID_TYPE==CRONOS_MULTIFLUID)
	n_omInt = gdata.fluids->get_N_OMINT();
	n_omIntAll = gdata.fluids->get_N_OMINT_ALL();

	// In multifluid case use first fluid:
	q_rho = gdata.fluids->fluids[0].get_q_rho_global();
	q_sx= gdata.fluids->fluids[0].get_q_sx_global();
	q_sy= gdata.fluids->fluids[0].get_q_sy_global();
	q_sz= gdata.fluids->fluids[0].get_q_sz_global();
	q_Bx= gdata.fluids->get_q_Bx();
	q_By= gdata.fluids->get_q_By();
	q_Bz= gdata.fluids->get_q_Bz();
	q_Eges= gdata.fluids->fluids[0].get_q_Eges_global();
	q_Eadd= gdata.fluids->fluids[0].get_q_Eadd_global();

	q_Ax = n_omIntAll+N_ADD;
	q_Ay = n_omIntAll+N_ADD+1;
	q_Az = n_omIntAll+N_ADD+2;

#else
	n_omInt = gdata.fluid.get_N_OMINT();
	n_omIntAll = gdata.fluid.get_N_OMINT_ALL();

	q_rho = gdata.fluid.get_q_rho();
	q_sx = gdata.fluid.get_q_sx();
	q_sy = gdata.fluid.get_q_sy();
	q_sz = gdata.fluid.get_q_sz();
	q_Bx = gdata.fluid.get_q_Bx();
	q_By = gdata.fluid.get_q_By();
	q_Bz = gdata.fluid.get_q_Bz();
	q_Eges = gdata.fluid.get_q_Eges();
	q_Eadd = gdata.fluid.get_q_Eadd();

	q_Ax = n_omInt+N_ADD;
	q_Ay = n_omInt+N_ADD+1;
	q_Az = n_omInt+N_ADD+2;
#endif


	// default values for mass fractions and mean-molecular weight
	massFraction_X = 1.;
	massFraction_Y = 0.;
	double massFraction_Z = 1.-massFraction_X - massFraction_Y;

	Quantity mu = CRONOS_CONSTANTS::AtomicMassUnit/(2.*massFraction_X + 0.75*massFraction_Y + 0.5*massFraction_Z);
	meanParticleMass = mu;


	// Indicate that pointer TrafoNorm is not yet initialised
	TrafoNorm = NULL;

}


int pow(int value, unsigned int exponent) {
	int result(1);
	for(unsigned int i=0; i<exponent; ++i) {
		result *= value;
	}
	return result;
}


void ProblemType::set_InfoData(Data &gdata) {

	int debug = static_cast<int>(value((char*)"debug"));
	if(debug > 0 && gdata.rank == 0) {
		cout << "======================================================" << endl;
		cout << " Writing info on... " << endl;
	}

	// Compute Output information:
	for(int i=6; i>=0; --i) {
		int val = debug - pow(2, static_cast<unsigned int>(i));
		if(val > -1) {
			debug -= pow(2,static_cast<unsigned int>(i));
			infodata[i] = true;
		} else {
			infodata[i] = false;
		}
	}
	if(infodata[0]) {
		if(gdata.rank == 0) {
			cout << "  ... Runge Kutta Stepping " << endl;
		}
	}
	if(infodata[1]) {
		if(gdata.rank == 0) {
			cout << "  ... value of cfl number " << endl;
		}
	}
	if(infodata[2]) {
		if(gdata.rank == 0) {
			cout << "  ... I/O of data " << endl;
		}
	}
	if(infodata[3]) {
		if(gdata.rank == 0) {
			cout << "  ... Energetics " << endl;
		}
	}
	if(infodata[4]) {
		if(gdata.rank == 0) {
			cout << "  ... divergence of B " << endl;
		}
	}
	if(infodata[5]) {
		if(gdata.rank == 0) {
			cout << "  ... Timing " << endl;
		}
	}
#if(CRSWITCH_DUAL_ENERGY == CRONOS_ON)
	if(infodata[6]) {
		if(gdata.rank == 0) {
			cout << "  ... Energy correction " << endl;
		}
	}
#endif
	if(debug > 0 && gdata.rank == 0) {
		cout << "======================================================" << endl;
	}

}

bool ProblemType::checkout(int num) {
	if(num > 7 || num < 0) {
		cerr << " Error: no such value " << endl;
		return false;
	} else {
		return infodata[num];
	}
}

void ProblemType::set_Info(bool info) {
	this->info = info;
}

void ProblemType::set_AsciiOut(bool asciiOut) {
	this->asciiOut = asciiOut;
}

bool ProblemType::get_Info() {
	return info;
}

bool ProblemType::get_AsciiOut() {
	return asciiOut;
}

int ProblemType::get_Type()
{
	return type;
}


string ProblemType::get_Name()
{
	return name;
}


void ProblemType::init_fields(Data &gdata)
{
	int ibeg[3] = {gdata.om[0].getLow(0), gdata.om[0].getLow(1),
	               gdata.om[0].getLow(2)};
	int iend[3] = {gdata.om[0].getHigh(0), gdata.om[0].getHigh(1),
	               gdata.om[0].getHigh(2)};
	init_fields(gdata, ibeg, iend);
}


void ProblemType::init_fields(Data &gdata, int q) {

	int ibeg[3] = {gdata.om[q].getLow(0), gdata.om[q].getLow(1),
	               gdata.om[q].getLow(2)};
	int iend[3] = {gdata.om[q].getHigh(0), gdata.om[q].getHigh(1),
	               gdata.om[q].getHigh(2)};

	// Distinction of different kinds of variables
	bool centered[3] = { true, true, true };
	
	string VarName = gdata.om[q].getName();

	if(VarName == "B_x") {
		centered[0] = false;
	}
	if(VarName == "B_y") {
		centered[1] = false;
	}
	if(VarName == "B_z") {
		centered[2] = false;
	}

	if(VarName == "A_x") {
		centered[1] = false;
		centered[2] = false;
	}
	if(VarName == "A_y") {
		centered[0] = false;
		centered[2] = false;
	}
	if(VarName == "A_z") {
		centered[0] = false;
		centered[1] = false;
	}
	
	// Careful -> the current gridding ONLY works for the vector
	// potential - not for the magnetic field!


	for(int iz=ibeg[2]; iz<=iend[2]; ++iz) {
		for(int iy=ibeg[1]; iy<=iend[1]; ++iy) {
			for(int ix=ibeg[0]; ix<=iend[0]; ++ix) {

				REAL xPos(0.), yPos(0.), zPos(0.);
				if(centered[0]) {
					xPos = gdata.getCen_x(ix);
				} else {
					xPos = gdata.getEdgL_x(ix);
				}
				if(centered[1]) {
					yPos = gdata.getCen_y(iy);
				} else {
					yPos = gdata.getEdgL_y(iy);
				}

				if(centered[2]) {
					zPos = gdata.getCen_z(iz);
				} else {
					zPos = gdata.getEdgL_z(iz);
				}

			}
		}
	}

}



int ProblemType::get_q(NumMatrix<REAL,3> &omField, string fluidName) {

	typedef map<std::string, int> stringType;
#if (USE_COROTATION == CRONOS_ON)
	string omName = omField.getName();
	if(omName == "v_x_Corot") omName = "v_x";
	if(omName == "v_y_Corot") omName = "v_y";
	if(omName == "v_z_Corot") omName = "v_z";
	string fieldName = fluidName + "/" + omName;
#else
	string fieldName = fluidName + "/" + omField.getName();
#endif

	size_t avail = fieldNames.count(fieldName);

	int q(-1);
	if(avail > 0) {
		stringType::iterator stringIter;
		stringIter = fieldNames.find(fieldName);
		q = stringIter->second;
	}
	
	return q;
	
}


void ProblemType::name_User(Data &gdata)
{
	/* Standard naming scheme - different versions can be supplied by the
	   user
	*/

	// New version for single and multifluid:
#if(FLUID_TYPE == CRONOS_MULTIFLUID)
	numFluids = gdata.fluids->get_numFluids();
	// Loop over all fluids
	for(int iFluid=0; iFluid<gdata.fluids->get_numFluids(); ++iFluid) {
		CronosFluid fluid = gdata.fluids->fluids[iFluid];
		int q_min = fluid.get_IndexGlobal(0);
		int q_max = fluid.get_N_OMINT() + q_min;

#else
	CronosFluid fluid = gdata.fluid;
	int q_min = 0;
	int q_max = N_OMINT;
#endif
	for (int q=q_min; q<q_max; ++q) {
//		cout << " Making a loop " << q << " " << q_min << " " << q_max << endl;
//		cout << fluid.get_fieldName(q-q_min) << endl;
		gdata.om[q].rename(fluid.get_fieldName(q-q_min));

		string fluidName = fluid.get_Name();
		string fullName = fluidName + "/" +  gdata.om[q].getName();
//		cout << " fullNames: " << q << " " << fullName << endl;
		fieldNames[fullName] = q;
//		cout << " Setting the name " << fieldNames[fullName] << " " << fullName << endl;
//				fieldNames[gdata.om[q].getName()] = q;
	}

#if(FLUID_TYPE == CRONOS_MULTIFLUID)
	}
#endif

//	for(int q=0; q < N_OMINT; ++q) {
//		if(q == q_rho) {
//			gdata.om[q].rename("rho");
//		} else if (q == q_sx) {
//			gdata.om[q].rename("v_x");
//		} else if (q == q_sy) {
//			gdata.om[q].rename("v_y");
//		} else if (q == q_sz) {
//			gdata.om[q].rename("v_z");
//		} else if (q == q_Bx) {
//			gdata.om[q].rename("B_x");
//		} else if (q == q_By) {
//			gdata.om[q].rename("B_y");
//		} else if (q == q_Bz) {
//			gdata.om[q].rename("B_z");
//		} else if (q == q_Eges) {
//			gdata.om[q].rename("Etherm");
//		} else if (q == q_Eadd) {
//#if(AUX_ENERGY == ENTROPY)
//			gdata.om[q].rename("Entropy");
//#else
//			gdata.om[q].rename("Etherm");
//#endif
//		}
//
//		fieldNames[gdata.om[q].getName()] = q;
//	}


	for(int q=0; q < N_SUBS; ++q) {
		if(q == 0) {
			gdata.om[n_omInt+N_ADD+q].rename("A_x");
		} else if (q == 1) {
			gdata.om[n_omInt+N_ADD+q].rename("A_y");
		} else if (q == 2) {
			gdata.om[n_omInt+N_ADD+q].rename("A_z");
		}
	}

#if (OMS_USER == TRUE)
#if(FLUID_TYPE == CRONOS_MULTIFLUID)
	int n_omIntUser = gdata.fluids->get_N_OMINT_USER();
#else
	int n_omIntUser = N_OMINT_USER;
#endif


	for(int q=0; q<n_omIntUser; ++q) {

		char name_user[255];
		sprintf(name_user, "om_user%i",q);
		gdata.om_user[q].rename(name_user);

	}

#endif

}


void ProblemType::writeMovie(const Data &gdata, Movie & mov) const 
{
	for(int q=0; q<N_OMINT; ++q) {
		mov.writeStdSlices(gdata, gdata.om[q], gdata.om[q].getName());
	}

#if (OMS_USER == TRUE)
	for(int q=0; q<N_OMINT_USER; ++q) {
		mov.writeStdSlices(gdata, gdata.om_user[q], gdata.om_user[q].getName());
	}
#endif

		mov.writeStdSlices(gdata, gdata.computeAbs(1,3), "vabs");
		mov.writeStdSlices(gdata, gdata.computeAbs(4,6), "Babs");
}



bool ProblemType::force_max(const int &q) {
	return false;
}

bool ProblemType::force_min(const int &q) {
	return false;
}

REAL ProblemType::max_Val(const int &q) {
	return 1.e99;
}

REAL ProblemType::min_Val(const int &q) {
	return -1.e99;
}


void ProblemType::computeFluct(Data &gdata, double &ekfluc, double &ebfluc)
{
	REAL dV   = gdata.get_CellVolume(0,0,0);

	ekfluc = 0.;
	ebfluc = 0.;
	for (int k = 0; k < gdata.mx[2]; k++){
		for (int j = 0; j < gdata.mx[1]; j++){
			for (int i = 0; i < gdata.mx[0]; i++){

#if (NON_LINEAR_GRID == CRONOS_ON)
				dV = gdata.get_CellVolume(i,j,k);
#endif

				ekfluc += (sqr(gdata.om[q_sx](i,j,k)) + sqr(gdata.om[q_sy](i,j,k)) +
				           sqr(gdata.om[q_sz](i,j,k)))*0.5*dV*gdata.om[q_rho](i,j,k);
#if (FLUID_TYPE == CRONOS_MHD)
				ebfluc += (sqr(0.5*(gdata.om[q_Bx](i,j,k) + gdata.om[q_Bx](i-1,j,k))) +
				           sqr(0.5*(gdata.om[q_By](i,j,k) + gdata.om[q_By](i,j-1,k))) +
				           sqr(0.5*(gdata.om[q_Bz](i,j,k) + gdata.om[q_Bz](i,j,k-1))))*0.5*dV;
#endif

			}
		}
	}
}


REAL ProblemType::eta(Data &gdata,
                      const REAL &ii, const REAL &jj, const REAL &kk)
{
	return eta0;
}

REAL ProblemType::nu(Data &gdata,
                     const REAL &ii, const REAL &jj, const REAL &kk)
{
	return nu0;
}

REAL ProblemType::c2_iso(const Data &gdata,
                         const REAL &/*ii*/, const REAL &/*jj*/, const REAL &/*kk*/) const
{
	return cs2;
}



void ProblemType::set_q(int q, NumMatrix<REAL,3> &omField) {

	// Check if the Name is already in Use:
	size_t used = fieldNames.count(omField.getName());

	if(used) {
		
		typedef map<std::string, int> stringType;
		stringType::iterator stringIter;
		stringIter = fieldNames.find(omField.getName());

		cerr << omField.getName() << " already in uses as ";
		cerr << stringIter->second;
		cerr << endl;

	} else {

		fieldNames[omField.getName()] = q;

	}

	return;

}

void ProblemType::get_Cons(Data &gdata, cronos::vector<REAL> &Pos,
                           phys_fields_1D &physField, 
                           phys_fields_1D &physFieldsUser,
                           int dir, REAL shift) {
	//! Get conservative point values
	/*! Routine to compute conserved point value variables from characteristic
	  ones. Here we assume that they are the same
	  @PARAM physField generic fields
	  @PARAM physFieldsUser user fields
	 */
	for(int q=0; q<N_OMINT_USER; ++q) {
		physFieldsUser.uConORIG[q] = physFieldsUser.uPriORIG[q];
	}
}


REAL ProblemType::grid_user_x(REAL ratio) {
	return ratio;
}

REAL ProblemType::grid_user_y(REAL ratio) {
	return ratio;
}

REAL ProblemType::grid_user_z(REAL ratio) {
	return ratio;
}

void ProblemType::bc_User(Data &gdata, NumMatrix<REAL,3> &omb,
                          int dir, int top, int q, int rim) {
	cerr << " Error: User boundaries at: " << endl;
	cerr << "  direction: " << dir << " top: " << top << endl;
	cerr << " set in cat file - but not defined in mod-file " << endl;
	cerr << " EXITING " << endl;;
#ifdef parallel
	MPI_Finalize();
#endif
      exit(-33);

}

int ProblemType::checkConvergence(Data &gdata) {
	//! Method can be used to check whether convergence has been reached
	return 0;
}
