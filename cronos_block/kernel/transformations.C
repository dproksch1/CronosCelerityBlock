#include "transformations.H"
#include "physical_constants.H"
#include "gridfunc.H"
#include <stdlib.h>


Transformations::Transformations(const CronosFluid &fluid, ProblemType &Problem, bool TPhys, int iFluid) {
	TempErr = 0;

	// Copying some constants from Problem-class to use locally
	this->q_rho = fluid.get_q_rho_global();
	this->q_sx = fluid.get_q_sx_global();
	this->q_sy = fluid.get_q_sy_global();
	this->q_sz = fluid.get_q_sz_global();
	this->q_Bx = fluid.get_q_Bx_global();
	this->q_By = fluid.get_q_By_global();
	this->q_Bz = fluid.get_q_Bz_global();
	this->q_Eges = fluid.get_q_Eges_global();
	this->q_Eadd = fluid.get_q_Eadd_global();

	this->q_rho_loc = fluid.get_q_rho();
	this->q_sx_loc = fluid.get_q_sx();
	this->q_sy_loc = fluid.get_q_sy();
	this->q_sz_loc = fluid.get_q_sz();
	this->q_Bx_loc = fluid.get_q_Bx();
	this->q_By_loc = fluid.get_q_By();
	this->q_Bz_loc = fluid.get_q_Bz();
	this->q_Eges_loc = fluid.get_q_Eges();
	this->q_Eadd_loc = fluid.get_q_Eadd();

	this->fluidType = fluid.get_fluid_type();
	this->iFluid = iFluid;
	this->magFluid = fluid.has_MagField();

	magFluid = false;

	this->TPhys = TPhys;

	if(ENERGETICS == FULL && this->TPhys) {
		TNorm = value((char*)"TempNorm");   // K
		if(Problem.TrafoNorm != NULL) {
			TNorm = Problem.TrafoNorm->get_num(Problem.TrafoNorm->TEMP, 1.0*CRONOS_CONSTANTS::Kelvin);
		}
	} else {
		TNorm = 1.;
	}

	if(ENERGETICS == FULL) {
		thermal = static_cast<int>(value((char*)"thermal"));
	}

#if (USE_COROTATION == CRONOS_ON)
	omegaZ = value((char*)"omegaZ");
#endif

	if(Problem.TrafoNorm != NULL) {
		kBoverMeanMolWeight_num = Problem.TrafoNorm->get_num(CRONOS_CONSTANTS::BoltzmannConstant / Problem.meanParticleMass);
	} else {
		kBoverMeanMolWeight_num = 1.;
	}
}

void Transformations::set_thermal(bool thermal) {
	this->thermal = thermal;
}


void Transformations::TransMomen2Vel(Data &gdata, gridFunc &gfunc,
                                     ProblemType &Problem)
{
	if(gdata.om[q_sx].getName() == "v_x" || 
	   gdata.om[q_sy].getName() == "v_y" ||
	   gdata.om[q_sz].getName() == "v_z") {
		throw CException(" Velocity instead of momentum ");
	}

	for (int q = q_sx; q <= q_sz; ++q) {
		gdata.om[q] /= gdata.om[q_rho];
	}

	gdata.om[q_sx].rename("v_x");
	gdata.om[q_sy].rename("v_y");
	gdata.om[q_sz].rename("v_z");
}



void Transformations::TransVel2Momen(Data &gdata, gridFunc &gfunc,
                                     ProblemType &Problem)
{
	if(gdata.om[q_sx].getName() == "s_x" || 
	   gdata.om[q_sy].getName() == "s_y" ||
	   gdata.om[q_sz].getName() == "s_z") {
		throw CException(" Momentum instead of velocity ");
	}

	for (int q = q_sx; q <= q_sz; ++q) {
		gdata.om[q] *= gdata.om[q_rho];
	}

	gdata.om[q_sx].rename("s_x");
	gdata.om[q_sy].rename("s_y");
	gdata.om[q_sz].rename("s_z");

}


void Transformations::TransEth2E(const Data &gdata, gridFunc &gfunc,
                                 ProblemType &Problem) const
{

	if(gdata.om[q_Eges].getName() != "Etherm") {
		cerr << " Energy is: " << gdata.om[q_Eges].getName() << " " << q_Eges << endl;
		throw CException(" Transformations::TransEth2E - om[q_Eges] is not set as thermal energy ");
	}

	if(Problem.gamma < 1.0000000001) {
		throw CException(" Must not be isothermal ");
	}

	if(gdata.om[q_sx].getName() == "v_x" && gdata.om[q_sy].getName() == "v_y" &&
	   gdata.om[q_sz].getName() == "v_z") {

		for(int k = -B+1; k<=gdata.mx[2]+B; ++k){
			for(int j = -B+1; j<=gdata.mx[1]+B; ++j){
				for(int i = -B+1; i<=gdata.mx[0]+B; ++i){				  
					//   add e_kin
					gdata.om[q_Eges](i,j,k) += 0.5*(sqr(gdata.om[q_sx](i,j,k)) + 
					                                sqr(gdata.om[q_sy](i,j,k)) +
					                                sqr(gdata.om[q_sz](i,j,k)))/gdata.om[q_rho](i,j,k);
				}
			}
		}
	}
	gdata.om[q_Eges].rename("Eges");
}


void Transformations::TransE2Eth(Data &gdata, gridFunc &gfunc,
                                 ProblemType &Problem, int n, bool DualEnergy)
{
	if(gdata.om[q_Eges].getName() != "Eges") {
		throw CException(" om[7] is not set as overall energy ");
	}
    
	if(Problem.gamma < 1.0000000001) {
		throw CException(" Must not be isothermal ");
	}
#if(CRSWITCH_DUAL_ENERGY == CRONOS_OFF)
	DualEnergy = false;
#endif

	int n_omIntAll = N_OMINT;

	// Saving overall energy:
	if(ENERGETICS == FULL && DualEnergy) {
		gdata.om[n_omIntAll] = gdata.om[q_Eges];
		gdata.om[n_omIntAll].rename(gdata.om[q_Eges].getName());
	}

	TempErr = 0;

	if(gdata.om[q_sx].getName() == "v_x" && gdata.om[q_sy].getName() == "v_y" &&
	   gdata.om[q_sz].getName() == "v_z") {
		
		// CRONOS_CRIT: Formerly done only in computational domain
		// The idea was that values in ghost cells are provided by BCs
		// If no BCs are given, however, the option is to do computation
		// everywhere and to ignore the BCs

		// This, however, does not work, when there are significant changes for the total
		// energy near the boundary...

		 for(int k = -B+1; k<=gdata.mx[2]+B; ++k){
		 	for(int j = -B+1; j<=gdata.mx[1]+B; ++j){
		 		for(int i = -B+1; i<=gdata.mx[0]+B; ++i){

					REAL E_kin = 0.5*(sqr(gdata.om[q_sx](i,j,k)) +
					                  sqr(gdata.om[q_sy](i,j,k)) +
					                  sqr(gdata.om[q_sz](i,j,k)))*gdata.om[q_rho](i,j,k);

					// double del_vx_dx = 0.5*(gdata.om[q_sx](i+1,j,k) -
					//                         gdata.om[q_sx](i-1,j,k));

					gdata.om[q_Eges](i,j,k) -= E_kin; //  sub e_kin
				}
			}
		}
		// Check of negative energy only in computational domain - not in ghost cells
		for(int k = 0; k<=gdata.mx[2]; ++k){
			for(int j = 0; j<=gdata.mx[1]; ++j){
				for(int i = 0; i<=gdata.mx[0]; ++i){

	        
#if(CRSWITCH_DUAL_ENERGY == CRONOS_OFF)
					if(!Problem.force_min(q_Eges) && 
					   gdata.om[q_Eges](i,j,k) < 0.){
						cerr << " Transformations::TransE2Eth " << endl;
//						cerr << q_rho << " " << q_sx << " " << q_Eges << endl;
						cerr << " Error! negative energy at cell ("
							  << i << " " << j << " " << k << "), pos ("
							  << gdata.getCen_x(i) << " "
							  << gdata.getCen_y(j) << " "
							  << gdata.getCen_z(k) << "); "
							  << gdata.om[q_Eges].getName() << " = "
							  << gdata.om[q_Eges](i,j,k) << endl;
						cerr << "(rho, V) = ("
							  << gdata.om[q_rho](i,j,k) << " "
							  << gdata.om[q_sx](i,j,k) << " "
							  << gdata.om[q_sy](i,j,k) << " "
							  << gdata.om[q_sz](i,j,k) << ")";
						if(magFluid) {
							cerr << "; " << endl << "B = ("
								  << gdata.om[q_Bx](i,j,k) << " "
								  << gdata.om[q_Bx](i-1,j,k) << " "
								  << gdata.om[q_By](i,j,k) << " "
								  << gdata.om[q_By](i,j-1,k) << " "
								  << gdata.om[q_Bz](i,j,k) << " "
								  << gdata.om[q_Bz](i,j,k-1) << ").";
						}
						cerr << endl;
						TempErr     += 1;
						exit(3);
					}
#else
					if(!Problem.force_min(q_Eges) && 
					   gdata.om[q_Eges](i,j,k) < 0.){
						TempErr += 1;
//						cout << i << " " << j << " " << k << " " << gdata.om[q_Eges](i,j,k) << endl;
					}
#endif

				}
			}
		}
	} else {
		// for(int k = -B+1; k<=gdata.mx[2]+B; ++k){
		// 	for(int j = -B+1; j<=gdata.mx[1]+B; ++j){
		// 		for(int i = -B+1; i<=gdata.mx[0]+B; ++i){
		for(int k = 0; k<=gdata.mx[2]; ++k){
			for(int j = 0; j<=gdata.mx[1]; ++j){
				for(int i = 0; i<=gdata.mx[0]; ++i){

					gdata.om[q_Eges](i,j,k) -= 0.5*(sqr(gdata.om[q_sx](i,j,k)) +    //  sub e_kin
					                                sqr(gdata.om[q_sy](i,j,k)) +
					                                sqr(gdata.om[q_sz](i,j,k)))/gdata.om[q_rho](i,j,k);

					if(!Problem.force_min(q_Eges) &&
					   gdata.om[q_Eges](i,j,k) < 0){
						TempErr     += 1;
//						cout << i << " " << j << " " << k << " " << gdata.om[q_Eges](i,j,k) << endl;

					}
				}
			}
		}
	}

#if(CRSWITCH_DUAL_ENERGY == CRONOS_OFF)
	if(TempErr > 0) {
		string ErrString = " ERROR: Temperature < 0 at ";
		char cerrs[255];
		sprintf(cerrs,"%i",TempErr);
		ErrString += cerrs;
		ErrString += " Gridpoints ";
		throw CException(ErrString);	
	}
#endif


	// cout << " Etherm(trafo): " << gdata.om[q_Eges](1,5,0) << " ";
	// cout << gdata.om[q_Eges](0,5,0) << " ";
	// cout << gdata.om[q_Eges](-1,5,0) << " ";
	// cout << endl;

	// Store thermal pressure if necessary
	if(gdata.storePressure) {
		for(int k = 0; k<=gdata.mx[2]; ++k){
			for(int j = 0; j<=gdata.mx[1]; ++j){
				for(int i = 0; i<=gdata.mx[0]; ++i){
					gdata.pTherm(i,j,k) = (Problem.gamma-1)*gdata.om[q_Eges](i,j,k);
				}
			}
		}
		gfunc.boundary(gdata, Problem, gdata.pTherm,3);
	}


	gdata.om[q_Eges].rename("Etherm");
	gfunc.boundary(gdata, Problem, gdata.om[q_Eges],B,q_Eges, iFluid);


	// cout << " Etherm(bc): " << gdata.om[q_Eges](1,5,0) << " ";
	// cout << gdata.om[q_Eges](0,5,0) << " ";
	// cout << gdata.om[q_Eges](-1,5,0) << " ";
	// cout << endl;

	
	if(DualEnergy) {
#if(CRSWITCH_DUAL_ENERGY == CRONOS_ON)
		EntropyCorrection(gdata, gfunc, Problem, n);
#endif
	} else {
		if(TempErr > 0) {
			throw CException(" T < 0");
		}
	}
}


void Transformations::TransE2Eth(Data &gdata, gridFunc &gfunc,
                                 ProblemType &Problem) {
	TransE2Eth(gdata, gfunc, Problem, 0, false);
}


void Transformations::TransT2Eth(const Data &gdata, gridFunc &gfunc,
                                 ProblemType &Problem) const
{

//	static const double kBoverMu_num = Problem.TrafoNorm->get_num(CRONOS_CONSTANTS::BoltzmannConstant / Problem.meanParticleMass);

	if(gdata.om[q_Eges].getName() != "Temp") {
		throw CException(" om[q_Eges] is not set as temperature ");
	}

	if(Problem.gamma < 1.0000000001) {
		throw CException(" Must not be isothermal ");
	}
//	 REAL fac = 1./((Problem.gamma - 1.)*TNorm);
//	 cout << " factor " << fac << " " << gdata.om[q_Eges](12,12,12) << endl;
//	REAL fac = kBoverMu_num / (Problem.gamma - 1);

	REAL fac = kBoverMeanMolWeight_num / (Problem.gamma - 1.);
	fac /= TNorm;

	for(int k = -B; k<=gdata.mx[2]+B; ++k){
		for(int j = -B; j<=gdata.mx[1]+B; ++j){
			for(int i = -B; i<=gdata.mx[0]+B; ++i){
        
				// T -> e_th
				gdata.om[q_Eges](i,j,k) *= gdata.om[q_rho](i,j,k)*fac;

			}
		}
	}
	gdata.om[q_Eges].rename("Etherm");
}


void Transformations::TransT2E(const Data &gdata, gridFunc &gfunc,
                               ProblemType &Problem) const
{
	if(gdata.om[q_Eges].getName() != "Temp") {
		throw CException(" om[q_Eges] is not set as temperature ");
	}
	TransT2Eth(gdata, gfunc, Problem);
	TransEth2E(gdata, gfunc, Problem);

}

REAL Transformations::TransEth2T(Data &gdata, gridFunc &gfunc,
                                 ProblemType &Problem) 
{
	if(gdata.om[q_Eges].getName() != "Etherm") {
		throw CException(" om[q_Eges] is not set as thermal energy ");
	}

	if(Problem.gamma < 1.0000000001) {
		throw CException(" Must not be isothermal ");
	}

	REAL fac = TNorm*(Problem.gamma-1.)/kBoverMeanMolWeight_num;

	for(int k = -B; k<=gdata.mx[2]+B; ++k){
		for(int j = -B; j<=gdata.mx[1]+B; ++j){
			for(int i = -B; i<=gdata.mx[0]+B; ++i){
				// e_th -> T
//				gdata.om[q_Eges](i,j,k) *= TNorm*(Problem.gamma-1.)/gdata.om[q_rho](i,j,k);

				gdata.om[q_Eges](i,j,k) *= fac/gdata.om[q_rho](i,j,k);

				if(!Problem.force_min(q_Eges) && gdata.om[q_Eges](i,j,k) < 0){
					TempErr += 1;
				}
			}
		}
	}

	gdata.om[q_Eges].rename("Temp");
	gfunc.boundary(gdata, Problem, gdata.om[q_Eges],B,q_Eges, iFluid);

	REAL TAve(1.);
	if(TPhys) {
		TAve = gdata.computeRMS(q_Eges);
	}
	return TAve;
}