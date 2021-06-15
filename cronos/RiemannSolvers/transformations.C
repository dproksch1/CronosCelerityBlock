#include "transformations.H"
#include "physical_constants.H"
#include <stdlib.h>


//Transformations::Transformations(ProblemType &Problem, bool TPhys) {
//	TempErr = 0;
//
//	// Copying some constants from Problem-class to use locally
//	this->q_rho = Problem.q_rho;
//	this->q_sx = Problem.q_sx;
//	this->q_sy = Problem.q_sy;
//	this->q_sz = Problem.q_sz;
//	this->q_Bx = Problem.q_Bx;
//	this->q_By = Problem.q_By;
//	this->q_Bz = Problem.q_Bz;
//	this->q_Eges = Problem.q_Eges;
//	this->q_Eadd = Problem.q_Eadd;
//	this->iFluid = 0;
//
//	this->TPhys = TPhys;
//#if(FLUID_TYPE == CRONOS_MHD)
//	magFluid = true;
//#elif(FLUID_TYPE == CRONOS_HYDRO)
//	magFluid = false;
//#endif
//
//	if(ENERGETICS == FULL && this->TPhys) {
//		TNorm = value((char*)"TempNorm");   // K
//		TNorm = Problem.TrafoNorm->get_num(Problem.TrafoNorm->TEMP, 1*CRONOS_CONSTANTS::Kelvin);
//	} else {
//		TNorm = 1.;
//	}
//
//	if(ENERGETICS == FULL) {
//		thermal = static_cast<int>(value((char*)"thermal"));
//	}
//
//#if (USE_COROTATION == CRONOS_ON)
//	omegaZ = value((char*)"omegaZ");
//#endif
//
//
//}

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

#if(FLUID_TYPE == CRONOS_MHD)
	magFluid = true;
#elif(FLUID_TYPE == CRONOS_HYDRO)
	magFluid = false;
#endif

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

//void Transformations::reset_Indices(ProblemType &Problem) {
//	//! Reset field indices (only necessary for multifluid simulations)
//	this->q_rho = Problem.q_rho;
//	this->q_sx = Problem.q_sx;
//	this->q_sy = Problem.q_sy;
//	this->q_sz = Problem.q_sz;
//	this->q_Bx = Problem.q_Bx;
//	this->q_By = Problem.q_By;
//	this->q_Bz = Problem.q_Bz;
//	this->q_Eges = Problem.q_Eges;
//	this->q_Eadd = Problem.q_Eadd;
//}

void Transformations::reset_Indices(CronosFluid &fluid) {
	//! Reset field indices (only necessary for multifluid simulations)
	this->q_rho = fluid.get_q_rho_global();
	this->q_sx = fluid.get_q_sx_global();
	this->q_sy = fluid.get_q_sy_global();
	this->q_sz = fluid.get_q_sz_global();
	this->q_Bx = fluid.get_q_Bx_global();
	this->q_By = fluid.get_q_By_global();
	this->q_Bz = fluid.get_q_Bz_global();
	this->q_Eges = fluid.get_q_Eges_global();
	this->q_Eadd= fluid.get_q_Eadd_global();
}

void Transformations::set_thermal(bool thermal) {
	this->thermal = thermal;
}


void Transformations::TransPrim2Cons(Data &gdata, gridFunc &gfunc,
                                     ProblemType &Problem) {
	//! Transform from primitive to conservative variables

	if (ENERGETICS == FULL) {
		if(gdata.om[q_Eges].getName() == "Etherm") {
			TransEth2E(gdata, gfunc, Problem);
		} else if(gdata.om[q_Eges].getName() == "Temp") {
			TransT2E(gdata, gfunc, Problem);
		}
	}

	if(gdata.om[q_sx].getName() == "v_x" && 
	   gdata.om[q_sy].getName() == "v_y" &&
	   gdata.om[q_sz].getName() == "v_z") {

#if (USE_ANGULAR_MOMENTUM == TRUE)
		TransVel2AngMom(gdata, gfunc, Problem);
#else
		TransVel2Momen(gdata, gfunc, Problem);
#endif

	}
}


void Transformations::TransCons2Prim(Data &gdata, gridFunc &gfunc,
                                     ProblemType &Problem) {
	//! Transform from conservative to primitive variables
#if (USE_ANGULAR_MOMENTUM == TRUE)
	TransAngMom2Vel(gdata, gfunc, Problem);
#else
	TransMomen2Vel(gdata, gfunc, Problem);
#endif

	if(ENERGETICS == FULL) {
		if(thermal) {
			TransE2Eth(gdata, gfunc, Problem);
		} else {
			TransE2T(gdata, gfunc, Problem);
		}
	}

}


void Transformations::TransE(Data &gdata, gridFunc &gfunc, ProblemType &Problem,
                             string NewType) {
	string OldType = gdata.om[q_Eges].getName();
	if(NewType == OldType) {
		return;
	}

	if(OldType == "Eges") {
		if(NewType == "Etherm") {
			TransE2Eth(gdata, gfunc, Problem);
		} else if (NewType == "Temp") {
			TransE2T(gdata, gfunc, Problem);
		} else {
			throw(" No such energy form implemented ");
		}
	} else if(OldType == "Etherm") {
		if(NewType == "Eges") {
			TransEth2E(gdata, gfunc, Problem);
		} else if (NewType == "Temp") {
			TransEth2T(gdata, gfunc, Problem);
		} else {
			throw(" No such energy form implemented ");
		}
	} else if(OldType == "Temp") {
		if(NewType == "Eges") {
			TransT2E(gdata, gfunc, Problem);
		} else if (NewType == "Etherm") {
			TransT2Eth(gdata, gfunc, Problem);
		} else {
			throw(" No such energy form implemented ");
		}
	}    
}


void Transformations::TransE2Eth(Data &gdata, gridFunc &gfunc,
                                 ProblemType &Problem) {
	TransE2Eth(gdata, gfunc, Problem, 0, false);
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

#if(FLUID_TYPE == CRONOS_MULTIFLUID)
	int n_omIntAll = gdata.fluids->get_N_OMINT_ALL();
#else
	int n_omIntAll = N_OMINT;
#endif

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
//		for(int k = 0; k<=gdata.mx[2]; ++k){
//			for(int j = 0; j<=gdata.mx[1]; ++j){
//				for(int i = 0; i<=gdata.mx[0]; ++i){

					REAL E_kin = 0.5*(sqr(gdata.om[q_sx](i,j,k)) +
					                  sqr(gdata.om[q_sy](i,j,k)) +
					                  sqr(gdata.om[q_sz](i,j,k)))*gdata.om[q_rho](i,j,k);

					// double del_vx_dx = 0.5*(gdata.om[q_sx](i+1,j,k) -
					//                         gdata.om[q_sx](i-1,j,k));

					gdata.om[q_Eges](i,j,k) -= E_kin; //  sub e_kin

#if (FLUID_TYPE != CRONOS_HYDRO)
#if(FLUID_TYPE == CRONOS_MULTIFLUID)
					if(magFluid) {
#endif
					
					// REAL E_mag = 0.5*(sqr(0.5*(gdata.om[q_Bx](i,j,k) +
					//                            gdata.om[q_Bx](i-1,j,k))) +
					//                   sqr(0.5*(gdata.om[q_By](i,j,k) +
					//                            gdata.om[q_By](i,j-1,k))) +
					//                   sqr(0.5*(gdata.om[q_Bz](i,j,k) +
					//                            gdata.om[q_Bz](i,j,k-1))));
					REAL E_mag = get_EMag(gdata, i, j, k);

					// double del_Bx_dx = (gdata.om[q_Bx](i,j,k) - 
					//                     gdata.om[q_Bx](i-1,j,k));
					// double del_By_dy = (gdata.om[q_By](i,j,k) -
					//                     gdata.om[q_By](i,j-1,k));

					gdata.om[q_Eges](i,j,k) -= E_mag; //  sub e_mag

#if(FLUID_TYPE == CRONOS_MULTIFLUID)
					}
#endif
#endif
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
#ifdef parallel
						cerr << " for rank: " << gdata.rank << endl;
#endif
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
#if (FLUID_TYPE != CRONOS_HYDRO)
#if(FLUID_TYPE == CRONOS_MULTIFLUID)
					if(magFluid) {
#endif
					//  sub e_mag
					gdata.om[q_Eges](i,j,k) -= get_EMag(gdata, i, j, k);
#if(FLUID_TYPE == CRONOS_MULTIFLUID)
					}
#endif
#endif

					if(!Problem.force_min(q_Eges) &&
					   gdata.om[q_Eges](i,j,k) < 0){
						TempErr     += 1;
//						cout << i << " " << j << " " << k << " " << gdata.om[q_Eges](i,j,k) << endl;

					}
				}
			}
		}
	}

#ifdef parallel  
	MPI_Barrier(gdata.comm3d);
	gdata.MpiSum(TempErr);
	MPI_Barrier(gdata.comm3d);
#endif

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

#if (FLUID_TYPE != CRONOS_HYDRO)
#if(FLUID_TYPE == CRONOS_MULTIFLUID)
					if(magFluid) {
#endif
					//   add e_mag
					REAL E_mag = get_EMag(gdata, i, j, k);
					// double del_Bx_dx = (gdata.om[q_Bx](i,j,k) -
					//                     gdata.om[q_Bx](i-1,j,k));

					gdata.om[q_Eges](i,j,k) += E_mag;
#if(FLUID_TYPE == CRONOS_MULTIFLUID)
					}
#endif
#endif
	        
					//   add e_kin
					gdata.om[q_Eges](i,j,k) += 0.5*(sqr(gdata.om[q_sx](i,j,k)) + 
					                                sqr(gdata.om[q_sy](i,j,k)) +
					                                sqr(gdata.om[q_sz](i,j,k)))*gdata.om[q_rho](i,j,k);
				}
			}
		}
	} else {
		for(int k = -B+1; k<=gdata.mx[2]+B; ++k){
			for(int j = -B+1; j<=gdata.mx[1]+B; ++j){
				for(int i = -B+1; i<=gdata.mx[0]+B; ++i){

#if (FLUID_TYPE != CRONOS_HYDRO)
#if(FLUID_TYPE == CRONOS_MULTIFLUID)
					if(magFluid) {
#endif
					//   add e_mag
					gdata.om[q_Eges](i,j,k) += get_EMag(gdata, i, j, k);
#if(FLUID_TYPE == CRONOS_MULTIFLUID)
					}
#endif
#endif
				  
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


REAL Transformations::TransE2T(Data &gdata, gridFunc &gfunc,
                               ProblemType &Problem)
{

	if(gdata.om[q_Eges].getName() != "Eges") {
		throw CException(" om[q_Eges] is not set as overall energy ");
	}
	TransE2Eth(gdata, gfunc, Problem);
	return TransEth2T(gdata, gfunc, Problem);

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


#if (USE_COROTATION == CRONOS_ON)
//void Transformations::setOmega(double omegaZ) {
//	//! Set value for angular velocity of co-rotating frame
//	this->omegaZ = omegaZ;
//}

void Transformations::TransCorotToInert(Data &gdata, gridFunc &gfunc, ProblemType &problem) {
	//! Transform velocity from corotating to inertial frame
	/*!
	 * Here we need to distinguish between different coordinate systems (see below)
	 * Beware: angular velocity is assumed to be in z-direction
	 * */

	if(gdata.om[q_sx].getName() == "v_x" ||
	   gdata.om[q_sy].getName() == "v_y" ||
	   gdata.om[q_sz].getName() == "v_z") {
		throw CException(" Velocity instead of corotating velocity");
	}


#if(GEOM == CARTESIAN)
	Pot & ux = gdata.om[q_sx];
	Pot & uy = gdata.om[q_sy];
#elif(GEOM == CYLINDRICAL)
	Pot & uPhi = gdata.om[q_sy];
#else
	Pot & uPhi = gdata.om[q_sz];
#endif

	// Loop over entire grid:
	for(int iz = -B; iz<=gdata.mx[2]+B; ++iz){
		for(int iy = -B; iy<=gdata.mx[1]+B; ++iy){
			REAL yPos = gdata.getCen_y(iy);
			for(int ix = -B; ix<=gdata.mx[0]+B; ++ix){
				// Need different equations for the different grid types.
#if(GEOM == CARTESIAN)
				REAL xPos = gdata.getCen_x(ix);

				ux(ix, iy, iz) -= yPos*omegaZ;
				uy(ix, iy, iz) += xPos*omegaZ;

#elif(GEOM == CYLINDRICAL)

				REAL rhoPos = gdata.getCen_x(ix);

				uPhi(ix, iy, iz) += rhoPos*omegaZ;

#else

				REAL rad = gdata.getCen_x(ix);

				uPhi(ix, iy, iz) += rad*sin(yPos)*omegaZ;

#endif
			}
		}
	}

	gdata.om[q_sx].rename("v_x");
	gdata.om[q_sy].rename("v_y");
	gdata.om[q_sz].rename("v_z");

}

void Transformations::TransInertToCorot(Data &gdata, gridFunc &gfunc, ProblemType &problem) {
	//! Transform velocity from inertial to corotating frame
	/*!
	 * Here we need to distinguish between different coordinate systems (see below)
	 * Beware: angular velocity is assumed to be in z-direction
	 * */

	if(gdata.om[q_sx].getName() == "v_x_Corot" ||
	   gdata.om[q_sy].getName() == "v_y_Corot" ||
	   gdata.om[q_sz].getName() == "v_z_Corot") {
		throw CException(" Corotating velocity instead of inertial frame velocity");
	}


#if(GEOM == CARTESIAN)
	Pot & ux = gdata.om[q_sx];
	Pot & uy = gdata.om[q_sy];
#elif(GEOM == CYLINDRICAL)
	Pot & uPhi = gdata.om[q_sy];
#else
	Pot & uPhi = gdata.om[q_sz];
#endif

	// Loop over entire grid:
	for(int iz = -B; iz<=gdata.mx[2]+B; ++iz){
		for(int iy = -B; iy<=gdata.mx[1]+B; ++iy){
			REAL yPos = gdata.getCen_y(iy);
			for(int ix = -B; ix<=gdata.mx[0]+B; ++ix){
				// Need different equations for the different grid types.
#if(GEOM == CARTESIAN)
				REAL xPos = gdata.getCen_x(ix);

				ux(ix, iy, iz) += yPos*omegaZ;
				uy(ix, iy, iz) -= xPos*omegaZ;

#elif(GEOM == CYLINDRICAL)

				REAL rhoPos = gdata.getCen_x(ix);

				uPhi(ix, iy, iz) -= rhoPos*omegaZ;

#else

				REAL rad = gdata.getCen_x(ix);

				uPhi(ix, iy, iz) -= rad*sin(yPos)*omegaZ;

#endif
			}
		}
	}

	gdata.om[q_sx].rename("v_x_Corot");
	gdata.om[q_sy].rename("v_y_Corot");
	gdata.om[q_sz].rename("v_z_Corot");

}

REAL Transformations::TransCorotToInert_x(Data &gdata, REAL vCorot, int ix, int iy, int iz) {
	//! Local transform from corotating to inertial frame
	/*!
	 * transformation of x-component
	 * */
#if(GEOM == CARTESIAN)
	REAL yPos = gdata.getCen_y(iy);
	return vCorot - yPos*omegaZ;
#else
	return vCorot;
#endif

}

REAL Transformations::TransCorotToInert_y(Data &gdata, REAL vCorot, int ix, int iy, int iz) {
	//! Local transform from corotating to inertial frame
	/*!
	 * transformation of y-component
	 * */
#if(GEOM == CARTESIAN)
	REAL xPos = gdata.getCen_x(ix);
	return vCorot + xPos*omegaZ;
#elif(GEOM == CYLINDRICAL)
	REAL r_cyl = gdata.getCen_x(ix);
	return vCorot + r_cyl*omegaZ;
#else
	return vCorot;
#endif

}

REAL Transformations::TransCorotToInert_z(Data &gdata, REAL vCorot, int ix, int iy, int iz) {
	//! Local transform from corotating to inertial frame
	/*!
	 * transformation of z-component
	 * */
#if(GEOM == CARTESIAN)
	return vCorot;
#elif(GEOM == CYLINDRICAL)
	return vCorot;
#else
	REAL rad = gdata.getCen_x(ix);
	REAL theta = gdata.getCen_y(iy);
	return vCorot + rad*sin(theta)*omegaZ;
#endif

}



REAL Transformations::TransInertToCorot_x(Data &gdata, REAL vInert, int ix, int iy, int iz) {
	//! Local transform from inertial to corotating frame
	/*!
	 * transformation of x-component
	 * */
#if(GEOM == CARTESIAN)
	REAL yPos = gdata.getCen_y(iy);
	return vInert + yPos*omegaZ;
#else
	return vInert;
#endif

}


REAL Transformations::TransInertToCorot_y(Data &gdata, REAL vInert, int ix, int iy, int iz) {
	//! Local transform from inertial to corotating frame
	/*!
	 * transformation of y-component
	 * */
#if(GEOM == CARTESIAN)
	REAL xPos = gdata.getCen_x(ix);
	return vInert- xPos*omegaZ;
#elif(GEOM == CYLINDRICAL)
	REAL r_cyl = gdata.getCen_x(ix);
	return vInert - r_cyl*omegaZ;
#else
	return vInert;
#endif

}

REAL Transformations::TransInertToCorot_z(Data &gdata, REAL vInert, int ix, int iy, int iz) {
	//! Local transform from inertial to corotating frame
	/*!
	 * transformation of z-component
	 * */
#if(GEOM == CARTESIAN)
	return vInert;
#elif(GEOM == CYLINDRICAL)
	return vInert;
#else
	REAL rad = gdata.getCen_x(ix);
	REAL theta = gdata.getCen_y(iy);
	return vInert - rad*sin(theta)*omegaZ;
#endif

}





void Transformations::src_Corotating(Data &gdata, ProblemType &problem, NumMatrix<REAL, 3> nom[]) {
	//! Add corotation source term -\rho (\vec{\Omega} \times \vec{u})
	/*!
	 * Beware: u_1 -> u_x and u_2 -> u_y for Cartesian coordinates, while
	 * u_1 -> u_rcyl and u_2 -> u_phi for cylindrical ones.
	 * */
  Pot & rho = gdata.om[q_rho];
#if(GEOM != SPHERICAL)
	Pot & u_1 = gdata.om[q_sx];
	Pot & u_2 = gdata.om[q_sy];
#else
	Pot & uRad = gdata.om[q_sx];
	Pot & uThe = gdata.om[q_sy];
	Pot & uPhi = gdata.om[q_sz];
#endif

	// Loop over numerical domain (without ghost cells)
	for(int iz = 0; iz<=gdata.mx[2]; ++iz){
		for(int iy = 0; iy<=gdata.mx[1]; ++iy){
#if(GEOM==SPHERICAL)
			REAL theta = gdata.getCen_y(iy);
#endif
			for(int ix = 0; ix<=gdata.mx[0]; ++ix){
#if(GEOM!=SPHERICAL)
			  nom[q_sx](ix,iy,iz) -= rho(ix,iy,iz)*u_2(ix,iy,iz)*omegaZ;
			  nom[q_sy](ix,iy,iz) += rho(ix,iy,iz)*u_1(ix,iy,iz)*omegaZ;
#else // Spherical case
			  nom[q_sx](ix,iy,iz) -= rho(ix,iy,iz)*sin(theta)*uPhi(ix,iy,iz)*omegaZ; // change of radial velocity
			  nom[q_sy](ix,iy,iz) -= rho(ix,iy,iz)*cos(theta)*uPhi(ix,iy,iz)*omegaZ; // change of theta component
			  nom[q_sz](ix,iy,iz) += rho(ix,iy,iz)*(sin(theta)*uRad(ix,iy,iz) +
									  cos(theta)*uThe(ix,iy,iz))*omegaZ; // change of phi component
#endif

			}
		}
	}

}



void Transformations::store_uInert(Data &gdata, cronos::vector<REAL> &Pos, phys_fields_1D &fields, int dir) {
	//! Store inertial frame velocity frame:
	fields.uInertial[0] = fields.uPri[q_sx_loc];
	fields.uInertial[1] = fields.uPri[q_sy_loc];
	fields.uInertial[2] = fields.uPri[q_sz_loc];

	int ix, iy, iz;
	ix = Pos.get(0);
	iy = Pos.get(1);
	iz = Pos.get(2);


	// Compute velocity in co-rotating frame
	for (int i = -2; i <= gdata.mx[dir]+1; ++i){
		if(dir==0) {
			ix = i;
		} else if(dir==1) {
			iy = i;
		} else {
			iz = i;
		}
		fields.uPri[q_sx_loc](i) = TransInertToCorot_x(gdata, fields.uInertial[0](i), ix, iy, iz);
		fields.uPri[q_sy_loc](i) = TransInertToCorot_y(gdata, fields.uInertial[1](i), ix, iy, iz);
		fields.uPri[q_sz_loc](i) = TransInertToCorot_z(gdata, fields.uInertial[2](i), ix, iy, iz);
		// fields.uPri[q_sx_loc](i) = TransCorotToInert_x(gdata, fields.uInertial[0](i), ix, iy, iz);
		// fields.uPri[q_sy_loc](i) = TransCorotToInert_y(gdata, fields.uInertial[1](i), ix, iy, iz);
		// fields.uPri[q_sz_loc](i) = TransCorotToInert_z(gdata, fields.uInertial[2](i), ix, iy, iz);
	}
}


void Transformations::store_uInert(Data &gdata, phys_fields_0D &fields, int ix, int iy, int iz) {
	//! Store inertial frame velocity frame:
	fields.uInertial[0] = fields.uPri[q_sx_loc];
	fields.uInertial[1] = fields.uPri[q_sy_loc];
	fields.uInertial[2] = fields.uPri[q_sz_loc];


	// Compute velocity in co-rotating frame
	fields.uPri[q_sx_loc] = TransInertToCorot_x(gdata, fields.uInertial[0], ix, iy, iz);
	fields.uPri[q_sy_loc] = TransInertToCorot_y(gdata, fields.uInertial[1], ix, iy, iz);
	fields.uPri[q_sz_loc] = TransInertToCorot_z(gdata, fields.uInertial[2], ix, iy, iz);
}
#endif


#if (USE_ANGULAR_MOMENTUM == TRUE)

void Transformations::TransAngMom2Momen(Data &gdata, gridFunc &gfunc,
                                        ProblemType &Problem)
{
	if(gdata.om[q_sx].getName() != "l_x" || 
	   gdata.om[q_sy].getName() != "l_y" ||
	   gdata.om[q_sz].getName() != "l_z") {
		throw CException(" Has to be angular momentum ");
	}

#if (GEOM == CYLINDRICAL) 
	for(int k = -B; k<=gdata.mx[2]+B; ++k){
		for(int j = -B; j<=gdata.mx[1]+B; ++j){
			for(int i = -B; i<=gdata.mx[0]+B; ++i){
				REAL r_cyl = std::abs(gdata.getCen_x(i));
				gdata.om[q_sy](i,j,k) /= r_cyl;
			}
		}
	}
#else
	throw CException(" Angular Momentum not implemented yet ");
#endif

	gdata.om[q_sx].rename("s_x");
	gdata.om[q_sy].rename("s_y");
	gdata.om[q_sz].rename("s_z");
}
  

void Transformations::TransMomen2AngMom(Data &gdata, gridFunc &gfunc,
                                        ProblemType &Problem)
{
	if(gdata.om[q_sx].getName() != "s_x" || 
	   gdata.om[q_sy].getName() != "s_y" ||
	   gdata.om[q_sz].getName() != "s_z") {
		throw CException(" Has to be momentum ");
	}
  
#if (GEOM == CYLINDRICAL) 
	for(int k = -B; k<=gdata.mx[2]+B; ++k){
		for(int j = -B; j<=gdata.mx[1]+B; ++j){
			for(int i = -B; i<=gdata.mx[0]+B; ++i){
				REAL r_cyl = gdata.getCen_x(i);
				gdata.om[q_sy](i,j,k) *= std::abs(r_cyl);
			}
		}
	}
#else
	throw CException(" Angular Momentum not implemented yet ");
#endif

	gdata.om[q_sx].rename("l_x");
	gdata.om[q_sy].rename("l_y");
	gdata.om[q_sz].rename("l_z");
}
  

void Transformations::TransAngMom2Vel(Data &gdata, gridFunc &gfunc,
                                      ProblemType &Problem)
{
	TransAngMom2Momen(gdata, gfunc, Problem);
	TransMomen2Vel(gdata, gfunc, Problem);
}


void Transformations::TransVel2AngMom(Data &gdata, gridFunc &gfunc,
                                      ProblemType &Problem)
{
	TransVel2Momen(gdata, gfunc, Problem);
	TransMomen2AngMom(gdata, gfunc, Problem);
}
#endif

#if(CRSWITCH_DUAL_ENERGY == CRONOS_ON)
void Transformations::EntropyCorrection(Data &gdata, gridFunc &gfunc,
                                        ProblemType &Problem, int rkstep) {
	// Entropy correction according to Balsara & Spicer (1999) JCP 148, 133


	if(gdata.om[q_sx].getName() != "v_x" || 
	   gdata.om[q_sy].getName() != "v_y" ||
	   gdata.om[q_sz].getName() != "v_z") {
		throw CException(" Must be velocity ");
	}


#if(FLUID_TYPE == CRONOS_MULTIFLUID)
	int n_omIntAll = gdata.fluids->get_N_OMINT_ALL();
#else
	int n_omIntAll = N_OMINT;
#endif
  
	if(gdata.om[n_omIntAll].getName() != "Eges" ||
	   gdata.om[q_Eges].getName() != "Etherm") {
		throw CException(" Need overall and thermal energy for dual energy evolution ");
	}
    
	Pot v_fast(gdata.mx);
	NumMatrix<double,3> switch1(Index::set(0,0,0),
	                            Index::set(gdata.mx[0],gdata.mx[1],gdata.mx[2]));

	Pot & Etherm = gdata.om[q_Eges];
	Pot & Eges = gdata.om[n_omIntAll];


	// Compute characteristic speeds
	for (int k = -1; k <= gdata.mx[2]+1; ++k) {
		for (int j = -1; j <= gdata.mx[1]+1; ++j) {
			for (int i = -1; i <= gdata.mx[0]+1; ++i) {

				REAL rhoinv = 1./gdata.om[q_rho](i,j,k);

				double ptherm = (Problem.gamma - 1.)*Etherm(i,j,k);
				
//#if (FLUID_TYPE == CRONOS_MHD)
				if(magFluid) {
					// REAL Bsq = ((sqr(gdata.om[q_Bx](i,j,k)) +
					//              sqr(gdata.om[q_Bx](i-1,j,k)) +
					//              gdata.om[q_Bx](i,j,k)*gdata.om[q_Bx](i-1,j,k)) +
					//             (sqr(gdata.om[q_By](i,j,k)) +
					//              sqr(gdata.om[q_By](i,j-1,k)) +
					//              gdata.om[q_By](i,j,k)*gdata.om[q_By](i,j-1,k)) +
					//             (sqr(gdata.om[q_Bz](i,j,k)) +
					//              sqr(gdata.om[q_Bz](i,j,k-1)) +
					//              gdata.om[q_Bz](i,j,k)*gdata.om[q_Bz](i,j,k-1)))*onethird;
					REAL Bsq = 2.*get_EMag(gdata, i, j, k);

					v_fast(i,j,k) = sqrt(Problem.gamma*ptherm*rhoinv +
							Bsq*rhoinv);
				} else {
					v_fast(i,j,k) = sqrt(Problem.gamma*ptherm*rhoinv);

				}
			}
		}
	}

	REAL alpha1(0.05), alpha2(0.1), alpha3(0.00125);



	for (int k = 0; k <= gdata.mx[2]; ++k) {
		for (int j = 0; j <= gdata.mx[1]; ++j) {
			for (int i = 0; i <= gdata.mx[0]; ++i) {

				REAL absval = (std::abs(Etherm(i+1,j,k) - Etherm(i-1,j,k)) +
				               std::abs(Etherm(i,j+1,k) - Etherm(i,j-1,k)) +
				               std::abs(Etherm(i,j,k+1) - Etherm(i,j,k-1)));

				REAL minval = min(Etherm(i  ,j,k  ),
				                  Etherm(i+1,j,k), Etherm(i-1,j,k  ),
				                  Etherm(i,j+1,k), Etherm(i,j-1,k  ),
				                  Etherm(i,j,k+1), Etherm(i,j,k-1));

				if(absval < alpha2*minval) {
					switch1(i,j,k) = 1;
				} else {
					switch1(i,j,k) = 0;
				}
			}
		}
	}




	int num = 0;
	for (int k = 0; k <= gdata.mx[2]; ++k) {
		for (int j = 0; j <= gdata.mx[1]; ++j) {
			for (int i = 0; i <= gdata.mx[0]; ++i) {

				bool sw[3] = {false, false, false};
				// Get indicators:
				if(gdata.om[q_Eges](i,j,k) < alpha1*Eges(i,j,k)) {
					sw[0] = true;
				}

				if(switch1(i,j,k) == 1) {
					sw[1] = true;
				}

				// 	REAL  maxval = std::max(v_fast(i,j,k), std::max(v_fast(i+1,j,k), std::max(v_fast(i-1,j,k), std::max(v_fast(i,j+1,k), std::max(v_fast(i,j-1,k), std::max(v_fast(i,j,k+1), v_fast(i,j,k)))))));
				REAL maxval = max(v_fast(i,j  ,k),
				                  v_fast(i+1,j,k), v_fast(i-1,j,k),
				                  v_fast(i,j+1,k), v_fast(i,j-1,k),
				                  v_fast(i,j,k+1), v_fast(i,j,k-1));
#if (GEOM == CARTESIAN) 
				REAL div_v = ((gdata.om[q_sx](i+1,j,k) - gdata.om[q_sx](i-1,j,k))/
				              (gdata.getCen_x(i+1) - gdata.getCen_x(i-1))             +
				              (gdata.om[q_sy](i,j+1,k) - gdata.om[q_sy](i,j-1,k))/
				              (gdata.getCen_y(j+1) - gdata.getCen_y(j-1))             +
				              (gdata.om[q_sz](i,j,k+1) - gdata.om[q_sz](i,j,k-1))/
				              (gdata.getCen_z(k+1) - gdata.getCen_z(k-1)));
				double del = std::min(gdata.getCen_dx(0,i),
				                      std::min(gdata.getCen_dx(1,j),
				                               gdata.getCen_dx(2,k)));

#elif (GEOM == CYLINDRICAL) 
				Pot & v_rad = gdata.om[q_sx];
				Pot & v_phi = gdata.om[q_sy];
				Pot & v_z   = gdata.om[q_sz];
				REAL r_cyl = gdata.getCen_x(i);
				REAL div_v = (((v_rad(i+1,j,k) - v_rad(i-1,j,k))/
				               (gdata.getCen_x(i+1) - gdata.getCen_x(i-1)) +
				               (v_rad(i,j,k)/r_cyl)) +
				              ((v_phi(i,j+1,k) - v_phi(i,j-1,k))/
				               (r_cyl*(gdata.getCen_y(j+1) - gdata.getCen_y(j-1)))) +
				              ((v_z(i,j,k+1) - v_z(i,j,k-1))/
				               (gdata.getCen_z(k+1) - gdata.getCen_z(k-1))));
				double del = std::min(gdata.getCen_dx(0,i),
				                      std::min(r_cyl*gdata.getCen_dx(1,j),
				                               gdata.getCen_dx(2,k)));

#elif (GEOM == SPHERICAL)
				Pot & v_rad = gdata.om[q_sx];
				Pot & v_tht = gdata.om[q_sy];
				Pot & v_phi = gdata.om[q_sz];
				REAL r_sph = gdata.getCen_x(i);
				REAL theta = gdata.getCen_y(j);
				REAL div_v = (((v_rad(i+1,j,k) - v_rad(i-1,j,k))/
				               (gdata.getCen_x(i+1) - gdata.getCen_x(i-1)) +
				               (2.*v_rad(i,j,k)/r_sph)) +
				              (((v_tht(i,j+1,k) - v_tht(i,j-1,k))/
				                (gdata.getCen_y(j+1) - gdata.getCen_y(j-1)) +
				                (v_tht(i,j,k)/tan(theta)))/r_sph) +
				              ((v_phi(i,j,k+1) - v_phi(i,j,k-1))/
				               (r_sph*sin(theta)*(gdata.getCen_z(k+1) -
				                                  gdata.getCen_z(k-1)))));
				double del = std::min(gdata.getCen_dx(0,i),
				                      std::min(r_sph*gdata.getCen_dx(1,j),
				                               r_sph*sin(theta)*gdata.getCen_dx(2,k)));

#endif
				if(-alpha3*maxval < del*div_v) {
					sw[2] = true;
				}

#if(AUX_ENERGY == ENTROPY)
				if((sw[0] && sw[1] && sw[2]) || Etherm(i,j,k) < 1.e-18) {
					num++;
					Etherm(i,j,k) = gdata.om[q_Eadd](i,j,k)*pow(gdata.om[q_rho](i,j,k),Problem.gamma-1.);
					// Compute thermal energy from pressure
//					Etherm(i,j,k) /= Problem.gamma-1.;
				} else {
					gdata.om[q_Eadd](i,j,k) = Etherm(i,j,k)*pow(gdata.om[q_rho](i,j,k),1.-Problem.gamma);
					// Compute pressure from thermal energy
//					gdata.om[q_Eadd](i,j,k) *= Problem.gamma-1.;
				}

				if(Etherm(i,j,k) < 1.e-18) {
					Etherm(i,j,k) = 1.e-18;
					gdata.om[q_Eadd](i,j,k) = Etherm(i,j,k)*pow(gdata.om[q_rho](i,j,k),1.-Problem.gamma);
				}

#else
				if((sw[0] && sw[1] && sw[2]) || Etherm(i,j,k) < 1.e-18) {
					num++;
					Etherm(i,j,k) = gdata.om[q_Eadd](i,j,k);
				} else {
					gdata.om[q_Eadd](i,j,k) = Etherm(i,j,k);
				}
#endif
			}
		}
	}

#ifdef parallel  
	MPI_Barrier(gdata.comm3d);
	gdata.MpiSum(num);
	MPI_Barrier(gdata.comm3d);
#endif
	if(Problem.checkout(6) && num > 1 &&
	   rkstep == TIME_SUBSTEPS-1) {
		if(gdata.rank == 0) {
			cout << " Using Energy correction in " << num << " of ";
#ifdef parallel
			cout << ((gdata.global_mx[0]+1)*
			         (gdata.global_mx[1]+1)*
			         (gdata.global_mx[2]+1));
#else
			cout << (gdata.mx[0]+1)*(gdata.mx[1]+1)*(gdata.mx[2]+1);
#endif
			cout << " cells -- to correct ";
			cout << TempErr << " errors " << endl;
		}
	}
  
	if(num > 0) {
		gfunc.boundary(gdata, Problem, gdata.om[q_Eges],B,q_Eges, iFluid);
	}

}
#endif


#if(CRSWITCH_DUAL_ENERGY == CRONOS_ON)
void Transformations::computeEntropyFromE(Data &gdata, gridFunc &gfunc,
                                          ProblemType &Problem)
{
	if(gdata.om[q_Eges].getName() != "Etherm") {
		throw CException(" Need thermal energy for dual energy computation ");
	}
	// Only if dual energy is chosen
#if(AUX_ENERGY == ENTROPY)
	// Get entropy from thermal energy
	for(int k = -B; k <= gdata.mx[2]+B; ++k) {
		for (int j = -B; j <= gdata.mx[1]+B; ++j) {
			for (int i = -B; i <= gdata.mx[0]+B; ++i) {
				gdata.om[q_Eadd](i,j,k) = 
					gdata.om[q_Eges](i,j,k)*pow(gdata.om[q_rho](i,j,k),1.-Problem.gamma);
			}
		}
	}
#else
	// Save thermal energy
	gdata.om[q_Eadd] = gdata.om[q_Eges];
#endif
}
#endif




REAL Transformations::max(const REAL &val1, const REAL &val2, const REAL &val3,
                          const REAL &val4, const REAL &val5, const REAL &val6,
                          const REAL &val7) {
	return std::max(val1,
	                std::max(val2,
	                         std::max(val3,
	                                  std::max(val4,
	                                           std::max(val5, 
	                                                    std::max(val6,val7))))));
}


REAL Transformations::min(const REAL &val1, const REAL &val2, const REAL &val3,
                          const REAL &val4, const REAL &val5, const REAL &val6,
                          const REAL &val7) {
	return std::min(val1,
	                std::min(val2,
	                         std::min(val3,
	                                  std::min(val4,
	                                           std::min(val5, 
	                                                    std::min(val6,val7))))));
}


void Transformations::get_ConsUser(Data &gdata, phys_fields_1D &fields,
                                   int num) {
	for(int q=0; q<num; ++q) {
		fields.uConORIG[q] = fields.uPriORIG[q];
	}
	
}

inline double Transformations::get_EMag(const Data &gdata,
                                        int ix, int iy, int iz) {
	//! Compute magnetic energy
	REAL E_mag = 0.5*(sqr(0.5*(gdata.om[q_Bx](ix  ,iy  ,iz  ) +
	                           gdata.om[q_Bx](ix-1,iy  ,iz  ))) +
	                  sqr(0.5*(gdata.om[q_By](ix  ,iy  ,iz  ) +
	                           gdata.om[q_By](ix  ,iy-1,iz  ))) +
	                  sqr(0.5*(gdata.om[q_Bz](ix  ,iy  ,iz  ) +
	                           gdata.om[q_Bz](ix  ,iy  ,iz-1))));

	return E_mag;
}

#include "MHD/trafo.C"
#include "HD/trafo_hd.C"

inline REAL Transformations::TransEth2E(REAL rhoinv, REAL psq, REAL Bsq, REAL ETherm) const {
#if (FLUID_TYPE == CRONOS_MHD)
	return TransEth2E_MHD(rhoinv, psq, Bsq, ETherm);
#elif (FLUID_TYPE == CRONOS_HYDRO)
	return TransEth2E_MHD(rhoinv, psq, Bsq, ETherm);
#else
	if(fluidType==CRONOS_MHD) {
		return TransEth2E_MHD(rhoinv, psq, Bsq, ETherm);
	} else {
		return TransEth2E_HD(rhoinv, psq, Bsq, ETherm);
	}
#endif
}

inline REAL Transformations::TransT2E(const ProblemType &Problem, REAL rhoinv, REAL psq, REAL Bsq, REAL ETherm) const {
#if (FLUID_TYPE == CRONOS_MHD)
	return TransT2E_MHD(Problem, rhoinv, psq, Bsq, ETherm);
#elif (FLUID_TYPE == CRONOS_HYDRO)
	return TransT2E_MHD(Problem, rhoinv, psq, Bsq, ETherm);
#else
	if(fluidType==CRONOS_MHD) {
		return TransT2E_MHD(Problem, rhoinv, psq, Bsq, ETherm);
	} else {
		return TransT2E_HD(Problem, rhoinv, psq, Bsq, ETherm);
	}
#endif
}

void Transformations::get_Cons(Queue queue, Data &gdata,
                               ProblemType &Problem,
                               EquationOfState  &eos,
                               cronos::vector<REAL> &Pos,
                               phys_fields_1D &fields,
                               int dir, REAL shift) {
#if (FLUID_TYPE == CRONOS_MHD)
	return get_Cons_MHD(gdata, Problem, eos, Pos, fields, dir, shift);
#elif (FLUID_TYPE == CRONOS_HYDRO)
	return get_Cons_HD(queue, gdata, Problem, eos, Pos, fields, dir, shift);
#else
	if(fluidType==CRONOS_MHD) {
		return get_Cons_MHD(gdata, Problem, eos, Pos, fields, dir, shift);
	} else {
		return get_Cons_HD(gdata, Problem, eos, Pos, fields, dir, shift);
	}
#endif
}

void Transformations::get_Cons(const Data &gdata, const ProblemType &Problem,
	const EquationOfState &eos, phys_fields_0D &fields, int ix, int iy, int iz, int face) {
#if (FLUID_TYPE == CRONOS_MHD)
	return get_Cons_MHD(gdata, Problem, eos, fields, ix, iy, iz, face);
#elif (FLUID_TYPE == CRONOS_HYDRO)
	return get_Cons_HD(gdata, Problem, eos, fields, ix, iy, iz, face);
#else
	if(fluidType==CRONOS_MHD) {
		return get_Cons_MHD(gdata, Problem, eos, fields, ix, iy, iz, face);
	} else {
		return get_Cons_HD(gdata, Problem, eos, fields, ix, iy, iz, face);
	}
#endif
}

//#if (FLUID_TYPE == CRONOS_MHD)
//#elif (FLUID_TYPE == CRONOS_HYDRO)
//#include "HD/trafo_hd.C"
//#endif


//REAL Transformations::TransCorotToInert(Data &gdata, REAL vCorot, int dir, int ix, int iy, int iz) {
//	//! Local transform from corotating to inertial frame
//	/*!
//	 * transformation of x-component
//	 * */
//#if(GEOM == CARTESIAN)
//	if(dir==0) {
//		REAL yPos = gdata.getCen_y(iy);
//		return vCorot - yPos*omegaZ;
//	} else if (dir==1) {
//		REAL xPos = gdata.getCen_x(ix);
//		return vCorot + xPos*omegaZ;
//	} else {
//		return vCorot;
//	}
//#elif(GEOM == CYLINDRICAL)
//	if(dir==1) { // only for v_phi
//
//	}
//#else
//	return vCorot;
//#endif
//
//}
//
//
//REAL Transformations::TransInertToCorot(Data &gdata, REAL vInert, int dir, int ix, int iy, int iz) {
//
//}
//
