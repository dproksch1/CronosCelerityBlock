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


void Transformations::TransPrim2Cons(Data &gdata, gridFunc &gfunc,
                                     ProblemType &Problem, Queue& queue) {
	//! Transform from primitive to conservative variables

	if (ENERGETICS == FULL) {
		if(gdata.om[q_Eges].getName() == "Etherm") {
			TransEth2E(gdata, gfunc, Problem, queue);
		} else if(gdata.om[q_Eges].getName() == "Temp") {
			TransT2E(gdata, gfunc, Problem, queue);
		}
	}

	if(gdata.om[q_sx].getName() == "v_x" && 
	   gdata.om[q_sy].getName() == "v_y" &&
	   gdata.om[q_sz].getName() == "v_z") {

#if (USE_ANGULAR_MOMENTUM == TRUE)
		TransVel2AngMom(gdata, gfunc, Problem);
#else
		TransVel2Momen(gdata, gfunc, Problem, queue);
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

void Transformations::TransVel2Momen(Data &gdata, gridFunc &gfunc,
                                     ProblemType &Problem, Queue &queue) {

	//dproksch TODO: implement tensor multiplication for Celerity buffer
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
					                                sqr(gdata.om[q_sz](i,j,k)))*gdata.om[q_rho](i,j,k);
				}
			}
		}
	} else {
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


void Transformations::TransEth2E(const Data &gdata, gridFunc &gfunc,
                                 ProblemType &Problem, Queue &queue) const
{

	if(gdata.om[q_Eges].getName() != "Etherm") {
		cerr << " Energy is: " << gdata.om[q_Eges].getName() << " " << q_Eges << endl;
		throw CException(" Transformations::TransEth2E - om[q_Eges] is not set as thermal energy ");
	}

	if(Problem.gamma < 1.0000000001) {
		throw CException(" Must not be isothermal ");
	}

	const int izStart = -B + 1;
	const int izEnd = gdata.mx[2] + B;
	const int iyStart = -B + 1;
	const int iyEnd = gdata.mx[1] + B;
	const int ixStart = -B + 1;
	const int ixEnd = gdata.mx[0] + B;

	auto range = Range<3>(izEnd-izStart, iyEnd-iyStart, ixEnd-ixStart);

	if(gdata.om[q_sx].getName() == "v_x" && gdata.om[q_sy].getName() == "v_y" &&
	   gdata.om[q_sz].getName() == "v_z") {

		queue.submit(celerity::allow_by_ref, [=, &gdata](celerity::handler& cgh) {

			celerity::accessor omSYCL_Eges{gdata.omSYCL[q_Eges], cgh, celerity::access::one_to_one{}, celerity::read_write};	
			celerity::accessor omSYCL_sx{gdata.omSYCL[q_sx], cgh, celerity::access::one_to_one{}, celerity::read_only};
			celerity::accessor omSYCL_sy{gdata.omSYCL[q_sy], cgh, celerity::access::one_to_one{}, celerity::read_only};
			celerity::accessor omSYCL_sz{gdata.omSYCL[q_sz], cgh, celerity::access::one_to_one{}, celerity::read_only};
			celerity::accessor omSYCL_rho{gdata.omSYCL[q_rho], cgh, celerity::access::one_to_one{}, celerity::read_only};

			cgh.parallel_for<class BufferInitializationKernel>(range, [=](celerity::item<3> item) {

				size_t iz = item.get_id(0) - izStart;
				size_t iy = item.get_id(1) - iyStart;
				size_t ix = item.get_id(2) - izStart;

				omSYCL_Eges[ix][iy][iz] += 0.5*(sqr(omSYCL_sx[ix][iy][iz]) + 
					                            sqr(omSYCL_sy[ix][iy][iz]) +
					                            sqr(omSYCL_sz[ix][iy][iz])) * omSYCL_rho[ix][iy][iz];
			});
		});

	} else {

		queue.submit(celerity::allow_by_ref, [=, &gdata](celerity::handler& cgh) {

			celerity::accessor omSYCL_Eges{gdata.omSYCL[q_Eges], cgh, celerity::access::one_to_one{}, celerity::read_write};	
			celerity::accessor omSYCL_sx{gdata.omSYCL[q_sx], cgh, celerity::access::one_to_one{}, celerity::read_only};
			celerity::accessor omSYCL_sy{gdata.omSYCL[q_sy], cgh, celerity::access::one_to_one{}, celerity::read_only};
			celerity::accessor omSYCL_sz{gdata.omSYCL[q_sz], cgh, celerity::access::one_to_one{}, celerity::read_only};
			celerity::accessor omSYCL_rho{gdata.omSYCL[q_rho], cgh, celerity::access::one_to_one{}, celerity::read_only};

			cgh.parallel_for<class BufferInitializationKernel>(range, [=](celerity::item<3> item) {

				size_t iz = item.get_id(0) - izStart;
				size_t iy = item.get_id(1) - iyStart;
				size_t ix = item.get_id(2) - izStart;

				omSYCL_Eges[ix][iy][iz] += 0.5*(sqr(omSYCL_sx[ix][iy][iz]) + 
					                            sqr(omSYCL_sy[ix][iy][iz]) +
					                            sqr(omSYCL_sz[ix][iy][iz])) / omSYCL_rho[ix][iy][iz];
			});
		});
	}
	gdata.om[q_Eges].rename("Eges");
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

					double E_kin = 0.5*(sqr(gdata.om[q_sx](i,j,k)) +
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
//	 double fac = 1./((Problem.gamma - 1.)*TNorm);
//	 cout << " factor " << fac << " " << gdata.om[q_Eges](12,12,12) << endl;
//	double fac = kBoverMu_num / (Problem.gamma - 1);

	double fac = kBoverMeanMolWeight_num / (Problem.gamma - 1.);
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


void Transformations::TransT2Eth(const Data &gdata, gridFunc &gfunc,
                                 ProblemType &Problem, Queue &queue) const
{

//	static const double kBoverMu_num = Problem.TrafoNorm->get_num(CRONOS_CONSTANTS::BoltzmannConstant / Problem.meanParticleMass);

	if(gdata.om[q_Eges].getName() != "Temp") {
		throw CException(" om[q_Eges] is not set as temperature ");
	}

	if(Problem.gamma < 1.0000000001) {
		throw CException(" Must not be isothermal ");
	}
//	 double fac = 1./((Problem.gamma - 1.)*TNorm);
//	 cout << " factor " << fac << " " << gdata.om[q_Eges](12,12,12) << endl;
//	double fac = kBoverMu_num / (Problem.gamma - 1);

	double fac = kBoverMeanMolWeight_num / (Problem.gamma - 1.);
	fac /= TNorm;

	const int izStart = -B;
	const int izEnd = gdata.mx[2] + B;
	const int iyStart = -B;
	const int iyEnd = gdata.mx[1] + B;
	const int ixStart = -B;
	const int ixEnd = gdata.mx[0] + B;

	auto range = Range<3>(izEnd-izStart, iyEnd-iyStart, ixEnd-ixStart);

	queue.submit(celerity::allow_by_ref, [=, &gdata](celerity::handler& cgh) {

		celerity::accessor omSYCL_Eges{gdata.omSYCL[q_Eges], cgh, celerity::access::one_to_one{}, celerity::read_write};	
		celerity::accessor omSYCL_rho{gdata.omSYCL[q_rho], cgh, celerity::access::one_to_one{}, celerity::read_only};

		cgh.parallel_for<class BufferInitializationKernel>(range, [=](celerity::item<3> item) {

			size_t iz = item.get_id(0) - izStart;
			size_t iy = item.get_id(1) - iyStart;
			size_t ix = item.get_id(2) - izStart;

			omSYCL_Eges[ix][iy][iz] *= fac * omSYCL_rho[ix][iy][iz];

		});
	});
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


void Transformations::TransT2E(const Data &gdata, gridFunc &gfunc,
                               ProblemType &Problem, Queue &queue) const
{
	if(gdata.om[q_Eges].getName() != "Temp") {
		throw CException(" om[q_Eges] is not set as temperature ");
	}
	TransT2Eth(gdata, gfunc, Problem, queue);
	TransEth2E(gdata, gfunc, Problem, queue);

}

double Transformations::TransEth2T(Data &gdata, gridFunc &gfunc,
                                 ProblemType &Problem) 
{
	if(gdata.om[q_Eges].getName() != "Etherm") {
		throw CException(" om[q_Eges] is not set as thermal energy ");
	}

	if(Problem.gamma < 1.0000000001) {
		throw CException(" Must not be isothermal ");
	}

	double fac = TNorm*(Problem.gamma-1.)/kBoverMeanMolWeight_num;

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

	double TAve(1.);
	if(TPhys) {
		TAve = gdata.computeRMS(q_Eges);
	}
	return TAve;
}


double Transformations::TransE2T(Data &gdata, gridFunc &gfunc,
                               ProblemType &Problem)
{

	if(gdata.om[q_Eges].getName() != "Eges") {
		throw CException(" om[q_Eges] is not set as overall energy ");
	}
	TransE2Eth(gdata, gfunc, Problem);
	return TransEth2T(gdata, gfunc, Problem);

}


void Transformations::get_Cons(const Data &gdata, const ProblemType &Problem,
	const EquationOfState &eos, phys_fields_0D &fields, int ix, int iy, int iz, int face)
{
	NumArray<double> Pos(3); //stays uninitialized

	int dir = face / 2;


	// Compute thermal pressure:
	if(ENERGETICS == FULL) {
		if(thermal) { // Compute thermal pressure from Etherm
			fields.ptherm = (Problem.gamma-1.)*fields.uPri(q_Eges_loc);
		} else { // Compute thermal pressure from Temperature
			fields.ptherm = fields.uPri(q_rho_loc)*fields.uPri(q_Eges_loc);
		}
	} else {
		fields.ptherm = eos.pressure(gdata, Problem, fields.uPri(q_rho_loc),
				Pos(0), Pos(1), Pos(2));
	}

	// Set / compute conservative variables:

	fields.uCon[q_rho_loc] = fields.uPri[q_rho_loc]; // Density

	fields.uCon[q_sx_loc] = fields.uPri[q_sx_loc]*fields.uPri[q_rho_loc]; // v_x -> s_x
	fields.uCon[q_sy_loc] = fields.uPri[q_sy_loc]*fields.uPri[q_rho_loc]; // v_y -> s_y
	fields.uCon[q_sz_loc] = fields.uPri[q_sz_loc]*fields.uPri[q_rho_loc]; // v_z -> s_z

	// When working with the angular momentum: do trafo.
#if (USE_ANGULAR_MOMENTUM == TRUE)

	int ishift = 0;
	if(face % 2 > 0) ishift = 1;
	int shift_vec[3] = {0,0,0};
	shift_vec[dir] = ishift;

	fields.uCon[q_sy_loc] *= gdata.h1(ix, iy, iz,
			shift_vec[0],shift_vec[1],shift_vec[2]);

#endif

	double Bsq(0.);
	// Thermal energy / temperature / overall energy
	if(ENERGETICS == FULL) {
		double rhoinv = 1./fields.uPri[q_rho_loc];

		double psq = (sqr(fields.uPri[q_sx_loc]) + sqr(fields.uPri[q_sy_loc]) +
				sqr(fields.uPri[q_sz_loc]))*sqr(fields.uPri[q_rho_loc]);

		if(thermal) {
			fields.uCon[q_Eges_loc] = TransEth2E(rhoinv,psq,Bsq,
					fields.uPri[q_Eges_loc]);
		} else {
			fields.uCon[q_Eges_loc] = TransT2E(Problem, rhoinv,psq,Bsq,
					fields.uPri[q_Eges_loc]);
		}

#if(CRSWITCH_DUAL_ENERGY == CRONOS_ON)
		// Computation of entropy for dual energy corrections
		double Etherm = fields.ptherm/(Problem.gamma - 1.);
		fields.uPri[q_Eadd_loc] = Etherm*pow(fields.uCon[q_rho_loc],
				1.-Problem.gamma);
		fields.uCon[q_Eadd_loc] = fields.uPri[q_Eadd_loc];
#endif

	}

	// Total pressure:
	fields.ptotal = fields.ptherm;

	// store inertial fram velocity and compute co-rotating frame velocity
#if (USE_COROTATION == CRONOS_ON)
	store_uInert(gdata, fields, ix, iy, iz);
#endif

}


#if (USE_COROTATION == CRONOS_ON)

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
			double yPos = gdata.getCen_y(iy);
			for(int ix = -B; ix<=gdata.mx[0]+B; ++ix){
				// Need different equations for the different grid types.
#if(GEOM == CARTESIAN)
				double xPos = gdata.getCen_x(ix);

				ux(ix, iy, iz) -= yPos*omegaZ;
				uy(ix, iy, iz) += xPos*omegaZ;

#elif(GEOM == CYLINDRICAL)

				double rhoPos = gdata.getCen_x(ix);

				uPhi(ix, iy, iz) += rhoPos*omegaZ;

#else

				double rad = gdata.getCen_x(ix);

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
			double yPos = gdata.getCen_y(iy);
			for(int ix = -B; ix<=gdata.mx[0]+B; ++ix){
				// Need different equations for the different grid types.
#if(GEOM == CARTESIAN)
				double xPos = gdata.getCen_x(ix);

				ux(ix, iy, iz) += yPos*omegaZ;
				uy(ix, iy, iz) -= xPos*omegaZ;

#elif(GEOM == CYLINDRICAL)

				double rhoPos = gdata.getCen_x(ix);

				uPhi(ix, iy, iz) -= rhoPos*omegaZ;

#else

				double rad = gdata.getCen_x(ix);

				uPhi(ix, iy, iz) -= rad*sin(yPos)*omegaZ;

#endif
			}
		}
	}

	gdata.om[q_sx].rename("v_x_Corot");
	gdata.om[q_sy].rename("v_y_Corot");
	gdata.om[q_sz].rename("v_z_Corot");

}

double Transformations::TransCorotToInert_x(Data &gdata, double vCorot, int ix, int iy, int iz) {
	//! Local transform from corotating to inertial frame
	/*!
	 * transformation of x-component
	 * */
#if(GEOM == CARTESIAN)
	double yPos = gdata.getCen_y(iy);
	return vCorot - yPos*omegaZ;
#else
	return vCorot;
#endif

}

double Transformations::TransCorotToInert_y(Data &gdata, double vCorot, int ix, int iy, int iz) {
	//! Local transform from corotating to inertial frame
	/*!
	 * transformation of y-component
	 * */
#if(GEOM == CARTESIAN)
	double xPos = gdata.getCen_x(ix);
	return vCorot + xPos*omegaZ;
#elif(GEOM == CYLINDRICAL)
	double r_cyl = gdata.getCen_x(ix);
	return vCorot + r_cyl*omegaZ;
#else
	return vCorot;
#endif

}

double Transformations::TransCorotToInert_z(Data &gdata, double vCorot, int ix, int iy, int iz) {
	//! Local transform from corotating to inertial frame
	/*!
	 * transformation of z-component
	 * */
#if(GEOM == CARTESIAN)
	return vCorot;
#elif(GEOM == CYLINDRICAL)
	return vCorot;
#else
	double rad = gdata.getCen_x(ix);
	double theta = gdata.getCen_y(iy);
	return vCorot + rad*sin(theta)*omegaZ;
#endif

}



double Transformations::TransInertToCorot_x(Data &gdata, double vInert, int ix, int iy, int iz) {
	//! Local transform from inertial to corotating frame
	/*!
	 * transformation of x-component
	 * */
#if(GEOM == CARTESIAN)
	double yPos = gdata.getCen_y(iy);
	return vInert + yPos*omegaZ;
#else
	return vInert;
#endif

}


double Transformations::TransInertToCorot_y(Data &gdata, double vInert, int ix, int iy, int iz) {
	//! Local transform from inertial to corotating frame
	/*!
	 * transformation of y-component
	 * */
#if(GEOM == CARTESIAN)
	double xPos = gdata.getCen_x(ix);
	return vInert- xPos*omegaZ;
#elif(GEOM == CYLINDRICAL)
	double r_cyl = gdata.getCen_x(ix);
	return vInert - r_cyl*omegaZ;
#else
	return vInert;
#endif

}

double Transformations::TransInertToCorot_z(Data &gdata, double vInert, int ix, int iy, int iz) {
	//! Local transform from inertial to corotating frame
	/*!
	 * transformation of z-component
	 * */
#if(GEOM == CARTESIAN)
	return vInert;
#elif(GEOM == CYLINDRICAL)
	return vInert;
#else
	double rad = gdata.getCen_x(ix);
	double theta = gdata.getCen_y(iy);
	return vInert - rad*sin(theta)*omegaZ;
#endif

}





void Transformations::src_Corotating(Data &gdata, ProblemType &problem, NumMatrix<double, 3> nom[]) {
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
			double theta = gdata.getCen_y(iy);
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
#endif //USE_COROTATION


double Transformations::TransEth2E(double rhoinv, double psq, double Bsq, double ETherm) const {
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

double Transformations::TransT2E(const ProblemType &Problem, double rhoinv, double psq, double Bsq, double ETherm) const {
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

double Transformations::TransEth2E_HD(double rhoinv, double psq, double Bsq, double ETherm)
{
	double Energy(ETherm + 0.5*psq*rhoinv);
	return Energy;
}

double Transformations::TransT2E_HD(ProblemType &Problem,
		double rhoinv, double psq, double Bsq, double Temp)
{
	double Energy(1./(rhoinv*(Problem.gamma-1.))*Temp);
	return TransEth2E(rhoinv, psq, Bsq, Energy);
}

double Transformations::TransEth2E_MHD(double rhoinv, double psq, double Bsq, double ETherm) const
{
	double Energy(ETherm + 0.5*psq*rhoinv + 0.5*Bsq);
	return Energy;
}

double Transformations::TransT2E_MHD(const ProblemType &Problem,
		double rhoinv, double psq, double Bsq, double Temp) const
{
	double Energy(1./(rhoinv*(Problem.gamma-1.))*Temp);
	return TransEth2E(rhoinv, psq, Bsq, Energy);
}