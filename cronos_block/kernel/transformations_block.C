#include "transformations_block.H"
#include "physical_constants.H"
#include <stdlib.h>

Transformations_Block::Transformations_Block(const CronosFluid &fluid, ProblemType &Problem, bool TPhys, int iFluid) {
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

void Transformations_Block::set_thermal(bool thermal) {
	this->thermal = thermal;
}

bool Transformations_Block::get_thermal() {
	return thermal;
}

//! @brief Transform energy density buffer from temperature to overall energy density
void Transformations_Block::TransT2E(Queue &queue, const Data &gdata,
									 GridFunc &gfunc, ProblemType &Problem) const
{
	if(gdata.om[q_Eges].getName() != "Temp") {
		throw CException(" om[q_Eges] is not set as temperature ");
	}
	TransT2Eth(queue, gdata, gfunc, Problem);
	TransEth2E(queue, gdata, gfunc, Problem);

}

//! @brief Transform grid data buffers from primitive to conservative form
void Transformations_Block::TransPrim2Cons(Queue& queue, Data &gdata,
									 GridFunc &gfunc, ProblemType &Problem) {
	//! Transform from primitive to conservative variables

	if (ENERGETICS == FULL) {
		if(gdata.om[q_Eges].getName() == "Etherm") {
			TransEth2E(queue, gdata, gfunc, Problem);
		} else if(gdata.om[q_Eges].getName() == "Temp") {
			TransT2E(queue, gdata, gfunc, Problem);
		}
	}

	if(gdata.om[q_sx].getName() == "v_x" && 
	   gdata.om[q_sy].getName() == "v_y" &&
	   gdata.om[q_sz].getName() == "v_z") {

#if (USE_ANGULAR_MOMENTUM == TRUE)
		//TransVel2AngMom(gdata, gfunc, Problem);
		throw CException("AngMomen not implemented in Transformation");
#else
		TransVel2Momen(queue, gdata, gfunc, Problem);
#endif

	}
}

//! @brief Transform grid data buffers from conservative to primitive form
void Transformations_Block::TransCons2Prim(Queue &queue, Data &gdata,
									 GridFunc &gfunc, ProblemType &Problem) {
	//! Transform from conservative to primitive variables
#if (USE_ANGULAR_MOMENTUM == TRUE)
	TransAngMom2Vel(gdata, gfunc, Problem);
#else
	TransMomen2Vel(queue, gdata, gfunc, Problem);
#endif

	if(ENERGETICS == FULL) {
		if(thermal) {
			TransE2Eth(queue, gdata, gfunc, Problem, 0, false);
		} else {
			//TransE2T(queue, gdata, gfunc, Problem);
		}
	}

}

//! @brief Transform directional buffers from momentum to velocity
void Transformations_Block::TransMomen2Vel(Queue &queue, Data &gdata,
									 GridFunc &gfunc, ProblemType &Problem)
{
	if(gdata.om[q_sx].getName() == "v_x" || 
	   gdata.om[q_sy].getName() == "v_y" ||
	   gdata.om[q_sz].getName() == "v_z") {
		throw CException(" Velocity instead of momentum ");
	}

	auto range = gdata.omSYCL[q_rho].get_range();

	for (int q = q_sx; q <= q_sz; ++q) {
		//pointwise kernel including ghost cells
		queue.submit(celerity::allow_by_ref, [=, &gdata](celerity::handler& cgh) {
			celerity::accessor omSYCL_rho{gdata.omSYCL[q_rho], cgh, celerity::access::one_to_one{}, celerity::read_only};
			celerity::accessor omSYCL_q{gdata.omSYCL[q], cgh, celerity::access::one_to_one{}, celerity::read_write};

			cgh.parallel_for<class Momen2VelTransformationKernel>(range, [=](celerity::item<3> item) {
				omSYCL_q[item.get_id(0)][item.get_id(1)][item.get_id(2)] /= omSYCL_rho[item.get_id(0)][item.get_id(1)][item.get_id(2)];
			});
		});
	}

	gdata.om[q_sx].rename("v_x");
	gdata.om[q_sy].rename("v_y");
	gdata.om[q_sz].rename("v_z");
}

//! @brief Transform directional buffers from velocity to momentum
void Transformations_Block::TransVel2Momen(Queue &queue, Data &gdata,
									 GridFunc &gfunc, ProblemType &Problem) {

	if(gdata.om[q_sx].getName() == "s_x" || 
	   gdata.om[q_sy].getName() == "s_y" ||
	   gdata.om[q_sz].getName() == "s_z") {
		throw CException(" Momentum instead of velocity ");
	}

	auto range = gdata.omSYCL[q_rho].get_range();

	for (int q = q_sx; q <= q_sz; ++q) {
		//pointwise kernel including ghost cells
		queue.submit(celerity::allow_by_ref, [=, &gdata](celerity::handler& cgh) {
			celerity::accessor omSYCL_rho{gdata.omSYCL[q_rho], cgh, celerity::access::one_to_one{}, celerity::read_only};
			celerity::accessor omSYCL_q{gdata.omSYCL[q], cgh, celerity::access::one_to_one{}, celerity::read_write};

			cgh.parallel_for<class Vel2MomenTransformationKernel>(range, [=](celerity::item<3> item) {
				omSYCL_q[item.get_id(0)][item.get_id(1)][item.get_id(2)] *= omSYCL_rho[item.get_id(0)][item.get_id(1)][item.get_id(2)];
			});
		});
	}

	gdata.om[q_sx].rename("s_x");
	gdata.om[q_sy].rename("s_y");
	gdata.om[q_sz].rename("s_z");
}

//! @brief Transform energy density buffer from thermal energy density to overall energy density
void Transformations_Block::TransEth2E(Queue &queue, const Data &gdata,
								 GridFunc &gfunc, ProblemType &Problem) const
{

	if(gdata.om[q_Eges].getName() != "Etherm") {
		cerr << " Energy is: " << gdata.om[q_Eges].getName() << " " << q_Eges << endl;
		throw CException(" Transformations_Block::TransEth2E - om[q_Eges] is not set as thermal energy ");
	}

	if(Problem.gamma < 1.0000000001) {
		throw CException(" Must not be isothermal ");
	}

	auto range = gdata.omSYCL[q_rho].get_range();

	if(gdata.om[q_sx].getName() == "v_x" && gdata.om[q_sy].getName() == "v_y" &&
	   gdata.om[q_sz].getName() == "v_z") {
		
		//pointwise kernel including ghost cells
		//     (non-parallel function starts at -B+1 for some reason)
		queue.submit(celerity::allow_by_ref, [=, &gdata](celerity::handler& cgh) {

			celerity::accessor omSYCL_Eges{gdata.omSYCL[q_Eges], cgh, celerity::access::one_to_one{}, celerity::read_write};	
			celerity::accessor omSYCL_sx{gdata.omSYCL[q_sx], cgh, celerity::access::one_to_one{}, celerity::read_only};
			celerity::accessor omSYCL_sy{gdata.omSYCL[q_sy], cgh, celerity::access::one_to_one{}, celerity::read_only};
			celerity::accessor omSYCL_sz{gdata.omSYCL[q_sz], cgh, celerity::access::one_to_one{}, celerity::read_only};
			celerity::accessor omSYCL_rho{gdata.omSYCL[q_rho], cgh, celerity::access::one_to_one{}, celerity::read_only};

			cgh.parallel_for<class Eth2ETransformationKernelMulVar>(range, [=](celerity::item<3> item) {

				size_t ix = item.get_id(0);
				size_t iy = item.get_id(1);
				size_t iz = item.get_id(2);

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

			cgh.parallel_for<class Eth2ETransformationKernelDivVar>(range, [=](celerity::item<3> item) {

				size_t ix = item.get_id(0);
				size_t iy = item.get_id(1);
				size_t iz = item.get_id(2);

				omSYCL_Eges[ix][iy][iz] += 0.5*(sqr(omSYCL_sx[ix][iy][iz]) + 
					                            sqr(omSYCL_sy[ix][iy][iz]) +
					                            sqr(omSYCL_sz[ix][iy][iz])) / omSYCL_rho[ix][iy][iz];
			});
		});
	}
	gdata.om[q_Eges].rename("Eges");
}

//! @brief Transform energy density buffer from overall energy density to thermal energy density
void Transformations_Block::TransE2Eth(Queue &queue, Data &gdata, GridFunc &gfunc,
								 ProblemType &Problem, int n, bool DualEnergy)
{
	if(Problem.gamma < 1.0000000001) {
		throw CException(" Must not be isothermal ");
	}	
	
	#if(CRSWITCH_DUAL_ENERGY == CRONOS_OFF)
		DualEnergy = false;
	#endif

	int n_omIntAll = N_OMINT;

	// // Saving overall energy:
	// if(ENERGETICS == FULL && DualEnergy) {

	// 	queue.submit(celerity::allow_by_ref, [=, &gdata](celerity::handler& cgh) {
	// 		celerity::accessor om_save_acc{gdata.omSYCL_save[0], cgh, celerity::access::one_to_one{}, celerity::read_write};
	// 		celerity::accessor om_Eges_acc{gdata.omSYCL[4], cgh, celerity::access::one_to_one{}, celerity::read_only};
	// 		cgh.parallel_for<class IntegrationKernel>(gdata.omSYCL_save[0].get_range(), [=](celerity::item<3> item){

	// 			size_t ix = item.get_id(0);
	// 			size_t iy = item.get_id(1);
	// 			size_t iz = item.get_id(2);

	// 			om_save_acc[ix][iy][iz] = om_Eges_acc[ix][iy][iz];
	// 		});
	// 	});

	// 	gdata.om[n_omIntAll].rename(gdata.om[q_Eges].getName());
	// }

	queue.submit(celerity::allow_by_ref, [=, &gdata](celerity::handler& cgh) {
		celerity::accessor om_rho_acc{gdata.omSYCL[0], cgh, celerity::access::one_to_one{}, celerity::read_only};
		celerity::accessor om_sx_acc{gdata.omSYCL[1], cgh, celerity::access::one_to_one{}, celerity::read_only};
		celerity::accessor om_sy_acc{gdata.omSYCL[2], cgh, celerity::access::one_to_one{}, celerity::read_only};
		celerity::accessor om_sz_acc{gdata.omSYCL[3], cgh, celerity::access::one_to_one{}, celerity::read_only};
		celerity::accessor om_Eges_acc{gdata.omSYCL[4], cgh, celerity::access::one_to_one{}, celerity::read_write};
		cgh.parallel_for<class IntegrationKernel>(gdata.omSYCL[0].get_range(), [=](celerity::item<3> item){

			size_t ix = item.get_id(0);
			size_t iy = item.get_id(1);
			size_t iz = item.get_id(2);

			double E_kin = 0.5*(sqr(om_sx_acc[ix][iy][iz]) +
				                  sqr(om_sy_acc[ix][iy][iz]) +
				                  sqr(om_sz_acc[ix][iy][iz]))*om_rho_acc[ix][iy][iz];

			om_Eges_acc[ix][iy][iz] -= E_kin;
		});
	});

	gdata.om[q_Eges].rename("Etherm");
	gfunc.boundary(queue, gdata, Problem, B, q_Eges, iFluid);
}

//! @brief Transform energy density buffer from overall energy density to thermal energy density
// void Transformations_Block::TransE2Eth(Data &gdata, GridFunc &gfunc,
//                                  ProblemType &Problem, Queue &queue) {
	
// 	if(gdata.om[q_Eges].getName() != "Eges") {
// 		throw CException(" om[7] is not set as overall energy ");
// 	}
    
// 	if(Problem.gamma < 1.0000000001) {
// 		throw CException(" Must not be isothermal ");
// 	}

// 	int n_omIntAll = N_OMINT;
// 	TempErr = 0;
// 	bool force_min = Problem.force_min(q_Eges);

// 	if(gdata.om[q_sx].getName() == "v_x" && gdata.om[q_sy].getName() == "v_y" &&
// 	   gdata.om[q_sz].getName() == "v_z") {
		
		
// 		//pointwise kernel including ghost cells
// 		//     (non-parallel function starts at -B+1 for some reason)
// 		queue.submit(celerity::allow_by_ref, [=, &gdata](celerity::handler& cgh) {

// 			celerity::accessor omSYCL_Eges{gdata.omSYCL[q_Eges], cgh, celerity::access::one_to_one{}, celerity::read_write};	
// 			celerity::accessor omSYCL_sx{gdata.omSYCL[q_sx], cgh, celerity::access::one_to_one{}, celerity::read_only};
// 			celerity::accessor omSYCL_sy{gdata.omSYCL[q_sy], cgh, celerity::access::one_to_one{}, celerity::read_only};
// 			celerity::accessor omSYCL_sz{gdata.omSYCL[q_sz], cgh, celerity::access::one_to_one{}, celerity::read_only};
// 			celerity::accessor omSYCL_rho{gdata.omSYCL[q_rho], cgh, celerity::access::one_to_one{}, celerity::read_only};

// 			cgh.parallel_for<class Eth2ETransformationKernelMulVar>(range, [=](celerity::item<3> item) {

// 				size_t ix = item.get_id(0);
// 				size_t iy = item.get_id(1);
// 				size_t iz = item.get_id(2);

// 				omSYCL_Eges[ix][iy][iz] -= 0.5*(sqr(omSYCL_sx[ix][iy][iz]) + 
// 					                            sqr(omSYCL_sy[ix][iy][iz]) +
// 					                            sqr(omSYCL_sz[ix][iy][iz])) * omSYCL_rho[ix][iy][iz];
// 			});
// 		});

// 		celerity::accessor omSYCL_Eges{gdata.omSYCL[q_Eges], cgh, celerity::read_write}

// 		// Check of negative energy only in computational domain - not in ghost cells
// 		for(int k = 0; k<=gdata.mx[2]; ++k){
// 			for(int j = 0; j<=gdata.mx[1]; ++j){
// 				for(int i = 0; i<=gdata.mx[0]; ++i){

	        
// #if(CRSWITCH_DUAL_ENERGY == CRONOS_OFF)
// 					if(!force_min && omSYCL_Eges[i][j][k] < 0.){
// 						cerr << " Transformations_Block::TransE2Eth " << endl;
// //						cerr << q_rho << " " << q_sx << " " << q_Eges << endl;
// 						cerr << " Error! negative energy at cell ("
// 							  << i << " " << j << " " << k << "), pos ("
// 							  << gdata.getCen_x(i) << " "
// 							  << gdata.getCen_y(j) << " "
// 							  << gdata.getCen_z(k) << "); "
// 							  << gdata.om[q_Eges].getName() << " = "
// 							  << omSYCL_Eges[i][j][k] << endl;
// 						TempErr     += 1;
// 						exit(3);
// 					}
// #else
// 					if(!Problem.force_min(q_Eges) && 
// 					   omSYCL_Eges[i][j][k] < 0.){
// 						TempErr += 1;
// //						cout << i << " " << j << " " << k << " " << gdata.om[q_Eges](i,j,k) << endl;
// 					}
// #endif

// 				}
// 			}
// 		}
// 	} else {
// 		// for(int k = -B+1; k<=gdata.mx[2]+B; ++k){
// 		// 	for(int j = -B+1; j<=gdata.mx[1]+B; ++j){
// 		// 		for(int i = -B+1; i<=gdata.mx[0]+B; ++i){

// 		//pointwise kernel EXCLUDING ghost cells

// 		const int izStart = n_ghost[2];
// 		const int izEnd = gdata.mx[2] + n_ghost[2];
// 		const int iyStart = n_ghost[1];
// 		const int iyEnd = gdata.mx[1] + n_ghost[1];
// 		const int ixStart = n_ghost[0];
// 		const int ixEnd = gdata.mx[0] + n_ghost[0];

// 		auto range = CelerityRange<3>(izEnd, iyEnd, ixEnd);

// 		queue.submit(celerity::allow_by_ref, [=, &gdata](celerity::handler& cgh) {

// 			celerity::accessor omSYCL_Eges{gdata.omSYCL[q_Eges], cgh, celerity::access::one_to_one{}, celerity::read_write};	
// 			celerity::accessor omSYCL_sx{gdata.omSYCL[q_sx], cgh, celerity::access::one_to_one{}, celerity::read_only};
// 			celerity::accessor omSYCL_sy{gdata.omSYCL[q_sy], cgh, celerity::access::one_to_one{}, celerity::read_only};
// 			celerity::accessor omSYCL_sz{gdata.omSYCL[q_sz], cgh, celerity::access::one_to_one{}, celerity::read_only};
// 			celerity::accessor omSYCL_rho{gdata.omSYCL[q_rho], cgh, celerity::access::one_to_one{}, celerity::read_only};

// 			cgh.parallel_for<class Eth2ETransformationKernelMulVar>(range, [=](celerity::item<3> item) {

// 				size_t iz = item.get_id(0) + izStart;
// 				size_t iy = item.get_id(1) + iyStart;
// 				size_t ix = item.get_id(2) + ixStart;

// 				omSYCL_Eges[ix][iy][iz] -= 0.5*(sqr(omSYCL_sx[ix][iy][iz]) + 
// 												sqr(omSYCL_sy[ix][iy][iz]) +
// 												sqr(omSYCL_sz[ix][iy][iz])) / omSYCL_rho[ix][iy][iz];

// 				// rework for kernel
// 				if(!force_min &&
// 					gdata.om[q_Eges](i,j,k) < 0){
// 					TempErr     += 1;
// //						cout << i << " " << j << " " << k << " " << gdata.om[q_Eges](i,j,k) << endl;

// 				}
// 			});
// 		});

// 		for(int k = 0; k<=gdata.mx[2]; ++k){
// 			for(int j = 0; j<=gdata.mx[1]; ++j){
// 				for(int i = 0; i<=gdata.mx[0]; ++i){

// 					gdata.om[q_Eges](i,j,k) -= 0.5*(sqr(gdata.om[q_sx](i,j,k)) +    //  sub e_kin
// 					                                sqr(gdata.om[q_sy](i,j,k)) +
// 					                                sqr(gdata.om[q_sz](i,j,k)))/gdata.om[q_rho](i,j,k);

// 					if(!force_min &&
// 					   gdata.om[q_Eges](i,j,k) < 0){
// 						TempErr     += 1;
// //						cout << i << " " << j << " " << k << " " << gdata.om[q_Eges](i,j,k) << endl;

// 					}
// 				}
// 			}
// 		}
// 	}
// }

//! @brief Transform energy density buffer from temperature to thermal energy density
void Transformations_Block::TransT2Eth(Queue &queue, const Data &gdata,
								 GridFunc &gfunc, ProblemType &Problem) const
{

	if(gdata.om[q_Eges].getName() != "Temp") {
		throw CException(" om[q_Eges] is not set as temperature ");
	}

	if(Problem.gamma < 1.0000000001) {
		throw CException(" Must not be isothermal ");
	}

	double fac = kBoverMeanMolWeight_num / (Problem.gamma - 1.);
	fac /= TNorm;

	auto range = gdata.omSYCL[q_rho].get_range();
	
	//pointwise kernel including ghost cells
	queue.submit(celerity::allow_by_ref, [=, &gdata](celerity::handler& cgh) {

		celerity::accessor omSYCL_Eges{gdata.omSYCL[q_Eges], cgh, celerity::access::one_to_one{}, celerity::read_write};	
		celerity::accessor omSYCL_rho{gdata.omSYCL[q_rho], cgh, celerity::access::one_to_one{}, celerity::read_only};

		cgh.parallel_for<class T2EthTransformationKernel>(range, [=](celerity::item<3> item) {

			omSYCL_Eges[item.get_id(0)][item.get_id(1)][item.get_id(2)] *= fac * omSYCL_rho[item.get_id(0)][item.get_id(2)][item.get_id(3)];

		});
	});
	gdata.om[q_Eges].rename("Etherm");
}

