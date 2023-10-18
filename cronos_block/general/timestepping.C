#include <iostream>
#include <chrono>
#include <thread>

#include "timestepping.H"

static bool nom_temp_decl = false;

TimeIntegrator::TimeIntegrator(const int n_saves):
	n_saves(n_saves)
{
	om_save = new NumMatrix<double, DIM>[n_saves];
}


TimeIntegrator::~TimeIntegrator()
{
	delete [] om_save;
	omSYCL_save.clear();
}


void TimeIntegrator::set_corrField(int qch)
{
	this->qch = qch;

}

void TimeIntegrator::set_IntRange(int ibeg[], int iend[])
{
	for(int d=0; d<DIM; ++d) {
		this->ibeg[d] = ibeg[d];
		this->iend[d] = iend[d];
	}

	for (int i=0; i<n_saves; i++) {
		om_save[i].resize(ibeg, iend);
	}
}

void TimeIntegrator::init_omBuffer(Queue &queue, const int mx[])
{
	for (int i = 0; i < 5; i++) {
		omSYCL_save.push_back(CelerityBuffer<double, 3>(CelerityRange<3>(mx[0]+6 +1, mx[1]+6+1, mx[2]+6+1)));
	}

	// auto save_range = omSYCL_save[0].get_range();

	// for (int i = 0; i < 5; i++) {
	// 	queue.submit([=](celerity::handler& cgh) {
	// 		celerity::accessor omSYCL_save_acc{this->omSYCL_save[i], cgh, celerity::access::one_to_one{}, celerity::write_only, celerity::no_init};
	// 		cgh.parallel_for<class BufferInitializationKernel>(save_range, [=](celerity::item<3> item) {
	// 			omSYCL_save_acc[item.get_id(0)][item.get_id(1)][item.get_id(2)] = 0.0;
	// 		});
	// 	});
	// }
}


void TimeIntegrator::save_data(const Pot om[], const int substep)
{
	for (int k = ibeg[2]; k <= iend[2]; ++k){
		for (int j = ibeg[1]; j <= iend[1]; ++j){
			for (int i = ibeg[0]; i <= iend[0]; ++i){
				om_save[substep](i,j,k) = om[qch](i,j,k);
			}
		}
	}
}




double TimeIntegrator::get_dt(Data & gdata, const int substep)
{
	if (substep == TIME_SUBSTEPS - 1)	return gdata.dt;
	else return 0.;
}





RKSteps::RKSteps():
		TimeIntegrator(1)
{
	this->oneThird  = 1./3.;
	this->twoThirds = 2./3.;
}


double RKSteps::get_dt(Data & gdata, const int substep)
{
#if (RK_STEPS == 3)
	if 		(substep == 0) 	return gdata.dt;
	else if (substep == 1) 	return -0.5*gdata.dt;
	else if (substep == 2) 	return 0.5*gdata.dt;
#else
	if 		(substep == 0)	return gdata.dt;
#endif
	else return 0.;
}


void RKSteps::Substep(const Data &gdata, ProblemType &Problem,
                      NumMatrix<double, 3> &nom,
                      Pot om[], const int substep) {

#if  (RK_STEPS == 3)
	if(substep == 0) {
		save_data(om, 0);

		for (int k = ibeg[2]; k <= iend[2]; ++k){
			for (int j = ibeg[1]; j <= iend[1]; ++j){
				for (int i = ibeg[0]; i <= iend[0]; ++i){
					om[qch](i,j,k) -= gdata.dt*nom(i,j,k);
				}
			}
		}
	} else if (substep == 1) {
		for (int k = ibeg[2]; k <= iend[2]; ++k) {
			for (int j = ibeg[1]; j <= iend[1]; ++j) {
				for (int i = ibeg[0]; i <= iend[0]; ++i) {
//					om[qch](i,j,k) = 0.75*om[qold](i,j,k)
					om[qch](i,j,k) = 0.75*om_save[0](i,j,k)
						+0.25*om[qch](i,j,k)
						-0.25*gdata.dt*nom(i,j,k);
				}
			}
		}
	} else {
		double twth   = 2./3.;
		double twthdt = twth*gdata.dt;
		
		for (int k = ibeg[2]; k <= iend[2]; ++k){
			for (int j = ibeg[1]; j <= iend[1]; ++j){
				for (int i = ibeg[0]; i <= iend[0]; ++i){
//					om[qch](i,j,k) = 1./3.*om[qold](i,j,k)
					om[qch](i,j,k) = 1./3.*om_save[0](i,j,k)
						+twth*om[qch](i,j,k)
						-twthdt*nom(i,j,k);
				}
			}
		}
	}
#else

//#if (FLUID_TYPE != CRONOS_MHD)
#if (FLUID_TYPE == CRONOS_HYDRO)
	//	if(Problem.mag == 0) {
	if(om[qch].getName() == "B_x" ||
	   om[qch].getName() == "B_y" ||
	   om[qch].getName() == "B_z" ||
	   om[qch].getName() == "A_x" ||
	   om[qch].getName() == "A_y" ||
	   om[qch].getName() == "A_z") {
		return;
	}
#endif
	//	}

	if (substep == 0) { // First Runge Kutta step
		save_data(om, 0);

		for (int k = ibeg[2]; k <= iend[2]; ++k){
			for (int j = ibeg[1]; j <= iend[1]; ++j){
				for (int i = ibeg[0]; i <= iend[0]; ++i){
					om[qch](i,j,k) -= gdata.dt*nom(i,j,k);
				}
			}
		}
	} else if (substep == 1) { // Second Runge Kutta step
		for (int k = ibeg[2]; k <= iend[2]; ++k) {
			for (int j = ibeg[1]; j <= iend[1]; ++j) {
				for (int i = ibeg[0]; i <= iend[0]; ++i) {
					om[qch](i,j,k) =
						 0.5*om_save[0](i,j,k)
//						 0.5*om[qold](i,j,k)
						+0.5*om[qch](i,j,k)
						-0.5*gdata.dt*nom(i,j,k);
				}
			}
		}
	}
#endif
}

void RKSteps::Substep(Queue &queue, const Data &gdata, CelerityRange<3> omRange,
		CelerityBuffer<nom_t, 3> nomSYCL,
		Pot om[], const int substep, size_t nom_max[3]) {

	double dt = gdata.dt;
	int B = -3;
	int rangeEnd[3] = {gdata.mx[0] + 3 + 1, gdata.mx[1] + 3 + 1, gdata.mx[2] + 3 + 1};
	//auto range = celerity::range<3>(iend[0]-ibeg[0], iend[1]-ibeg[1], iend[2]-ibeg[2]);

#if  (RK_STEPS == 3)
	double nom_temp[nom_max[0]][nom_max[1]][nom_max[2]];
	queue.submit(celerity::allow_by_ref, [=, &nom_temp](celerity::handler& cgh) {
		celerity::accessor nomSYCL_acc{nomSYCL, cgh, celerity::access::all{}, celerity::read_only_host_task};
		cgh.host_task(celerity::on_master_node, [=, &nom_temp]{
			for (int i = -B; i < nom_max[0]+B; i++) {
				for (int j = -B; j < nom_max[1]+B; j++) {
					for (int k = -B; k < nom_max[2]+B; k++) {
						nom_temp[i][j][k] = nomSYCL_acc[i][j*nom_max[2] + k][qch];
					}
				}
			}
		});
	});

	if(substep == 0) {
		save_data(om, 0);
		for (int i = -B; i < nom_max[0]+B; i++) {
			for (int j = -B; j < nom_max[1]+B; j++) {
				for (int k = -B; k < nom_max[2]+B; k++) {
					om[qch](i + B,j + B,k + B) -= dt * nom_temp[i][j][k];
				}
			}
		}
	} else if (substep == 1) {
		for (int i = -B; i < nom_max[0]+B; i++) {
			for (int j = -B; j < nom_max[1]+B; j++) {
				for (int k = -B; k < nom_max[2]+B; k++) {
					om[qch](i + B,j + B,k + B) = 0.75*om_save[0](i + B,j + B,k + B) + 0.25*om[qch](i + B,j + B,k + B)
								 				- dt * nom_temp[i][j][k];
				}
			}
		}
	} else {
		double twth   = 2./3.;
		double twthdt = twth*dt;
		
		for (int i = -B; i < nom_max[0]+B; i++) {
			for (int j = -B; j < nom_max[1]+B; j++) {
				for (int k = -B; k < nom_max[2]+B; k++) {
					om[qch](i + B,j + B,k + B) = 1./3.*om_save[0](i + B,j + B,k + B) + twth*om[qch](i + B,j + B,k + B)
								  - twthdt * nom_temp[i][j][k];
				}
			}
		}
	}
#else

#if (FLUID_TYPE == CRONOS_HYDRO)
	//	if(Problem.mag == 0) {
	if(om[qch].getName() == "B_x" ||
	   om[qch].getName() == "B_y" ||
	   om[qch].getName() == "B_z" ||
	   om[qch].getName() == "A_x" ||
	   om[qch].getName() == "A_y" ||
	   om[qch].getName() == "A_z") {
		return;
	}
#endif
	//	}

	if (substep == 0) { // First Runge Kutta step

		save_data(om, 0);

		queue.submit(celerity::allow_by_ref, [=, &om](celerity::handler& cgh) {
			celerity::accessor nomSYCL_acc{nomSYCL, cgh, celerity::access::all{}, celerity::read_only_host_task};
			cgh.host_task(celerity::on_master_node, [=, &om]{
				for (int i = -B; i < nom_max[0]+B; i++) {
					for (int j = -B; j < nom_max[1]+B; j++) {
						for (int k = -B; k < nom_max[2]+B; k++) {
							om[qch](i + B,j + B,k + B) -= dt * nomSYCL_acc[i][j][k].mat[qch];
						}
					}
				}
			});
		});

	} else if (substep == 1) { // Second Runge Kutta step

		queue.submit(celerity::allow_by_ref, [=, &om](celerity::handler& cgh) {
			celerity::accessor nomSYCL_acc{nomSYCL, cgh, celerity::access::all{}, celerity::read_only_host_task};
			cgh.host_task(celerity::on_master_node, [=, &om]{
				for (int i = -B; i < nom_max[0]+B; i++) {
					for (int j = -B; j < nom_max[1]+B; j++) {
						for (int k = -B; k < nom_max[2]+B; k++) {
							om[qch](i + B,j + B,k + B) = 0.5*om_save[0](i + B,j + B,k + B) + 0.5*om[qch](i + B,j + B,k + B) 
								- 0.5 * dt * nomSYCL_acc[i][j][k].mat[qch];
						}
					}
				}
			});
		});
			
	}
#endif
}

//new one
void RKSteps::Substep(Queue &queue, const Data &gdata, CelerityRange<3> omRange,
		CelerityBuffer<nom_t, 3> nomSYCL,
		double dt, const int substep, size_t nom_max[3]) {

	int B = -3;
	//auto range = celerity::range<3>(iend[0]-ibeg[0], iend[1]-ibeg[1], iend[2]-ibeg[2]);

#if  (RK_STEPS == 3)
	
#else

	//	}

	if (substep == 0) { // First Runge Kutta step

		queue.submit(celerity::allow_by_ref, [=, &gdata](celerity::handler& cgh) {
			celerity::accessor om_rho_acc{gdata.omSYCL[0], cgh, celerity::access::one_to_one{}, celerity::read_only};
			celerity::accessor om_save_rho_acc{this->omSYCL_save[0], cgh, celerity::access::one_to_one{}, celerity::write_only, celerity::no_init};

			cgh.parallel_for<class IntegrationKernel_0>(nomSYCL.get_range(), [=](celerity::item<3> item){

				size_t ix = item.get_id(0);
				size_t iy = item.get_id(1);
				size_t iz = item.get_id(2);

				om_save_rho_acc[ix][iy][iz] = om_rho_acc[ix][iy][iz];
			});
		});

		queue.submit(celerity::allow_by_ref, [=, &gdata](celerity::handler& cgh) {
			celerity::accessor om_sx_acc{gdata.omSYCL[1], cgh, celerity::access::one_to_one{}, celerity::read_only};
			celerity::accessor om_save_sx_acc{this->omSYCL_save[1], cgh, celerity::access::one_to_one{}, celerity::write_only, celerity::no_init};

			cgh.parallel_for<class IntegrationKernel_0>(nomSYCL.get_range(), [=](celerity::item<3> item){

				size_t ix = item.get_id(0);
				size_t iy = item.get_id(1);
				size_t iz = item.get_id(2);

				om_save_sx_acc[ix][iy][iz] = om_sx_acc[ix][iy][iz];
			});
		});

		queue.submit(celerity::allow_by_ref, [=, &gdata](celerity::handler& cgh) {
			celerity::accessor om_sy_acc{gdata.omSYCL[2], cgh, celerity::access::one_to_one{}, celerity::read_only};
			celerity::accessor om_save_sy_acc{this->omSYCL_save[2], cgh, celerity::access::one_to_one{}, celerity::write_only, celerity::no_init};

			cgh.parallel_for<class IntegrationKernel_0>(nomSYCL.get_range(), [=](celerity::item<3> item){

				size_t ix = item.get_id(0);
				size_t iy = item.get_id(1);
				size_t iz = item.get_id(2);

				om_save_sy_acc[ix][iy][iz] = om_sy_acc[ix][iy][iz];
			});
		});

		queue.submit(celerity::allow_by_ref, [=, &gdata](celerity::handler& cgh) {
			celerity::accessor om_sz_acc{gdata.omSYCL[3], cgh, celerity::access::one_to_one{}, celerity::read_only};
			celerity::accessor om_save_sz_acc{this->omSYCL_save[3], cgh, celerity::access::one_to_one{}, celerity::write_only, celerity::no_init};

			cgh.parallel_for<class IntegrationKernel_0>(nomSYCL.get_range(), [=](celerity::item<3> item){

				size_t ix = item.get_id(0);
				size_t iy = item.get_id(1);
				size_t iz = item.get_id(2);

				om_save_sz_acc[ix][iy][iz] = om_sz_acc[ix][iy][iz];
			});
		});

		queue.submit(celerity::allow_by_ref, [=, &gdata](celerity::handler& cgh) {
			celerity::accessor om_Eges_acc{gdata.omSYCL[4], cgh, celerity::access::one_to_one{}, celerity::read_only};
			celerity::accessor om_save_Eges_acc{this->omSYCL_save[4], cgh, celerity::access::one_to_one{}, celerity::write_only, celerity::no_init};

			cgh.parallel_for<class IntegrationKernel_0>(nomSYCL.get_range(), [=](celerity::item<3> item){

				size_t ix = item.get_id(0);
				size_t iy = item.get_id(1);
				size_t iz = item.get_id(2);

				om_save_Eges_acc[ix][iy][iz] = om_Eges_acc[ix][iy][iz];
			});
		});

		queue.submit(celerity::allow_by_ref, [=, &gdata](celerity::handler& cgh) {
			celerity::accessor nomSYCL_acc{nomSYCL, cgh, celerity::access::one_to_one{}, celerity::read_only};
			celerity::accessor om_rho_acc{gdata.omSYCL[0], cgh, celerity::access::one_to_one{}, celerity::read_write};

			cgh.parallel_for<class IntegrationKernel_0>(nomSYCL.get_range(), [=](celerity::item<3> item){

				size_t ix = item.get_id(0);
				size_t iy = item.get_id(1);
				size_t iz = item.get_id(2);

				om_rho_acc[ix][iy][iz] -= dt * nomSYCL_acc[ix][iy][iz].mat[0];
			});
		});

		queue.submit(celerity::allow_by_ref, [=, &gdata](celerity::handler& cgh) {
			celerity::accessor nomSYCL_acc{nomSYCL, cgh, celerity::access::one_to_one{}, celerity::read_only};
			celerity::accessor om_rho_acc{gdata.omSYCL[0], cgh, celerity::access::one_to_one{}, celerity::read_only};
			celerity::accessor om_sx_acc{gdata.omSYCL[1], cgh, celerity::access::one_to_one{}, celerity::read_write};
			cgh.parallel_for<class IntegrationKernel_0>(nomSYCL.get_range(), [=](celerity::item<3> item){

				size_t ix = item.get_id(0);
				size_t iy = item.get_id(1);
				size_t iz = item.get_id(2);

				om_sx_acc[ix][iy][iz] = (om_sx_acc[ix][iy][iz] - dt * nomSYCL_acc[ix][iy][iz].mat[1]) / om_rho_acc[ix][iy][iz];
			});
		});

		queue.submit(celerity::allow_by_ref, [=, &gdata](celerity::handler& cgh) {
			celerity::accessor nomSYCL_acc{nomSYCL, cgh, celerity::access::one_to_one{}, celerity::read_only};
			celerity::accessor om_rho_acc{gdata.omSYCL[0], cgh, celerity::access::one_to_one{}, celerity::read_only};
			celerity::accessor om_sy_acc{gdata.omSYCL[2], cgh, celerity::access::one_to_one{}, celerity::read_write};

			cgh.parallel_for<class IntegrationKernel_0>(nomSYCL.get_range(), [=](celerity::item<3> item){

				size_t ix = item.get_id(0);
				size_t iy = item.get_id(1);
				size_t iz = item.get_id(2);

				om_sy_acc[ix][iy][iz] = (om_sy_acc[ix][iy][iz] - dt * nomSYCL_acc[ix][iy][iz].mat[2]) / om_rho_acc[ix][iy][iz];
			});
		});

		queue.submit(celerity::allow_by_ref, [=, &gdata](celerity::handler& cgh) {
			celerity::accessor nomSYCL_acc{nomSYCL, cgh, celerity::access::one_to_one{}, celerity::read_only};
			celerity::accessor om_rho_acc{gdata.omSYCL[0], cgh, celerity::access::one_to_one{}, celerity::read_only};
			celerity::accessor om_sz_acc{gdata.omSYCL[3], cgh, celerity::access::one_to_one{}, celerity::read_write};

			cgh.parallel_for<class IntegrationKernel_0>(nomSYCL.get_range(), [=](celerity::item<3> item){

				size_t ix = item.get_id(0);
				size_t iy = item.get_id(1);
				size_t iz = item.get_id(2);

				om_sz_acc[ix][iy][iz] = (om_sz_acc[ix][iy][iz] - dt * nomSYCL_acc[ix][iy][iz].mat[3]) / om_rho_acc[ix][iy][iz];
			});
		});

		queue.submit(celerity::allow_by_ref, [=, &gdata](celerity::handler& cgh) {
			celerity::accessor nomSYCL_acc{nomSYCL, cgh, celerity::access::one_to_one{}, celerity::read_write};
			celerity::accessor om_Eges_acc{gdata.omSYCL[4], cgh, celerity::access::one_to_one{}, celerity::read_write};

			cgh.parallel_for<class IntegrationKernel_0>(nomSYCL.get_range(), [=](celerity::item<3> item){

				size_t ix = item.get_id(0);
				size_t iy = item.get_id(1);
				size_t iz = item.get_id(2);

				om_Eges_acc[ix][iy][iz] -= dt * nomSYCL_acc[ix][iy][iz].mat[4];

				for (int d = 0; d < N_OMINT; d++) {
					nomSYCL_acc[ix][iy][iz].mat[d] = 0;
				}
			});
		});

		// queue.submit(celerity::allow_by_ref, [=, &gdata](celerity::handler& cgh) {
		// 	celerity::accessor nomSYCL_acc{nomSYCL, cgh, celerity::access::one_to_one{}, celerity::read_write};
		// 	celerity::accessor om_rho_acc{gdata.omSYCL[0], cgh, celerity::access::one_to_one{}, celerity::read_write};
		// 	celerity::accessor om_sx_acc{gdata.omSYCL[1], cgh, celerity::access::one_to_one{}, celerity::read_write};
		// 	celerity::accessor om_sy_acc{gdata.omSYCL[2], cgh, celerity::access::one_to_one{}, celerity::read_write};
		// 	celerity::accessor om_sz_acc{gdata.omSYCL[3], cgh, celerity::access::one_to_one{}, celerity::read_write};
		// 	celerity::accessor om_Eges_acc{gdata.omSYCL[4], cgh, celerity::access::one_to_one{}, celerity::read_write};
		// 	celerity::accessor om_save_rho_acc{this->omSYCL_save[0], cgh, celerity::access::one_to_one{}, celerity::write_only, celerity::no_init};
		// 	celerity::accessor om_save_sx_acc{this->omSYCL_save[1], cgh, celerity::access::one_to_one{}, celerity::write_only, celerity::no_init};
		// 	celerity::accessor om_save_sy_acc{this->omSYCL_save[2], cgh, celerity::access::one_to_one{}, celerity::write_only, celerity::no_init};
		// 	celerity::accessor om_save_sz_acc{this->omSYCL_save[3], cgh, celerity::access::one_to_one{}, celerity::write_only, celerity::no_init};
		// 	celerity::accessor om_save_Eges_acc{this->omSYCL_save[4], cgh, celerity::access::one_to_one{}, celerity::write_only, celerity::no_init};

		// 	cgh.parallel_for<class IntegrationKernel_0>(nomSYCL.get_range(), [=](celerity::item<3> item){

		// 		size_t ix = item.get_id(0);
		// 		size_t iy = item.get_id(1);
		// 		size_t iz = item.get_id(2);

		// 		om_save_rho_acc[ix][iy][iz] = om_rho_acc[ix][iy][iz];
		// 		om_save_sx_acc[ix][iy][iz] = om_sx_acc[ix][iy][iz];
		// 		om_save_sy_acc[ix][iy][iz] = om_sy_acc[ix][iy][iz];
		// 		om_save_sz_acc[ix][iy][iz] = om_sz_acc[ix][iy][iz];
		// 		om_save_Eges_acc[ix][iy][iz] = om_Eges_acc[ix][iy][iz];

		// 		om_rho_acc[ix][iy][iz] -= dt * nomSYCL_acc[ix][iy][iz].mat[0];
		// 		om_sx_acc[ix][iy][iz] = (om_sx_acc[ix][iy][iz] - dt * nomSYCL_acc[ix][iy][iz].mat[1]) / om_rho_acc[ix][iy][iz];
		// 		om_sy_acc[ix][iy][iz] = (om_sy_acc[ix][iy][iz] - dt * nomSYCL_acc[ix][iy][iz].mat[2]) / om_rho_acc[ix][iy][iz];
		// 		om_sz_acc[ix][iy][iz] = (om_sz_acc[ix][iy][iz] - dt * nomSYCL_acc[ix][iy][iz].mat[3]) / om_rho_acc[ix][iy][iz];
				
		// 		// double E_kin = 0.5*(sqr(om_sx_acc[ix][iy][iz]) +
		// 		// 	                  sqr(om_sy_acc[ix][iy][iz]) +
		// 		// 	                  sqr(om_sz_acc[ix][iy][iz]))*om_rho_acc[ix][iy][iz];				

		// 		om_Eges_acc[ix][iy][iz] -= dt * nomSYCL_acc[ix][iy][iz].mat[4];

		// 		for (int d = 0; d < N_OMINT; d++) {
		// 			nomSYCL_acc[ix][iy][iz].mat[d] = 0;
		// 		}
		// 	});
		// });

	} else if (substep == 1) { // Second Runge Kutta step

		queue.submit(celerity::allow_by_ref, [=, &gdata](celerity::handler& cgh) {
			celerity::accessor nomSYCL_acc{nomSYCL, cgh, celerity::access::one_to_one{}, celerity::read_only};
			celerity::accessor om_rho_acc{gdata.omSYCL[0], cgh, celerity::access::one_to_one{}, celerity::read_write};
			celerity::accessor om_save_rho_acc{this->omSYCL_save[0], cgh, celerity::access::one_to_one{}, celerity::read_only};
			cgh.parallel_for<class IntegrationKernel_1>(nomSYCL.get_range(), [=](celerity::item<3> item){

				size_t ix = item.get_id(0);
				size_t iy = item.get_id(1);
				size_t iz = item.get_id(2);

				om_rho_acc[ix][iy][iz] = 0.5*om_save_rho_acc[ix][iy][iz] + 0.5*om_rho_acc[ix][iy][iz] 
										- 0.5 * dt * nomSYCL_acc[ix][iy][iz].mat[0];
			});
		});

		queue.submit(celerity::allow_by_ref, [=, &gdata](celerity::handler& cgh) {
			celerity::accessor nomSYCL_acc{nomSYCL, cgh, celerity::access::one_to_one{}, celerity::read_only};
			celerity::accessor om_rho_acc{gdata.omSYCL[0], cgh, celerity::access::one_to_one{}, celerity::read_only};
			celerity::accessor om_sx_acc{gdata.omSYCL[1], cgh, celerity::access::one_to_one{}, celerity::read_write};
			celerity::accessor om_save_sx_acc{this->omSYCL_save[1], cgh, celerity::access::one_to_one{}, celerity::read_only};
			cgh.parallel_for<class IntegrationKernel_1>(nomSYCL.get_range(), [=](celerity::item<3> item){

				size_t ix = item.get_id(0);
				size_t iy = item.get_id(1);
				size_t iz = item.get_id(2);

				om_sx_acc[ix][iy][iz] = (0.5*om_save_sx_acc[ix][iy][iz] + 0.5*om_sx_acc[ix][iy][iz] 
										- 0.5 * dt * nomSYCL_acc[ix][iy][iz].mat[1]) / om_rho_acc[ix][iy][iz];

			});
		});

		queue.submit(celerity::allow_by_ref, [=, &gdata](celerity::handler& cgh) {
			celerity::accessor nomSYCL_acc{nomSYCL, cgh, celerity::access::one_to_one{}, celerity::read_only};
			celerity::accessor om_rho_acc{gdata.omSYCL[0], cgh, celerity::access::one_to_one{}, celerity::read_only};
			celerity::accessor om_sy_acc{gdata.omSYCL[2], cgh, celerity::access::one_to_one{}, celerity::read_write};
			celerity::accessor om_save_sy_acc{this->omSYCL_save[2], cgh, celerity::access::one_to_one{}, celerity::read_only};
			cgh.parallel_for<class IntegrationKernel_1>(nomSYCL.get_range(), [=](celerity::item<3> item){

				size_t ix = item.get_id(0);
				size_t iy = item.get_id(1);
				size_t iz = item.get_id(2);

				om_sy_acc[ix][iy][iz] = (0.5*om_save_sy_acc[ix][iy][iz] + 0.5*om_sy_acc[ix][iy][iz] 
										- 0.5 * dt * nomSYCL_acc[ix][iy][iz].mat[2]) / om_rho_acc[ix][iy][iz];

			});
		});

		queue.submit(celerity::allow_by_ref, [=, &gdata](celerity::handler& cgh) {
			celerity::accessor nomSYCL_acc{nomSYCL, cgh, celerity::access::one_to_one{}, celerity::read_only};
			celerity::accessor om_rho_acc{gdata.omSYCL[0], cgh, celerity::access::one_to_one{}, celerity::read_only};
			celerity::accessor om_sz_acc{gdata.omSYCL[3], cgh, celerity::access::one_to_one{}, celerity::read_write};
			celerity::accessor om_save_sz_acc{this->omSYCL_save[3], cgh, celerity::access::one_to_one{}, celerity::read_only};
			cgh.parallel_for<class IntegrationKernel_1>(nomSYCL.get_range(), [=](celerity::item<3> item){

				size_t ix = item.get_id(0);
				size_t iy = item.get_id(1);
				size_t iz = item.get_id(2);

				om_sz_acc[ix][iy][iz] = (0.5*om_save_sz_acc[ix][iy][iz] + 0.5*om_sz_acc[ix][iy][iz] 
										- 0.5 * dt * nomSYCL_acc[ix][iy][iz].mat[3]) / om_rho_acc[ix][iy][iz];

			});
		});

		queue.submit(celerity::allow_by_ref, [=, &gdata](celerity::handler& cgh) {
			celerity::accessor nomSYCL_acc{nomSYCL, cgh, celerity::access::one_to_one{}, celerity::read_write};
			celerity::accessor om_Eges_acc{gdata.omSYCL[4], cgh, celerity::access::one_to_one{}, celerity::read_write};
			celerity::accessor om_save_Eges_acc{this->omSYCL_save[4], cgh, celerity::access::one_to_one{}, celerity::read_only};
			cgh.parallel_for<class IntegrationKernel_1>(nomSYCL.get_range(), [=](celerity::item<3> item){

				size_t ix = item.get_id(0);
				size_t iy = item.get_id(1);
				size_t iz = item.get_id(2);

				om_Eges_acc[ix][iy][iz] = (0.5*om_save_Eges_acc[ix][iy][iz] + 0.5*om_Eges_acc[ix][iy][iz] 
										- 0.5 * dt * nomSYCL_acc[ix][iy][iz].mat[4]);// - dt * nomSYCL_acc[ix][iy][iz].mat[4] - E_kin;

				for (int d = 0; d < N_OMINT; d++) {
					nomSYCL_acc[ix][iy][iz].mat[d] = 0;
				}
			});
		});

		// queue.submit(celerity::allow_by_ref, [=, &gdata](celerity::handler& cgh) {
		// 	celerity::accessor nomSYCL_acc{nomSYCL, cgh, celerity::access::one_to_one{}, celerity::read_write};
		// 	celerity::accessor om_rho_acc{gdata.omSYCL[0], cgh, celerity::access::one_to_one{}, celerity::read_only};
		// 	celerity::accessor om_sx_acc{gdata.omSYCL[1], cgh, celerity::access::one_to_one{}, celerity::read_write};
		// 	celerity::accessor om_sy_acc{gdata.omSYCL[2], cgh, celerity::access::one_to_one{}, celerity::read_write};
		// 	celerity::accessor om_sz_acc{gdata.omSYCL[3], cgh, celerity::access::one_to_one{}, celerity::read_write};
		// 	celerity::accessor om_Eges_acc{gdata.omSYCL[4], cgh, celerity::access::one_to_one{}, celerity::read_write};
		// 	celerity::accessor om_save_rho_acc{this->omSYCL_save[0], cgh, celerity::access::one_to_one{}, celerity::read_only};
		// 	celerity::accessor om_save_sx_acc{this->omSYCL_save[1], cgh, celerity::access::one_to_one{}, celerity::read_only};
		// 	celerity::accessor om_save_sy_acc{this->omSYCL_save[2], cgh, celerity::access::one_to_one{}, celerity::read_only};
		// 	celerity::accessor om_save_sz_acc{this->omSYCL_save[3], cgh, celerity::access::one_to_one{}, celerity::read_only};
		// 	celerity::accessor om_save_Eges_acc{this->omSYCL_save[4], cgh, celerity::access::one_to_one{}, celerity::read_only};
		// 	cgh.parallel_for<class IntegrationKernel_1>(nomSYCL.get_range(), [=](celerity::item<3> item){

		// 		size_t ix = item.get_id(0);
		// 		size_t iy = item.get_id(1);
		// 		size_t iz = item.get_id(2);

		// 		// om_rho_acc[ix][iy][iz] = 0.5*om_save_rho_acc[ix][iy][iz] + 0.5*om_rho_acc[ix][iy][iz] 
		// 		// 						- 0.5 * dt * nomSYCL_acc[ix][iy][iz].mat[0];

		// 		om_sx_acc[ix][iy][iz] = (0.5*om_save_sx_acc[ix][iy][iz] + 0.5*om_sx_acc[ix][iy][iz] 
		// 								- 0.5 * dt * nomSYCL_acc[ix][iy][iz].mat[1]) / om_rho_acc[ix][iy][iz];

		// 		om_sy_acc[ix][iy][iz] = (0.5*om_save_sy_acc[ix][iy][iz] + 0.5*om_sy_acc[ix][iy][iz] 
		// 								- 0.5 * dt * nomSYCL_acc[ix][iy][iz].mat[2]) / om_rho_acc[ix][iy][iz];

		// 		om_sz_acc[ix][iy][iz] = (0.5*om_save_sz_acc[ix][iy][iz] + 0.5*om_sz_acc[ix][iy][iz] 
		// 								- 0.5 * dt * nomSYCL_acc[ix][iy][iz].mat[3]) / om_rho_acc[ix][iy][iz];

		// 		// double E_kin = 0.5*(sqr(om_sx_acc[ix][iy][iz]) +
		// 		// 	                  sqr(om_sy_acc[ix][iy][iz]) +
		// 		// 	                  sqr(om_sz_acc[ix][iy][iz]))*om_rho_acc[ix][iy][iz];

		// 		om_Eges_acc[ix][iy][iz] = (0.5*om_save_Eges_acc[ix][iy][iz] + 0.5*om_Eges_acc[ix][iy][iz] 
		// 								- 0.5 * dt * nomSYCL_acc[ix][iy][iz].mat[4]);// - dt * nomSYCL_acc[ix][iy][iz].mat[4] - E_kin;

		// 		for (int d = 0; d < N_OMINT; d++) {
		// 			nomSYCL_acc[ix][iy][iz].mat[d] = 0;
		// 		}
		// 	});
		// });
			
	}

	gdata.om[1].rename("v_x");
	gdata.om[2].rename("v_y");
	gdata.om[3].rename("v_z");
#endif
}

void RKSteps::Substep(Queue &queue, const Data &gdata, CelerityRange<3> omRange,
		CelerityBuffer<nom_t, 3> nomSYCL, int q,
		double dt, const int substep, size_t nom_max[3]) {

	int B = -3;
	//auto range = celerity::range<3>(iend[0]-ibeg[0], iend[1]-ibeg[1], iend[2]-ibeg[2]);

#if  (RK_STEPS == 3)
	
#else

	//	}
cout << "q == " << q << endl;
	if (substep == 0) { // First Runge Kutta step

		if (q == 0) {cout << "s0q0" << endl;
			queue.submit(celerity::allow_by_ref, [=, &gdata](celerity::handler& cgh) {
				celerity::accessor nomSYCL_acc{nomSYCL, cgh, celerity::access::one_to_one{}, celerity::read_only};
				celerity::accessor om_rho_out_acc{gdata.omSYCL_out[0], cgh, celerity::access::one_to_one{}, celerity::read_write};
				celerity::accessor om_rho_acc{gdata.omSYCL[0], cgh, celerity::access::one_to_one{}, celerity::read_write};
				celerity::accessor om_save_rho_acc{this->omSYCL_save[0], cgh, celerity::access::one_to_one{}, celerity::write_only, celerity::no_init};

				cgh.parallel_for<class IntegrationKernel_0_0>(nomSYCL.get_range(), [=](celerity::item<3> item){

					size_t ix = item.get_id(0);
					size_t iy = item.get_id(1);
					size_t iz = item.get_id(2);

					om_save_rho_acc[ix][iy][iz] = om_rho_acc[ix][iy][iz];

					om_rho_out_acc[ix][iy][iz] = om_rho_acc[ix][iy][iz] - dt * nomSYCL_acc[ix][iy][iz].mat[0];

				});
			});

			gdata.om[1].rename("v_x");
			gdata.om[2].rename("v_y");
			gdata.om[3].rename("v_z");
		} else if (q == 1) {cout << "s0q1" << endl;
			queue.submit(celerity::allow_by_ref, [=, &gdata](celerity::handler& cgh) {
				celerity::accessor nomSYCL_acc{nomSYCL, cgh, celerity::access::one_to_one{}, celerity::read_only};
				celerity::accessor om_rho_acc{gdata.omSYCL_out[0], cgh, celerity::access::one_to_one{}, celerity::read_only};
				celerity::accessor om_sx_out_acc{gdata.omSYCL_out[1], cgh, celerity::access::one_to_one{}, celerity::read_write};
				celerity::accessor om_sx_acc{gdata.omSYCL[1], cgh, celerity::access::one_to_one{}, celerity::read_write};
				celerity::accessor om_save_sx_acc{this->omSYCL_save[1], cgh, celerity::access::one_to_one{}, celerity::write_only, celerity::no_init};
				
				cgh.parallel_for<class IntegrationKernel_0_1>(nomSYCL.get_range(), [=](celerity::item<3> item){

					size_t ix = item.get_id(0);
					size_t iy = item.get_id(1);
					size_t iz = item.get_id(2);

					om_save_sx_acc[ix][iy][iz] = om_sx_acc[ix][iy][iz];

					double rho = om_rho_acc[ix][iy][iz] - dt * nomSYCL_acc[ix][iy][iz].mat[0];
					om_sx_out_acc[ix][iy][iz] = (om_sx_acc[ix][iy][iz] - dt * nomSYCL_acc[ix][iy][iz].mat[1]) / 
												om_rho_acc[ix][iy][iz];
				
				});
			});
		} else if (q == 2) {cout << "s0q2" << endl;
			queue.submit(celerity::allow_by_ref, [=, &gdata](celerity::handler& cgh) {
				celerity::accessor nomSYCL_acc{nomSYCL, cgh, celerity::access::one_to_one{}, celerity::read_only};
				celerity::accessor om_rho_acc{gdata.omSYCL_out[0], cgh, celerity::access::one_to_one{}, celerity::read_only};
				celerity::accessor om_sy_out_acc{gdata.omSYCL_out[2], cgh, celerity::access::one_to_one{}, celerity::read_write};
				celerity::accessor om_sy_acc{gdata.omSYCL[2], cgh, celerity::access::one_to_one{}, celerity::read_write};
				celerity::accessor om_save_sy_acc{this->omSYCL_save[2], cgh, celerity::access::one_to_one{}, celerity::write_only, celerity::no_init};
				
				cgh.parallel_for<class IntegrationKernel_0_2>(nomSYCL.get_range(), [=](celerity::item<3> item){

					size_t ix = item.get_id(0);
					size_t iy = item.get_id(1);
					size_t iz = item.get_id(2);

					om_save_sy_acc[ix][iy][iz] = om_sy_acc[ix][iy][iz];

					om_sy_out_acc[ix][iy][iz] = (om_sy_acc[ix][iy][iz] - dt * nomSYCL_acc[ix][iy][iz].mat[2]) / 
												om_rho_acc[ix][iy][iz];
				
				});
			});
		} else if (q == 3) {cout << "s0q3" << endl;
			queue.submit(celerity::allow_by_ref, [=, &gdata](celerity::handler& cgh) {
				celerity::accessor nomSYCL_acc{nomSYCL, cgh, celerity::access::one_to_one{}, celerity::read_only};
				celerity::accessor om_rho_acc{gdata.omSYCL_out[0], cgh, celerity::access::one_to_one{}, celerity::read_only};
				celerity::accessor om_sz_out_acc{gdata.omSYCL_out[3], cgh, celerity::access::one_to_one{}, celerity::read_write};
				celerity::accessor om_sz_acc{gdata.omSYCL[3], cgh, celerity::access::one_to_one{}, celerity::read_write};
				celerity::accessor om_save_sz_acc{this->omSYCL_save[3], cgh, celerity::access::one_to_one{}, celerity::write_only, celerity::no_init};
				
				cgh.parallel_for<class IntegrationKernel_0_3>(nomSYCL.get_range(), [=](celerity::item<3> item){

					size_t ix = item.get_id(0);
					size_t iy = item.get_id(1);
					size_t iz = item.get_id(2);

					om_save_sz_acc[ix][iy][iz] = om_sz_acc[ix][iy][iz];

					om_sz_out_acc[ix][iy][iz] = (om_sz_acc[ix][iy][iz] - dt * nomSYCL_acc[ix][iy][iz].mat[3]) / 
												om_rho_acc[ix][iy][iz];
				
				});
			});
		} else if (q == 4) {cout << "s0q4" << endl;
			queue.submit(celerity::allow_by_ref, [=, &gdata](celerity::handler& cgh) {
				celerity::accessor nomSYCL_acc{nomSYCL, cgh, celerity::access::one_to_one{}, celerity::read_only};
				celerity::accessor om_Eges_out_acc{gdata.omSYCL_out[4], cgh, celerity::access::one_to_one{}, celerity::read_write};
				celerity::accessor om_Eges_acc{gdata.omSYCL[4], cgh, celerity::access::one_to_one{}, celerity::read_write};
				celerity::accessor om_save_Eges_acc{this->omSYCL_save[4], cgh, celerity::access::one_to_one{}, celerity::write_only};

				cgh.parallel_for<class IntegrationKernel_0_4>(nomSYCL.get_range(), [=](celerity::item<3> item){

					size_t ix = item.get_id(0);
					size_t iy = item.get_id(1);
					size_t iz = item.get_id(2);

					om_save_Eges_acc[ix][iy][iz] = om_Eges_acc[ix][iy][iz];
					om_Eges_out_acc[ix][iy][iz] = om_Eges_acc[ix][iy][iz] - dt * nomSYCL_acc[ix][iy][iz].mat[4];

				});
			});
		}

	} else if (substep == 1) { // Second Runge Kutta step

		if (q == 0) {cout << "s1q0" << endl;
			queue.submit(celerity::allow_by_ref, [=, &gdata](celerity::handler& cgh) {
				celerity::accessor nomSYCL_acc{nomSYCL, cgh, celerity::access::one_to_one{}, celerity::read_only};
				celerity::accessor om_rho_acc{gdata.omSYCL[0], cgh, celerity::access::one_to_one{}, celerity::read_write};
				celerity::accessor om_save_rho_acc{this->omSYCL_save[0], cgh, celerity::access::one_to_one{}, celerity::read_only};
				cgh.parallel_for<class IntegrationKernel_1_0>(nomSYCL.get_range(), [=](celerity::item<3> item){

					size_t ix = item.get_id(0);
					size_t iy = item.get_id(1);
					size_t iz = item.get_id(2);

					om_rho_acc[ix][iy][iz] = 0.5*om_save_rho_acc[ix][iy][iz] + 0.5*om_rho_acc[ix][iy][iz] 
											- 0.5 * dt * nomSYCL_acc[ix][iy][iz].mat[0];
				});
			});

			gdata.om[1].rename("v_x");
			gdata.om[2].rename("v_y");
			gdata.om[3].rename("v_z");
		} else if (q == 1) {cout << "s1q1" << endl;
			queue.submit(celerity::allow_by_ref, [=, &gdata](celerity::handler& cgh) {
				celerity::accessor nomSYCL_acc{nomSYCL, cgh, celerity::access::one_to_one{}, celerity::read_only};
				celerity::accessor om_rho_acc{gdata.omSYCL[0], cgh, celerity::access::one_to_one{}, celerity::read_only};
				celerity::accessor om_sx_acc{gdata.omSYCL[1], cgh, celerity::access::one_to_one{}, celerity::read_write};
				celerity::accessor om_save_rho_acc{this->omSYCL_save[0], cgh, celerity::access::one_to_one{}, celerity::read_only};
				celerity::accessor om_save_sx_acc{this->omSYCL_save[1], cgh, celerity::access::one_to_one{}, celerity::read_only};
				cgh.parallel_for<class IntegrationKernel_1_1>(nomSYCL.get_range(), [=](celerity::item<3> item){

					size_t ix = item.get_id(0);
					size_t iy = item.get_id(1);
					size_t iz = item.get_id(2);

					double rho = 0.5*om_save_rho_acc[ix][iy][iz] + 0.5*om_rho_acc[ix][iy][iz] 
											- 0.5 * dt * nomSYCL_acc[ix][iy][iz].mat[0];

					om_sx_acc[ix][iy][iz] = (0.5*om_save_sx_acc[ix][iy][iz] + 0.5*om_sx_acc[ix][iy][iz] 
											- 0.5 * dt * nomSYCL_acc[ix][iy][iz].mat[1]) / om_rho_acc[ix][iy][iz];

				});
			});
		} else if (q == 2) {cout << "s1q2" << endl;
			queue.submit(celerity::allow_by_ref, [=, &gdata](celerity::handler& cgh) {
				celerity::accessor nomSYCL_acc{nomSYCL, cgh, celerity::access::one_to_one{}, celerity::read_only};
				celerity::accessor om_rho_acc{gdata.omSYCL[0], cgh, celerity::access::one_to_one{}, celerity::read_only};
				celerity::accessor om_sy_acc{gdata.omSYCL[2], cgh, celerity::access::one_to_one{}, celerity::read_write};
				celerity::accessor om_save_rho_acc{this->omSYCL_save[0], cgh, celerity::access::one_to_one{}, celerity::read_only};
				celerity::accessor om_save_sy_acc{this->omSYCL_save[2], cgh, celerity::access::one_to_one{}, celerity::read_only};
				cgh.parallel_for<class IntegrationKernel_1_2>(nomSYCL.get_range(), [=](celerity::item<3> item){

					size_t ix = item.get_id(0);
					size_t iy = item.get_id(1);
					size_t iz = item.get_id(2);

					double rho = 0.5*om_save_rho_acc[ix][iy][iz] + 0.5*om_rho_acc[ix][iy][iz] 
											- 0.5 * dt * nomSYCL_acc[ix][iy][iz].mat[0];

					om_sy_acc[ix][iy][iz] = (0.5*om_save_sy_acc[ix][iy][iz] + 0.5*om_sy_acc[ix][iy][iz] 
											- 0.5 * dt * nomSYCL_acc[ix][iy][iz].mat[2]) / om_rho_acc[ix][iy][iz];

				});
			});
		} else if (q == 3) {cout << "s1q3" << endl;
			queue.submit(celerity::allow_by_ref, [=, &gdata](celerity::handler& cgh) {
				celerity::accessor nomSYCL_acc{nomSYCL, cgh, celerity::access::one_to_one{}, celerity::read_only};
				celerity::accessor om_rho_acc{gdata.omSYCL[0], cgh, celerity::access::one_to_one{}, celerity::read_only};
				celerity::accessor om_sz_acc{gdata.omSYCL[3], cgh, celerity::access::one_to_one{}, celerity::read_write};
				celerity::accessor om_save_rho_acc{this->omSYCL_save[0], cgh, celerity::access::one_to_one{}, celerity::read_only};
				celerity::accessor om_save_sz_acc{this->omSYCL_save[3], cgh, celerity::access::one_to_one{}, celerity::read_only};
				cgh.parallel_for<class IntegrationKernel_1_3>(nomSYCL.get_range(), [=](celerity::item<3> item){

					size_t ix = item.get_id(0);
					size_t iy = item.get_id(1);
					size_t iz = item.get_id(2);

					double rho = 0.5*om_save_rho_acc[ix][iy][iz] + 0.5*om_rho_acc[ix][iy][iz] 
											- 0.5 * dt * nomSYCL_acc[ix][iy][iz].mat[0];

					om_sz_acc[ix][iy][iz] = (0.5*om_save_sz_acc[ix][iy][iz] + 0.5*om_sz_acc[ix][iy][iz] 
											- 0.5 * dt * nomSYCL_acc[ix][iy][iz].mat[3]) / om_rho_acc[ix][iy][iz];

				});
			});
		} else if (q == 4) {cout << "s1q4" << endl;
			queue.submit(celerity::allow_by_ref, [=, &gdata](celerity::handler& cgh) {
				celerity::accessor nomSYCL_acc{nomSYCL, cgh, celerity::access::one_to_one{}, celerity::read_only};
				celerity::accessor om_Eges_acc{gdata.omSYCL[4], cgh, celerity::access::one_to_one{}, celerity::read_write};
				celerity::accessor om_save_Eges_acc{this->omSYCL_save[4], cgh, celerity::access::one_to_one{}, celerity::read_only};
				cgh.parallel_for<class IntegrationKernel_1_4>(nomSYCL.get_range(), [=](celerity::item<3> item){

					size_t ix = item.get_id(0);
					size_t iy = item.get_id(1);
					size_t iz = item.get_id(2);

					om_Eges_acc[ix][iy][iz] = (0.5*om_save_Eges_acc[ix][iy][iz] + 0.5*om_Eges_acc[ix][iy][iz] 
											- 0.5 * dt * nomSYCL_acc[ix][iy][iz].mat[4]);// - dt * nomSYCL_acc[ix][iy][iz].mat[4] - E_kin;

				});
			});
		}
			
	}

#endif
}


VanLeerIntegrator::VanLeerIntegrator():
		TimeIntegrator(1)
{
}


double VanLeerIntegrator::get_dt(Data & gdata, const int substep)
{
	if 		(substep == 0)	return 0.5 * gdata.dt;
	else if	(substep == 1)	return 0.5 * gdata.dt;
	else return 0.;

}

void VanLeerIntegrator::Substep(const Data &gdata, ProblemType &Problem,
                      NumMatrix<double, 3> &nom,
                      Pot om[], const int substep) {
	if (substep == 0) {
		save_data(om, 0);
		for (int k = ibeg[2]; k <= iend[2]; ++k){
			for (int j = ibeg[1]; j <= iend[1]; ++j){
				for (int i = ibeg[0]; i <= iend[0]; ++i){
					om[qch](i,j,k) -= 0.5*gdata.dt*nom(i,j,k);
				}
			}
		}

	} else if (substep == 1) {
		for (int k = ibeg[2]; k <= iend[2]; ++k) {
			for (int j = ibeg[1]; j <= iend[1]; ++j) {
				for (int i = ibeg[0]; i <= iend[0]; ++i) {
					om[qch](i,j,k) =
//							om[qold](i,j,k)
							om_save[0](i,j,k)
							- gdata.dt*nom(i,j,k);
				}
			}
		}
	}
}

void VanLeerIntegrator::Substep(Queue &queue, const Data &gdata, CelerityRange<3> omRange,
                      CelerityBuffer<nom_t, 3> nomSYCL,
                      Pot om[], const int substep, size_t nom_max[3]) {}

void VanLeerIntegrator::Substep(Queue &queue, const Data &gdata, CelerityRange<3> omRange,
                      CelerityBuffer<nom_t, 3> nomSYCL,
                      double dt, const int substep, size_t nom_max[3]) {}

void VanLeerIntegrator::Substep(Queue &queue, const Data &gdata, CelerityRange<3> omRange,
		CelerityBuffer<nom_t, 3> nomSYCL, int q,
		double dt, const int substep, size_t nom_max[3]) {}