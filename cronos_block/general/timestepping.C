#include <iostream>
#include <chrono>
#include <thread>

#include "timestepping.H"
#include "utils.H"

TimeIntegrator::TimeIntegrator(const int n_saves):
	n_saves(n_saves)
{
	om_save = new NumMatrix<double, DIM>[n_saves];
}


TimeIntegrator::~TimeIntegrator()
{
	delete [] om_save;
}


void TimeIntegrator::set_corrField(int qch)
{
	this->qch = qch;

}

void TimeIntegrator::init_omSYCL_save(const int mx[DIM])
{
	omSYCL_save.clear();
	omSYCL_save.push_back(CelerityBuffer<double, 3>(Range<3>(mx[0]+6 +1, mx[1]+6+1, mx[2]+6+1)));
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

void TimeIntegrator::save_data(Queue &queue, std::vector<CelerityBuffer<double, 3>> omSYCL, const int substep, const int B)
{
	if (qch == 0) {
		queue.submit(celerity::allow_by_ref, [=](celerity::handler& cgh) {
			celerity::accessor omSYCL_acc{omSYCL[qch], cgh, celerity::access::all{}, celerity::read_only_host_task};
			celerity::accessor omSYCL_save_acc{omSYCL_save[substep], cgh, celerity::access::all{}, celerity::write_only_host_task};
			cgh.host_task(celerity::on_master_node, [=]{
				for (int k = ibeg[2]-B; k <= iend[2]-B; ++k){
					for (int j = ibeg[1]-B; j <= iend[1]-B; ++j){
						for (int i = ibeg[0]-B; i <= iend[0]-B; ++i){
							omSYCL_save_acc[i][j][k] = omSYCL_acc[i][j][k];
						}
					}
				}
			});
		});
	} else if (qch > 0 && qch < 4) {
		queue.submit(celerity::allow_by_ref, [=](celerity::handler& cgh) {
			celerity::accessor omSYCL_acc{omSYCL[qch], cgh, celerity::access::all{}, celerity::read_only_host_task};
			celerity::accessor omSYCL_rho_acc{omSYCL[gpu::fluidConst_q_rho], cgh, celerity::access::all{}, celerity::read_only_host_task};
			celerity::accessor omSYCL_save_acc{omSYCL_save[substep], cgh, celerity::access::all{}, celerity::write_only_host_task};
			cgh.host_task(celerity::on_master_node, [=]{
				for (int k = ibeg[2]-B; k <= iend[2]-B; ++k){
					for (int j = ibeg[1]-B; j <= iend[1]-B; ++j){
						for (int i = ibeg[0]-B; i <= iend[0]-B; ++i){
							omSYCL_save_acc[i][j][k] = omSYCL_acc[i][j][k] * omSYCL_rho_acc[i][j][k];
						}
					}
				}
			});
		});
	} else {
		queue.submit(celerity::allow_by_ref, [=](celerity::handler& cgh) {
			celerity::accessor omSYCL_acc{omSYCL[qch], cgh, celerity::access::all{}, celerity::read_only_host_task};
			celerity::accessor omSYCL_rho_acc{omSYCL[gpu::fluidConst_q_rho], cgh, celerity::access::all{}, celerity::read_only_host_task};
			celerity::accessor omSYCL_sx_acc{omSYCL[gpu::fluidConst_q_sx], cgh, celerity::access::all{}, celerity::read_only_host_task};
			celerity::accessor omSYCL_sy_acc{omSYCL[gpu::fluidConst_q_sy], cgh, celerity::access::all{}, celerity::read_only_host_task};
			celerity::accessor omSYCL_sz_acc{omSYCL[gpu::fluidConst_q_sz], cgh, celerity::access::all{}, celerity::read_only_host_task};
			celerity::accessor omSYCL_save_acc{omSYCL_save[substep], cgh, celerity::access::all{}, celerity::write_only_host_task};
			cgh.host_task(celerity::on_master_node, [=]{
				for (int k = ibeg[2]-B; k <= iend[2]-B; ++k){
					for (int j = ibeg[1]-B; j <= iend[1]-B; ++j){
						for (int i = ibeg[0]-B; i <= iend[0]-B; ++i){
							if (ENERGETICS == FULL) {
								double Ekin = 0.5*(sqr(omSYCL_sx_acc[i][j][k]) + 
														sqr(omSYCL_sy_acc[i][j][k]) +
														sqr(omSYCL_sz_acc[i][j][k])) * omSYCL_rho_acc[i][j][k];
								omSYCL_save_acc[i][j][k] = omSYCL_acc[i][j][k] + Ekin;
							}
						}
					}
				}
			});
		});
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
		CelerityBuffer<double, 3> nomSYCL, 
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
	double nom_temp[nom_max[0]][nom_max[1]][nom_max[2]];
	{
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
	});}

	queue.slow_full_sync();

	if (substep == 0) { // First Runge Kutta step

		save_data(om, 0);
		
		for (int i = -B; i < nom_max[0]+B; i++) {
			for (int j = -B; j < nom_max[1]+B; j++) {
				for (int k = -B; k < nom_max[2]+B; k++) {
					if(i == 1 + B && j == 1 + B && k == 1 + B) {
						
					}
					om[qch](i + B,j + B,k + B) -= dt * nom_temp[i][j][k];
				}
			}
		}

	} else if (substep == 1) { // Second Runge Kutta step

		for (int i = -B; i < nom_max[0]+B; i++) {
			for (int j = -B; j < nom_max[1]+B; j++) {
				for (int k = -B; k < nom_max[2]+B; k++) {
					om[qch](i + B,j + B,k + B) = 0.5*om_save[0](i + B,j + B,k + B) + 0.5*om[qch](i + B,j + B,k + B) 
								- 0.5 * dt * nom_temp[i][j][k];
				}
			}
		}
	}
#endif
}

void RKSteps::Substep(Queue &queue, const Data &gdata, CelerityRange<3> omRange,
		CelerityBuffer<double, 3> nomSYCL, Pot om[], const int substep, size_t nom_max1, size_t nom_max[3]) {

	double dt = gdata.dt;
	int B = -3;
	int rangeEnd[3] = {gdata.mx[0] + 3 + 1, gdata.mx[1] + 3 + 1, gdata.mx[2] + 3 + 1};
	int qch = this->qch;
	size_t nom_max0 = nom_max[0];
	size_t nom_max2 = nom_max[2];

#if  (RK_STEPS == 3)
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

	if (substep == 0) { // First Runge Kutta step
			save_data(queue, gdata.omSYCL, 0, B);
	}

	bool isVelocity = (gdata.om[gpu::fluidConst_q_sx].getName() == "v_x" && gdata.om[gpu::fluidConst_q_sy].getName() == "v_y" && gdata.om[gpu::fluidConst_q_sz].getName() == "v_z");

	if (qch == 0) {
		queue.submit(celerity::allow_by_ref, [=, &gdata](celerity::handler& cgh) {

			celerity::accessor om_acc{gdata.omSYCL[qch], cgh, celerity::access::one_to_one{}, celerity::read_write};
			celerity::accessor nomSYCL_acc{nomSYCL, cgh, celerity::access::all{}, celerity::read_only};
			celerity::accessor omSYCL_save_acc{omSYCL_save[0], cgh, celerity::access::one_to_one{}, celerity::read_only};

			cgh.parallel_for<class RKStepKernel>(omRange, [=](celerity::item<3> item) {
				size_t i = item.get_id(0);
				size_t j = item.get_id(1);
				size_t k = item.get_id(2);

				if (i >= -B && i < nom_max0 + B && j >= -B && j < nom_max1 + B && k >= -B && k < nom_max2 + B) {
					if (substep == 0) {
						om_acc[i][j][k] -= dt * nomSYCL_acc[i][j*nom_max1 + k][qch];
					} else {
						om_acc[i][j][k] = 0.5*omSYCL_save_acc[i][j][k] + 0.5*om_acc[i][j][k] 
									- 0.5 * dt * nomSYCL_acc[i][j*nom_max1 + k][qch];
					}
				}
			});
		});
	} else if (qch > 0 && qch < 4) { 
		queue.submit(celerity::allow_by_ref, [=, &gdata](celerity::handler& cgh) {

			celerity::accessor om_acc{gdata.omSYCL[qch], cgh, celerity::access::one_to_one{}, celerity::read_write};
			celerity::accessor om_rho_acc{gdata.omSYCL[0], cgh, celerity::access::one_to_one{}, celerity::read_only};
			celerity::accessor nomSYCL_acc{nomSYCL, cgh, celerity::access::all{}, celerity::read_only};
			celerity::accessor omSYCL_save_acc{omSYCL_save[0], cgh, celerity::access::one_to_one{}, celerity::read_only};

			cgh.parallel_for<class RKStepKernel>(omRange, [=](celerity::item<3> item) {
				size_t i = item.get_id(0);
				size_t j = item.get_id(1);
				size_t k = item.get_id(2);

				if (i >= -B && i < nom_max0 + B && j >= -B && j < nom_max1 + B && k >= -B && k < nom_max2 + B) {
					if (substep == 0) {
						if (isVelocity) {
							om_acc[i][j][k] = om_acc[i][j][k] * om_rho_acc[i][j][k] - (dt * nomSYCL_acc[i][j*nom_max1 + k][qch]) ;
						} else {
							om_acc[i][j][k] -= dt * nomSYCL_acc[i][j*nom_max1 + k][qch];
						}
					} else {
						if (isVelocity) {
							om_acc[i][j][k] = 0.5*omSYCL_save_acc[i][j][k] + 0.5*om_acc[i][j][k]*om_rho_acc[i][j][k] 
									- 0.5 * dt * nomSYCL_acc[i][j*nom_max1 + k][qch];
						} else {
							om_acc[i][j][k] = 0.5*omSYCL_save_acc[i][j][k] + 0.5*om_acc[i][j][k] 
									- 0.5 * dt * nomSYCL_acc[i][j*nom_max1 + k][qch];
						}
					}
				}
			});
		});
	} else if (qch == 4) {
		queue.submit(celerity::allow_by_ref, [=, &gdata](celerity::handler& cgh) {

			celerity::accessor om_acc{gdata.omSYCL[4], cgh, celerity::access::one_to_one{}, celerity::read_write};
			celerity::accessor om_rho_acc{gdata.omSYCL[0], cgh, celerity::access::one_to_one{}, celerity::read_only};
			celerity::accessor om_sx_acc{gdata.omSYCL[1], cgh, celerity::access::one_to_one{}, celerity::read_only};
			celerity::accessor om_sy_acc{gdata.omSYCL[2], cgh, celerity::access::one_to_one{}, celerity::read_only};
			celerity::accessor om_sz_acc{gdata.omSYCL[3], cgh, celerity::access::one_to_one{}, celerity::read_only};
			celerity::accessor nomSYCL_acc{nomSYCL, cgh, celerity::access::all{}, celerity::read_only};
			celerity::accessor omSYCL_save_acc{omSYCL_save[0], cgh, celerity::access::one_to_one{}, celerity::read_only};

			cgh.parallel_for<class RKStepKernel>(omRange, [=](celerity::item<3> item) {
				size_t i = item.get_id(0);
				size_t j = item.get_id(1);
				size_t k = item.get_id(2);

				if (i >= -B && i < nom_max0 + B && j >= -B && j < nom_max1 + B && k >= -B && k < nom_max2 + B) {
					if (substep == 0) {
							double Ekin = 0.5*(sqr(om_sx_acc[i][j][k]) + 
										sqr(om_sy_acc[i][j][k]) +
										sqr(om_sz_acc[i][j][k])) * om_rho_acc[i][j][k];
							om_acc[i][j][k] = (om_acc[i][j][k] + Ekin) - dt * nomSYCL_acc[i][j*nom_max1 + k][qch];
					} else {
						if(isVelocity && ENERGETICS == FULL) {
							double om_data = om_acc[i][j][k] + 0.5*(sqr(om_sx_acc[i][j][k]) + 
															sqr(om_sy_acc[i][j][k]) +
															sqr(om_sz_acc[i][j][k])) * om_rho_acc[i][j][k];
							om_acc[i][j][k] = 0.5*omSYCL_save_acc[i][j][k] + 0.5*om_data 
										- 0.5 * dt * nomSYCL_acc[i][j*nom_max1 + k][qch];
						} else {
							om_acc[i][j][k] = 0.5*omSYCL_save_acc[i][j][k] + 0.5*om_acc[i][j][k] 
									- 0.5 * dt * nomSYCL_acc[i][j*nom_max1 + k][qch];
						}
					}
				}
			});
		});
	}
    queue.slow_full_sync();
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
                      CelerityBuffer<double, 3> nomSYCL,
                      Pot om[], const int substep, size_t nom_max[3]) {}

void VanLeerIntegrator::Substep(Queue &queue, const Data &gdata, CelerityRange<3> omRange,
		CelerityBuffer<double, 3> nomSYCL, Pot om[], const int substep, size_t nom_max1, size_t nom_max[3]) {}


void TransformOmSycl_Cons2Prim(Queue &queue, const Data &gdata, CelerityRange<3> omRange, size_t nom_max[3]) {

	size_t nom_max0 = nom_max[0];
	size_t nom_max1 = nom_max[1];
	size_t nom_max2 = nom_max[2];
	int B = -3;

	queue.submit(celerity::allow_by_ref, [=, &gdata](celerity::handler& cgh) {

		celerity::accessor om_rho_acc{gdata.omSYCL[0], cgh, celerity::access::one_to_one{}, celerity::read_write};
		celerity::accessor om_sx_acc{gdata.omSYCL[1], cgh, celerity::access::one_to_one{}, celerity::read_write};
		celerity::accessor om_sy_acc{gdata.omSYCL[2], cgh, celerity::access::one_to_one{}, celerity::read_write};
		celerity::accessor om_sz_acc{gdata.omSYCL[3], cgh, celerity::access::one_to_one{}, celerity::read_write};
		celerity::accessor om_Eges_acc{gdata.omSYCL[4], cgh, celerity::access::one_to_one{}, celerity::read_write};

		cgh.parallel_for<class RKStepKernel>(omRange, [=](celerity::item<3> item) {
			size_t i = item.get_id(0);
			size_t j = item.get_id(1);
			size_t k = item.get_id(2);

			if (i >= -B && i < nom_max0 + B && j >= -B && j < nom_max1 + B && k >= -B && k < nom_max2 + B) {

				om_sx_acc[i][j][k] /= om_rho_acc[i][j][k];
				om_sy_acc[i][j][k] /= om_rho_acc[i][j][k];
				om_sz_acc[i][j][k] /= om_rho_acc[i][j][k];

				om_Eges_acc[i][j][k] -= 0.5*(sqr(om_sx_acc[i][j][k]) + 
									sqr(om_sy_acc[i][j][k]) +
									sqr(om_sz_acc[i][j][k])) * om_rho_acc[i][j][k];
			}
		});
	});
	queue.slow_full_sync();

}
