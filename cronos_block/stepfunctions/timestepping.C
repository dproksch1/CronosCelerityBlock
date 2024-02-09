#include <iostream>
#include <chrono>
#include <thread>

#include "timestepping.H"

static bool nom_temp_decl = false;

TimeIntegrator::TimeIntegrator(const int n_saves):
	n_saves(n_saves)
{}


TimeIntegrator::~TimeIntegrator()
{
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

}

void TimeIntegrator::init_omBuffer(Queue &queue, const int mx[])
{
	for (int i = 0; i < 5; i++) {
		omSYCL_save.push_back(CelerityBuffer<double, 3>(CelerityRange<3>(mx[0]+6 +1, mx[1]+6+1, mx[2]+6+1)));
	}
}

RKSteps::RKSteps():
		TimeIntegrator(1)
{
	this->oneThird  = 1./3.;
	this->twoThirds = 2./3.;
}


void RKSteps::Substep(Queue &queue, const Data &gdata, CelerityRange<3> omRange,
		CelerityBuffer<nom_t, 3> nomSYCL,
		double dt, const int substep, size_t nom_max[3]) {

	int B = -3;

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
			
	}

	gdata.om[1].rename("v_x");
	gdata.om[2].rename("v_y");
	gdata.om[3].rename("v_z");
#endif
}


VanLeerIntegrator::VanLeerIntegrator():
		TimeIntegrator(1)
{
}

void VanLeerIntegrator::Substep(Queue &queue, const Data &gdata, CelerityRange<3> omRange,
                      CelerityBuffer<nom_t, 3> nomSYCL,
                      double dt, const int substep, size_t nom_max[3]) {}