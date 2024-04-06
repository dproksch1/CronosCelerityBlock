#include "gridgen.H"
#include <iomanip>
#include <cmath>
#include <sys/time.h>
#include <vector>
#include <string.h>
//#include "DissipationMHD.H"
#include "reconst_block.H"
#include "changes.H"
#include "standalone_usage.H"
#include "compute_step.H"
#include <unistd.h>

#include "reconst_block.H"
#include "transformations_block.H"
#include "utils.H"
#include "queue.H"
#include "lazy_assertion.H"

using namespace std;
  
/*! 
 * @brief structured Celerity version of singlestep
 * @details This version of singlestep solves all directions simultaneously and uses Celerity kernels
 */
double HyperbolicSolver::singlestep(Data &gdata, GridFunc &gfunc,
                                  ProblemType &Problem, int n, Queue& queue)
{

// ----------------------------------------------------------------
// Checking for negative Temperatures
// ----------------------------------------------------------------

	if(gdata.rank == 0) {
		if(ENERGETICS == FULL && TempNeg !=0) {
			cout << " Negative Temperatures: " << TempNeg << endl;
			TempNeg = 0;
		}
	}

// ----------------------------------------------------------------
//      Setup Buffer for RK-Step
// ----------------------------------------------------------------

	gettimeofday(&tick, 0);
	
	auto omRange = gdata.omSYCL[0].get_range();
	size_t nom_max[3] = {omRange.get(0), omRange.get(1), omRange.get(2)};

	if (gdata.tstep <= 0 && n == 0) {
		for (int q = 0; q < N_OMINT; q++) {

			queue.submit(celerity::allow_by_ref, [=, &gdata](celerity::handler& cgh) {
				celerity::accessor omSYCL_acc{gdata.omSYCL[q], cgh, celerity::access::one_to_one{}, celerity::write_only_host_task};
				cgh.host_task(omRange, [=, &gdata](celerity::partition<3> part){

					const auto& sr = part.get_subrange();

					for (int i = sr.offset[0]; i < sr.offset[0] + sr.range[0]; i++) {
						for (int j = sr.offset[1]; j < sr.offset[1] + sr.range[1]; j++) {
							for (int k = sr.offset[2]; k < sr.offset[2] + sr.range[2]; k++) {
								omSYCL_acc[i][j][k] = gdata.om[q](i-3,j-3,k-3);
							}
						}
					}
				});
			});
		}

		queue.submit(celerity::allow_by_ref, [=, &gdata](celerity::handler& cgh) {
			celerity::accessor nomSYCL_acc{gdata.nomSYCL, cgh, celerity::access::one_to_one{}, celerity::write_only, celerity::no_init};
			cgh.parallel_for<class BufferInitializationKernel>(gdata.nomSYCL.get_range(), [=](celerity::item<3> item) {
				for (int d = 0; d < N_OMINT; d++) {
					nomSYCL_acc[item.get_id(0)][item.get_id(1)][item.get_id(2)].mat[d] = 0;
				}
			});
		});

		TimeIntegratorGeneric[0]->init_omBuffer(queue, gdata.mx);
	}


// ----------------------------------------------------------------
// Start the clock:
// ----------------------------------------------------------------

	gettimeofday(&tin11, 0);
	gettimeofday(&tin1, 0);
	cstart = clock();

// ---------------------------------------------------------------	      
//  Setup Carbuncle Flag
//----------------------------------------------------------------

	// Compute the carbuncle-flag if necessary
	if(gdata.use_carbuncleFlag) {
		compute_carbuncleFlag(queue, gdata, Riemann[DirX]->getAlphaCarbuncle());
	}

// ---------------------------------------------------------------	      
//      Buffer and Constant Initialization
//----------------------------------------------------------------

	double cfl_lin(0.), cfl_eta(0.);

	// range of ghost cells
	int n_ghost[3];
	for(int dir=0; dir<3; ++dir) {
			n_ghost[dir] = 3;
	}

	int n_omInt = gdata.fluid.get_N_OMINT();

	int fluidConst[] {gdata.fluid.get_q_rho(), gdata.fluid.get_q_sx(), gdata.fluid.get_q_sy(), gdata.fluid.get_q_sz(),
						gdata.fluid.get_q_Eges(), gdata.fluid.get_q_Eadd(), gdata.fluid.get_q_Bx(), gdata.fluid.get_q_By(),
						gdata.fluid.get_q_Bz()};

	const double problem_cs2 = Problem.get_cs2();
	const double problem_gamma = Problem.gamma;
	const double denominator = pow(Problem.rho0,1.-Problem.gamma);
	const double half_beta = (Problem.mag ? value((char*)"Plasma_beta")/2. : 1);
	const int fluidType = gdata.fluid.get_fluid_type();
	const bool thermal = Trafo->get_thermal();
	const bool use_carbuncle = gdata.use_carbuncleFlag;
	const double Theta = value_exists((char*)"Theta") ? value((char*)"Theta") : 1;

	double idx[DIM];
	for (int i = 0; i < DIM; i++) {
		idx[i] = gdata.idx[i];
	}

	auto range = CelerityRange<3>(gdata.mx[0] + 2*n_ghost[0] + 1, gdata.mx[1] + 2*n_ghost[1] + 1, gdata.mx[2] + 2*n_ghost[2] + 1);
	int rangeEnd[3] = {gdata.mx[0] + n_ghost[0] + 1, gdata.mx[1] + n_ghost[1] + 1, gdata.mx[2] + n_ghost[2] + 1};


// ---------------------------------------------------------------	      
//      Compute Kernel
//----------------------------------------------------------------

	queue.submit(celerity::allow_by_ref, [=, &gdata](celerity::handler& cgh) {

		auto rd = celerity::reduction(gdata.cflSYCL, cgh, cl::sycl::maximum<double>{},
								cl::sycl::property::reduction::initialize_to_identity{});
		celerity::accessor nom_acc{gdata.nomSYCL, cgh, celerity::access::one_to_one{}, celerity::read_write};
		celerity::accessor om_rho_acc{gdata.omSYCL[0], cgh, celerity::access::neighborhood{2,2,2}, celerity::read_only};
		celerity::accessor om_sx_acc{gdata.omSYCL[1], cgh, celerity::access::neighborhood{2,2,2}, celerity::read_only};
		celerity::accessor om_sy_acc{gdata.omSYCL[2], cgh, celerity::access::neighborhood{2,2,2}, celerity::read_only};
		celerity::accessor om_sz_acc{gdata.omSYCL[3], cgh, celerity::access::neighborhood{2,2,2}, celerity::read_only};
		celerity::accessor om_Eges_acc{gdata.omSYCL[4], cgh, celerity::access::neighborhood{2,2,2}, celerity::read_only};
		celerity::accessor carbuncle_flag{gdata.carbuncleFlagSYCL[0], cgh, celerity::access::neighborhood{2,2,2}, celerity::read_only};

		cgh.parallel_for<class ComputeKernel>(range, rd, [=](celerity::item<3> item, auto& max_cfl_lin) {

			size_t ix = item.get_id(0);
			size_t iy = item.get_id(1);
			size_t iz = item.get_id(2);
			double cfl_lin = -100.0;

			if (ix >= n_ghost[0] && ix < rangeEnd[0] && iy >= n_ghost[1] && iy < rangeEnd[1] && iz >= n_ghost[2] && iz < rangeEnd[2]) {

				double numFlux[DirMax][N_OMINT] = {};
				double num_ptotal[DirMax] = {};

				gpu::computeStep(om_rho_acc, om_sx_acc, om_sy_acc, om_sz_acc, om_Eges_acc, ix, iy, iz, &cfl_lin, numFlux, num_ptotal, carbuncle_flag, thermal, problem_gamma, problem_cs2,
													denominator, half_beta, fluidType, Theta, use_carbuncle, idx, fluidConst);

				double numFlux_Dir[DirMax][N_OMINT] = {};
				double num_ptotal_Dir[DirMax] = {};

				gpu::computeStep(om_rho_acc, om_sx_acc, om_sy_acc, om_sz_acc, om_Eges_acc, ix + 1, iy, iz, &cfl_lin, numFlux_Dir, num_ptotal_Dir,
													carbuncle_flag, thermal, problem_gamma, problem_cs2, denominator, half_beta, fluidType, Theta, use_carbuncle, idx, fluidConst);
				gpu::get_Changes(nom_acc, ix, iy, iz, DirX, numFlux[DirX], num_ptotal[DirX], numFlux_Dir[DirX], num_ptotal_Dir[DirX],
										N_OMINT, nom_max[2], idx);

				gpu::computeStep(om_rho_acc, om_sx_acc, om_sy_acc, om_sz_acc, om_Eges_acc,  ix, iy + 1, iz, &cfl_lin, numFlux_Dir, num_ptotal_Dir, carbuncle_flag, thermal, problem_gamma, problem_cs2,
													denominator, half_beta, fluidType, Theta, use_carbuncle, idx, fluidConst);
				gpu::get_Changes(nom_acc, ix, iy, iz, DirY, numFlux[DirY], num_ptotal[DirY], numFlux_Dir[DirY], num_ptotal_Dir[DirY],
										N_OMINT, nom_max[2], idx);

				gpu::computeStep(om_rho_acc, om_sx_acc, om_sy_acc, om_sz_acc, om_Eges_acc, ix, iy, iz + 1, &cfl_lin, numFlux_Dir, num_ptotal_Dir, carbuncle_flag, thermal, problem_gamma, problem_cs2,
													denominator, half_beta, fluidType, Theta, use_carbuncle, idx, fluidConst);
				gpu::get_Changes(nom_acc, ix, iy, iz, DirZ, numFlux[DirZ], num_ptotal[DirZ], numFlux_Dir[DirZ], num_ptotal_Dir[DirZ],
										N_OMINT, nom_max[2], idx);
			}

			max_cfl_lin.combine(cfl_lin);
			
		});
	});

// ----------------------------------------------------------------
//   Compute Courant number
// ----------------------------------------------------------------

	gdata.fetch_cfl(queue);
	cfl_lin = gdata.cfl;

	gettimeofday(&tock, 0);
	cstep = clock();

	double cfl = compute_cfl(gdata, Problem, cfl_eta, cfl_lin, n);
	cout << "cfl: " << cfl << endl;

// ----------------------------------------------------------------
//   User defined source terms:
// ----------------------------------------------------------------

	Problem.src_User(gdata, gdata.nom, gdata.nom_user);

// ----------------------------------------------------------------
//   Transform to conservative variables:
// ----------------------------------------------------------------

	Trafo->TransPrim2Cons(queue, gdata, gfunc, Problem);
	Problem.TransPrim2Cons(queue, gdata);

// ----------------------------------------------------------------
//   Determine domain to be integrated and apply changes:
// ----------------------------------------------------------------

	TimeIntegratorGeneric[0]->Substep(queue, gdata, omRange, gdata.nomSYCL, gdata.dt, n);

#if (OMS_USER == TRUE)

	for (int q=0; q<n_omIntUser; ++q) {

		TimeIntegratorUser[q]->Substep(gdata, Problem, gdata.nom_user[q],
		                              gdata.om_user, n);
		
	}

#endif

// ----------------------------------------------------------------
//   Determine time at intermediate steps:
// ----------------------------------------------------------------

#if (RK_STEPS == 3) 
	// Time must come before boundary for time dependent boundaries
	if (n == 0) gdata.time += gdata.dt;      // First trial step with t* = t+dt
	if (n == 1) gdata.time -= 0.5*gdata.dt;  // Second trial step with t* = t+0.5dt
	if (n == 2) gdata.time += 0.5*gdata.dt;      // Last Runge-Kutta step - next time step
#else
	if (n == 0) gdata.time += gdata.dt;
#endif

// ----------------------------------------------------------------
//   Transformation to primitive variables and apply bcs
// ----------------------------------------------------------------

#if(ENERGETICS == FULL)
	int q_max = q_Eges;
#else
	int q_max = n_omInt;
#endif

	// Boundary conditions
	for(int q=0; q<q_max; ++q) {
		gfunc.boundary(queue,gdata, Problem, B, q);
	}

	if(ENERGETICS == FULL) {

		Trafo->TransE2Eth(queue, gdata, gfunc, Problem, 0, true);
 
		for(int q=q_Eges; q<n_omInt; ++q) {
			gfunc.boundary(queue,gdata, Problem, B, q);
		}

	}

	Problem.TransCons2Prim(queue, gdata);

#if (OMS_USER == TRUE)

	for(int q=0; q<n_omIntUser; ++q) {
		gfunc.boundary(gdata, Problem, gdata.om_user[q],B,q);
	}
#endif

// ----------------------------------------------------------------
//   Runtime estimates
// ----------------------------------------------------------------

	gettimeofday(&tock2, 0);
	cend = clock();

	timeval tstepOld = tstep;
	if(n == RK_STEPS-1) {
		gettimeofday(&tstep, 0);
	}

	double delt = ((tin11.tv_sec + tin11.tv_usec/1.e6) - 
	             (tick.tv_sec + tick.tv_usec/1.e6));

	double delt2 = ((tock.tv_sec + tock.tv_usec/1.e6) - 
				(tin11.tv_sec + tin11.tv_usec/1.e6));

	double delt3 = ((tock2.tv_sec + tock2.tv_usec/1.e6) - 
	              (tock.tv_sec + tock.tv_usec/1.e6));

	double delt4 = delt2 + delt3;

	if(gdata.rank == 0) {
		if(Problem.get_Info() && Problem.checkout(5)) {
			cout << "------------------------------------------------------" << endl;
			cout << " Time needed for substeps:  " << delt << " " << delt2 << " " << delt3 << endl;
			cout << "             for full step: " << delt4 << endl;
		}

	}

// ----------------------------------------------------------------
//   compute divB
// ----------------------------------------------------------------

	// Computing div B only at end of timestep
	if(Problem.mag && n == RK_STEPS-1){
		compute_divB(gdata, gfunc, Problem);
	}

	return cfl;
}
