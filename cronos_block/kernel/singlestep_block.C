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
#include "transformations.H"
#include "RiemannSolverHD.H"
#include "PhysFluxesHD.H"
#include "utils.H"

#if(STANDALONE_USAGE == TRUE)
#include "utils.H"
#else
#include "queue.H"
#endif

using namespace std;

double HyperbolicSolver::singlestep(Data &gdata, gridFunc &gfunc,
                                  ProblemType &Problem, int n, Queue& queue)
{
  //! Block structured version of singlestep
  /*! 
   * This version of singlestep solves all directions simultaneously
   * with future applications like AMR in mind
   */

  if(eos == nullptr) {
	eos = std::make_unique<EquationOfState>(Problem);
  }
  if(sources == NULL) {
	sources = std::make_unique<SourceTerms>(gdata, Problem, IntegrateA, thermal);
  }

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
// Start the clock:
// ----------------------------------------------------------------

	gettimeofday(&tick, 0);
	gettimeofday(&tin1, 0);
	cstart = clock();  

// ----------------------------------------------------------------
// Relevant declarations
// ----------------------------------------------------------------

//	nom = new NumMatrix<double,3> [n_omInt];
//	for (int q = 0; q < n_omInt; ++q) {
//#if (FLUID_TYPE == CRONOS_MHD)
//		if((gdata.om[q].getName() == "B_x" ||
//		    gdata.om[q].getName() == "B_y" ||
//		    gdata.om[q].getName() == "B_z")) {
//			nom[q].resize(Index::set(0,0,0),
//			              Index::set(0,0,0));
//		} else {
//#endif
//			nom[q].resize(Index::set(0,0,0),
//			              Index::set(gdata.mx[0],gdata.mx[1],gdata.mx[2]));
//#if (FLUID_TYPE == CRONOS_MHD)
//		}
//}
//#endif
	for (int q = 0; q < n_omInt; ++q) {
		gdata.nom[q].clear();
	}

	// Prepare user fields if necessary
	#if (OMS_USER == TRUE)

//	nom_user = new NumMatrix<double,3> [n_omIntUser];
//	for (int q = 0; q < n_omIntUser; ++q) {
//		nom_user[q].resize(Index::set(0,0,0),
//		                   Index::set(gdata.mx[0],gdata.mx[1],gdata.mx[2]));
//		nom_user[q].clear();
//	}
	for (int q = 0; q < n_omIntUser; ++q) {
		gdata.nom_user[q].clear();
	}

#endif

// ----------------------------------------------------------------
//      Setup Buffer for RK-Step
// ----------------------------------------------------------------

	auto omRange = gdata.omSYCL[0].get_range();
	size_t nom_max[3] = {omRange.get(0), omRange.get(1), omRange.get(2)};

	for (int q = 0; q < N_OMINT; q++) {
		double om_temp[206][7][7];

		for (int i = 0; i < nom_max[0]; i++) {
			for (int j = 0; j < nom_max[1]; j++) {
				for (int k = 0; k < nom_max[2]; k++) {
					om_temp[i][j][k] = gdata.om[q](i-3,j-3,k-3);
				}
			}
		}

		queue.submit(celerity::allow_by_ref, [=, &gdata](celerity::handler& cgh) {
			celerity::accessor omSYCL_acc{gdata.omSYCL[q], cgh, celerity::access::all{}, celerity::write_only_host_task};
			cgh.host_task(celerity::on_master_node, [=, &gdata]{
				for (int i = 0; i < nom_max[0]; i++) {
					for (int j = 0; j < nom_max[1]; j++) {
						for (int k = 0; k < nom_max[2]; k++) {
							omSYCL_acc[i][j][k] = om_temp[i][j][k];
						}
					}
				}
			});
		});
	}

// ---------------------------------------------------------------	      
//      Trafo of variables to conservative form
//----------------------------------------------------------------

#if(CRSWITCH_DUAL_ENERGY == CRONOS_ON)
	if(n_omInt > 7) {
		Trafo->computeEntropyFromE(gdata, gfunc, Problem);
	}
#endif

#if (USE_COROTATION == CRONOS_ON)
	// Transform to inertial frame velocity
	Trafo->TransCorotToInert(gdata, gfunc, Problem);
#endif

	Trafo->TransPrim2Cons(gdata, gfunc, Problem);

	// User defined transform
	Problem.TransPrim2Cons(gdata);

// ---------------------------------------------------------------	      
//      Saving old variables for Runge-Kutta (conservative form)
//----------------------------------------------------------------

	if (n == 0) {
		for (int q = 0; q < n_omInt; ++q){
#if(FLUID_TYPE == CRONOS_HYDRO)
			gdata.om[n_Omega+q] = gdata.om[q];         // save om_n
#elif (FLUID_TYPE == CRONOS_MHD)
			if((q < q_Bx || q > q_Bz) || !IntegrateA) {
				gdata.om[n_Omega+q] = gdata.om[q];         // save om_n
			} else {
				gdata.om[n_Omega+q] = gdata.om[n_omInt+N_ADD+q-q_Bx];
			}
#endif
		}

#if (OMS_USER == TRUE)
		// Saving user-variables (if necessary)
		for (int q = 0; q < n_omIntUser; ++q) {
			gdata.om_user[n_OmegaUser+q] = gdata.om_user[q];
		}
#endif

	}

// ---------------------------------------------------------------	      
//      Trafo primitive variables for reconstruction
//----------------------------------------------------------------

	Trafo->TransCons2Prim(gdata, gfunc, Problem);

	// User defined transform
	Problem.TransCons2Prim(gdata);

	// Check for nans in primitive variables
	gdata.CheckNan(1);

	// Compute the carbuncle-flag if necessary
	if(gdata.use_carbuncleFlag) {
		Riemann[DirX]->compute_carbuncleFlag(gdata);
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

	Reconstruction_Block reconst(gdata, 0, gdata.fluid);

	int n_omInt = gdata.fluid.get_N_OMINT();

	const int izStart = -n_ghost[2] + 1;
	const int izEnd = gdata.mx[2] + n_ghost[2] - 1;
	const int iyStart = -n_ghost[1] + 1;
	const int iyEnd = gdata.mx[1] + n_ghost[1] - 1;
	const int ixStart = -n_ghost[0] + 1;
	const int ixEnd = gdata.mx[0] + n_ghost[0] - 1;
	
	//buffers	
	celerity::buffer<double, 1> max_buf{celerity::range{1}};
	//CelerityBuffer<double, 3> nomSYCL {celerity::range<3>(omRange.get(0), omRange.get(1) * omRange.get(2), N_OMINT)};
	CelerityBuffer<nom_t, 3> nomSYCL {celerity::range<3>(omRange.get(0), omRange.get(1), omRange.get(2))};
	CelerityBuffer<double, 3> uPriSYCL{celerity::range<3>(omRange.get(0) * omRange.get(1) * omRange.get(2), gpu::FaceMax, N_OMINT)};

	// std::vector<CelerityBuffer<double, 3>> nomSYCL;
	// for (int i = 0; i < gpu::FaceMax; i++) {
	// 	nomSYCL.push_back(CelerityBuffer<double, 3>(Range<3>(gdata.mx[0]+6 +1, gdata.mx[1]+6+1, gdata.mx[2]+6+1)));
	// }

	queue.submit(celerity::allow_by_ref, [=, &gdata](celerity::handler& cgh) {
		celerity::accessor uPri_acc{uPriSYCL, cgh, celerity::access::all{}, celerity::write_only, celerity::no_init};
		cgh.parallel_for<class PointwiseReconstructionKernel>(uPriSYCL.get_range(), [=](celerity::item<3> item) {
			uPri_acc[item.get_id(0)][item.get_id(1)][item.get_id(2)] = 0;
		});
	});

	queue.submit([=](celerity::handler& cgh) {
		celerity::accessor nomSYCL_acc{nomSYCL, cgh, celerity::access::one_to_one{}, celerity::write_only, celerity::no_init};
		// celerity::accessor nom_rho_acc{nomSYCL[0], cgh, celerity::access::one_to_one{}, celerity::write_only, celerity::no_init};
		// celerity::accessor nom_sx_acc{nomSYCL[1], cgh, celerity::access::one_to_one{}, celerity::write_only, celerity::no_init};
		// celerity::accessor nom_sy_acc{nomSYCL[2], cgh, celerity::access::one_to_one{}, celerity::write_only, celerity::no_init};
		// celerity::accessor nom_sz_acc{nomSYCL[3], cgh, celerity::access::one_to_one{}, celerity::write_only, celerity::no_init};
		// celerity::accessor nom_Eges_acc{nomSYCL[4], cgh, celerity::access::one_to_one{}, celerity::write_only, celerity::no_init};
		cgh.parallel_for<class BufferInitializationKernel>(nomSYCL/*[0]*/.get_range(), [=](celerity::item<3> item) {
			//nomSYCL_acc[item.get_id(0)][item.get_id(1)][item.get_id(2)] = 0;
			// nom_rho_acc[item.get_id(0)][item.get_id(1)][item.get_id(2)] = 0;
			// nom_sx_acc[item.get_id(0)][item.get_id(1)][item.get_id(2)] = 0;
			// nom_sy_acc[item.get_id(0)][item.get_id(1)][item.get_id(2)] = 0;
			// nom_sz_acc[item.get_id(0)][item.get_id(1)][item.get_id(2)] = 0;
			// nom_Eges_acc[item.get_id(0)][item.get_id(1)][item.get_id(2)] = 0;
			nomSYCL_acc[item.get_id(0)][item.get_id(1)][item.get_id(2)].rho = 0;
			nomSYCL_acc[item.get_id(0)][item.get_id(1)][item.get_id(2)].sx = 0;
			nomSYCL_acc[item.get_id(0)][item.get_id(1)][item.get_id(2)].sy = 0;
			nomSYCL_acc[item.get_id(0)][item.get_id(1)][item.get_id(2)].sz = 0;
			nomSYCL_acc[item.get_id(0)][item.get_id(1)][item.get_id(2)].Eges = 0;
		});
	});

	int fluidConst[] {gdata.fluid.get_q_rho(), gdata.fluid.get_q_sx(), gdata.fluid.get_q_sy(), gdata.fluid.get_q_sz(),
						gdata.fluid.get_q_Eges(), gdata.fluid.get_q_Eadd(), gdata.fluid.get_q_Bx(), gdata.fluid.get_q_By(),
						gdata.fluid.get_q_Bz()};

	const double problem_cs2 = Problem.get_cs2();
	const double problem_gamma = Problem.gamma;
	const double denominator = pow(Problem.rho0,1.-Problem.gamma);
	const double half_beta = (Problem.mag ? value((char*)"Plasma_beta")/2. : 1);
	const int fluidType = gdata.fluid.get_fluid_type();
	const bool thermal = Trafo->get_thermal();
	
	double idx[DIM];
	for (int i = 0; i < DIM; i++) {
		idx[i] = gdata.idx[i];
	}

// ---------------------------------------------------------------	      
//      Pointwise Reconstruction
//----------------------------------------------------------------

	auto range = Range<3>(gdata.mx[0] + 2*n_ghost[0] + 1, gdata.mx[1] + 2*n_ghost[1] + 1, gdata.mx[2] + 2*n_ghost[2] + 1);
	int rangeEnd[3] = {gdata.mx[0] + n_ghost[0] + 1, gdata.mx[1] + n_ghost[1] + 1, gdata.mx[2] + n_ghost[2] + 1};
	
	// for (int q = 0; q < N_OMINT	; ++q) {
	// 	queue.submit(celerity::allow_by_ref, [=, &gdata](celerity::handler& cgh) {

	// 		celerity::accessor om_acc{gdata.omSYCL[q], cgh, celerity::access::neighborhood{2,2,2}, celerity::read_only};
	// 		celerity::accessor uPri_acc{uPriSYCL, cgh, celerity::access::all{}, celerity::write_only, celerity::no_init};

	// 		cgh.parallel_for<class PointwiseReconstructionKernel>(range, [=](celerity::item<3> item) {

	// 			size_t ix = item.get_id(0);
	// 			size_t iy = item.get_id(1);
	// 			size_t iz = item.get_id(2);
	// 			size_t ixyz = (ix * range.get(1) + iy) * range.get(2) + iz;

	// 			if (ix >= n_ghost[0] && ix <= rangeEnd[0] && iy >= n_ghost[1] && iy <= rangeEnd[1] && iz >= n_ghost[2] && iz <= rangeEnd[2]) {
	// 				if (ix >= n_ghost[0] && iy >= n_ghost[1] && iz >= n_ghost[2]) {
	// 					auto [uPriWest, uPriSouth, uPriBottom] = gpu::computeWSB(om_acc, ix, iy, iz);
	// 					auto [uPriEast, uPriNorth, uPriTop] = gpu::computeENT(om_acc, ix, iy, iz);
	// 					uPri_acc[ixyz][gpu::FaceWest][q] = uPriWest;
	// 					uPri_acc[ixyz][gpu::FaceEast][q] = uPriEast;
	// 					uPri_acc[ixyz][gpu::FaceNorth][q] = uPriNorth;
	// 					uPri_acc[ixyz][gpu::FaceSouth][q] = uPriSouth;
	// 					uPri_acc[ixyz][gpu::FaceBottom][q] = uPriBottom;
	// 					uPri_acc[ixyz][gpu::FaceTop][q] = uPriTop;}
	// 				// } else {
	// 				// 	auto [uPriWest, uPriEast, uPriSouth, uPriNorth, uPriBottom, uPriTop] = gpu::computeLowerCorner(om_acc, ix, iy, iz);
	// 				// 	uPri_acc[ixyz][gpu::FaceWest][q] = uPriWest;
	// 				// 	uPri_acc[ixyz][gpu::FaceEast][q] = uPriEast;
	// 				// 	uPri_acc[ixyz][gpu::FaceNorth][q] = uPriNorth;
	// 				// 	uPri_acc[ixyz][gpu::FaceSouth][q] = uPriSouth;
	// 				// 	uPri_acc[ixyz][gpu::FaceBottom][q] = uPriBottom;
	// 				// 	uPri_acc[ixyz][gpu::FaceTop][q] = uPriTop;
	// 				// }
	// 			}
	// 		});
	// 	});
	// }

	queue.submit(celerity::allow_by_ref, [=, &gdata](celerity::handler& cgh) {

		celerity::accessor uPri_acc{uPriSYCL, cgh, celerity::access::all{}, celerity::read_write};
		celerity::accessor om_rho_acc{gdata.omSYCL[0], cgh, celerity::access::neighborhood{2,2,2}, celerity::read_only};
		celerity::accessor om_sx_acc{gdata.omSYCL[1], cgh, celerity::access::neighborhood{2,2,2}, celerity::read_only};
		celerity::accessor om_sy_acc{gdata.omSYCL[2], cgh, celerity::access::neighborhood{2,2,2}, celerity::read_only};
		celerity::accessor om_sz_acc{gdata.omSYCL[3], cgh, celerity::access::neighborhood{2,2,2}, celerity::read_only};
		celerity::accessor om_Eges_acc{gdata.omSYCL[4], cgh, celerity::access::neighborhood{2,2,2}, celerity::read_only};

		cgh.parallel_for<class PointwiseReconstructionKernel>(range, [=](celerity::item<3> item) {

			size_t ix = item.get_id(0);
			size_t iy = item.get_id(1);
			size_t iz = item.get_id(2);
			size_t ixyz = (ix * range.get(1) + iy) * range.get(2) + iz;

			if (ix >= n_ghost[0] && ix <= rangeEnd[0] && iy >= n_ghost[1] && iy <= rangeEnd[1] && iz >= n_ghost[2] && iz <= rangeEnd[2]) {
				auto [uPriWest0, uPriSouth0, uPriBottom0] = gpu::computeWSB(om_rho_acc, ix, iy, iz);
				auto [uPriEast0, uPriNorth0, uPriTop0] = gpu::computeENT(om_rho_acc, ix, iy, iz);
				uPri_acc[ixyz][gpu::FaceWest][0] = uPriWest0;
				uPri_acc[ixyz][gpu::FaceEast][0] = uPriEast0;
				uPri_acc[ixyz][gpu::FaceSouth][0] = uPriSouth0;
				uPri_acc[ixyz][gpu::FaceNorth][0] = uPriNorth0;
				uPri_acc[ixyz][gpu::FaceBottom][0] = uPriBottom0;
				uPri_acc[ixyz][gpu::FaceTop][0] = uPriTop0;

				auto [uPriWest1, uPriSouth1, uPriBottom1] = gpu::computeWSB(om_sx_acc, ix, iy, iz);
				auto [uPriEast1, uPriNorth1, uPriTop1] = gpu::computeENT(om_sx_acc, ix, iy, iz);
				uPri_acc[ixyz][gpu::FaceWest][1] = uPriWest1;
				uPri_acc[ixyz][gpu::FaceEast][1]= uPriEast1;
				uPri_acc[ixyz][gpu::FaceNorth][1] = uPriNorth1;
				uPri_acc[ixyz][gpu::FaceSouth][1] = uPriSouth1;
				uPri_acc[ixyz][gpu::FaceBottom][1] = uPriBottom1;
				uPri_acc[ixyz][gpu::FaceTop][1] = uPriTop1;

				auto [uPriWest2, uPriSouth2, uPriBottom2] = gpu::computeWSB(om_sy_acc, ix, iy, iz);
				auto [uPriEast2, uPriNorth2, uPriTop2] = gpu::computeENT(om_sy_acc, ix, iy, iz);
				uPri_acc[ixyz][gpu::FaceWest][2] = uPriWest2;
				uPri_acc[ixyz][gpu::FaceEast][2] = uPriEast2;
				uPri_acc[ixyz][gpu::FaceNorth][2] = uPriNorth2;
				uPri_acc[ixyz][gpu::FaceSouth][2] = uPriSouth2;
				uPri_acc[ixyz][gpu::FaceBottom][2] = uPriBottom2;
				uPri_acc[ixyz][gpu::FaceTop][2] = uPriTop2;

				auto [uPriWest3, uPriSouth3, uPriBottom3] = gpu::computeWSB(om_sz_acc, ix, iy, iz);
				auto [uPriEast3, uPriNorth3, uPriTop3] = gpu::computeENT(om_sz_acc, ix, iy, iz);
				uPri_acc[ixyz][gpu::FaceWest][3] = uPriWest3;
				uPri_acc[ixyz][gpu::FaceEast][3] = uPriEast3;
				uPri_acc[ixyz][gpu::FaceNorth][3] = uPriNorth3;
				uPri_acc[ixyz][gpu::FaceSouth][3] = uPriSouth3;
				uPri_acc[ixyz][gpu::FaceBottom][3] = uPriBottom3;
				uPri_acc[ixyz][gpu::FaceTop][3] = uPriTop3;

				auto [uPriWest4, uPriSouth4, uPriBottom4] = gpu::computeWSB(om_Eges_acc, ix, iy, iz);
				auto [uPriEast4, uPriNorth4, uPriTop4] = gpu::computeENT(om_Eges_acc, ix, iy, iz);
				uPri_acc[ixyz][gpu::FaceWest][4] = uPriWest4;
				uPri_acc[ixyz][gpu::FaceEast][4] = uPriEast4;
				uPri_acc[ixyz][gpu::FaceNorth][4] = uPriNorth4;
				uPri_acc[ixyz][gpu::FaceSouth][4] = uPriSouth4;
				uPri_acc[ixyz][gpu::FaceBottom][4] = uPriBottom4;
				uPri_acc[ixyz][gpu::FaceTop][4] = uPriTop4;
			}
			
		});
	});

	// queue.slow_full_sync();
	

// ---------------------------------------------------------------	      
//      Reduction Kernel
//----------------------------------------------------------------

	queue.submit([=](celerity::handler& cgh) {

		auto rd = celerity::reduction(max_buf, cgh, cl::sycl::maximum<double>{},
								cl::sycl::property::reduction::initialize_to_identity{});
		celerity::accessor nom_acc{nomSYCL, cgh, celerity::access::all{}, celerity::read_write};
		// celerity::accessor nom_rho_acc{nomSYCL[0], cgh, celerity::access::one_to_one{}, celerity::read_write};
		// celerity::accessor nom_sx_acc{nomSYCL[1], cgh, celerity::access::one_to_one{}, celerity::read_write};
		// celerity::accessor nom_sy_acc{nomSYCL[2], cgh, celerity::access::one_to_one{}, celerity::read_write};
		// celerity::accessor nom_sz_acc{nomSYCL[3], cgh, celerity::access::one_to_one{}, celerity::read_write};
		// celerity::accessor nom_Eges_acc{nomSYCL[4], cgh, celerity::access::one_to_one{}, celerity::read_write};
		celerity::accessor uPri_acc{uPriSYCL, cgh, celerity::access::all{}, celerity::read_only};

		cgh.parallel_for<class ReductionKernel>(range, rd, [=](celerity::item<3> item, auto& max_cfl_lin) {

			size_t ix = item.get_id(0);
			size_t iy = item.get_id(1);
			size_t iz = item.get_id(2);
			double cfl_lin = -100.0;

			if (ix >= n_ghost[0] && ix < rangeEnd[0] && iy >= n_ghost[1] && iy < rangeEnd[1] && iz >= n_ghost[2] && iz < rangeEnd[2]) {

				double numFlux[DirMax][N_OMINT] = {};
				double num_ptotal[DirMax] = {};

				size_t ixyz = (ix * range.get(1) + iy) * range.get(2) + iz;
				gpu::computeStep(uPri_acc, ix, iy, iz, ixyz, &cfl_lin, numFlux, num_ptotal, thermal, problem_gamma, problem_cs2,
													denominator, half_beta, fluidType, idx, fluidConst);

				double numFlux_Dir[DirMax][N_OMINT] = {};
				double num_ptotal_Dir[DirMax] = {};

				gpu::computeStep(uPri_acc, ix + 1, iy, iz, ixyz + range.get(1) * range.get(2), &cfl_lin, numFlux_Dir, num_ptotal_Dir,
													thermal, problem_gamma, problem_cs2, denominator, half_beta, fluidType, idx, fluidConst);
				gpu::get_Changes(nom_acc, ix, iy, iz, DirX, numFlux[DirX], num_ptotal[DirX], numFlux_Dir[DirX], num_ptotal_Dir[DirX],
										N_OMINT, nom_max[2], idx);
				// gpu::get_Changes(nom_rho_acc, 0, ix, iy, iz, DirX, numFlux[DirX], num_ptotal[DirX], numFlux_Dir[DirX], num_ptotal_Dir[DirX],
				// 						N_OMINT, nom_max[2], idx);
				// gpu::get_Changes(nom_sx_acc, 1, ix, iy, iz, DirX, numFlux[DirX], num_ptotal[DirX], numFlux_Dir[DirX], num_ptotal_Dir[DirX],
				// 						N_OMINT, nom_max[2], idx);
				// gpu::get_Changes(nom_sy_acc, 2, ix, iy, iz, DirX, numFlux[DirX], num_ptotal[DirX], numFlux_Dir[DirX], num_ptotal_Dir[DirX],
				// 						N_OMINT, nom_max[2], idx);
				// gpu::get_Changes(nom_sz_acc, 3, ix, iy, iz, DirX, numFlux[DirX], num_ptotal[DirX], numFlux_Dir[DirX], num_ptotal_Dir[DirX],
				// 						N_OMINT, nom_max[2], idx);
				// gpu::get_Changes(nom_Eges_acc, 4, ix, iy, iz, DirX, numFlux[DirX], num_ptotal[DirX], numFlux_Dir[DirX], num_ptotal_Dir[DirX],
				// 						N_OMINT, nom_max[2], idx);


				gpu::computeStep(uPri_acc, ix, iy + 1, iz, ixyz + range.get(2), &cfl_lin, numFlux_Dir, num_ptotal_Dir, thermal, problem_gamma, problem_cs2,
													denominator, half_beta, fluidType, idx, fluidConst);
				gpu::get_Changes(nom_acc, ix, iy, iz, DirY, numFlux[DirY], num_ptotal[DirY], numFlux_Dir[DirY], num_ptotal_Dir[DirY],
										N_OMINT, nom_max[2], idx);
				// gpu::get_Changes(nom_rho_acc, 0, ix, iy, iz, DirY, numFlux[DirY], num_ptotal[DirY], numFlux_Dir[DirY], num_ptotal_Dir[DirY],
				// 						N_OMINT, nom_max[2], idx);
				// gpu::get_Changes(nom_sx_acc, 1, ix, iy, iz, DirY, numFlux[DirY], num_ptotal[DirY], numFlux_Dir[DirY], num_ptotal_Dir[DirY],
				// 						N_OMINT, nom_max[2], idx);
				// gpu::get_Changes(nom_sy_acc, 2, ix, iy, iz, DirY, numFlux[DirY], num_ptotal[DirY], numFlux_Dir[DirY], num_ptotal_Dir[DirY],
				// 						N_OMINT, nom_max[2], idx);
				// gpu::get_Changes(nom_sz_acc, 3, ix, iy, iz, DirY, numFlux[DirY], num_ptotal[DirY], numFlux_Dir[DirY], num_ptotal_Dir[DirY],
				// 						N_OMINT, nom_max[2], idx);
				// gpu::get_Changes(nom_Eges_acc, 4, ix, iy, iz, DirY, numFlux[DirY], num_ptotal[DirY], numFlux_Dir[DirY], num_ptotal_Dir[DirY],
				// 						N_OMINT, nom_max[2], idx);

				gpu::computeStep(uPri_acc, ix, iy, iz + 1, ixyz + 1, &cfl_lin, numFlux_Dir, num_ptotal_Dir, thermal, problem_gamma, problem_cs2,
													denominator, half_beta, fluidType, idx, fluidConst);
				gpu::get_Changes(nom_acc, ix, iy, iz, DirZ, numFlux[DirZ], num_ptotal[DirZ], numFlux_Dir[DirZ], num_ptotal_Dir[DirZ],
										N_OMINT, nom_max[2], idx);
				// gpu::get_Changes(nom_rho_acc, 0, ix, iy, iz, DirZ, numFlux[DirZ], num_ptotal[DirZ], numFlux_Dir[DirZ], num_ptotal_Dir[DirZ],
				// 						N_OMINT, nom_max[2], idx);
				// gpu::get_Changes(nom_sx_acc, 1, ix, iy, iz, DirZ, numFlux[DirZ], num_ptotal[DirZ], numFlux_Dir[DirZ], num_ptotal_Dir[DirZ],
				// 						N_OMINT, nom_max[2], idx);
				// gpu::get_Changes(nom_sy_acc, 2, ix, iy, iz, DirZ, numFlux[DirZ], num_ptotal[DirZ], numFlux_Dir[DirZ], num_ptotal_Dir[DirZ],
				// 						N_OMINT, nom_max[2], idx);
				// gpu::get_Changes(nom_sz_acc, 3, ix, iy, iz, DirZ, numFlux[DirZ], num_ptotal[DirZ], numFlux_Dir[DirZ], num_ptotal_Dir[DirZ],
				// 						N_OMINT, nom_max[2], idx);
				// gpu::get_Changes(nom_Eges_acc, 4, ix, iy, iz, DirZ, numFlux[DirZ], num_ptotal[DirZ], numFlux_Dir[DirZ], num_ptotal_Dir[DirZ],
				// 						N_OMINT, nom_max[2], idx);
			}

			max_cfl_lin.combine(cfl_lin);
			
		});
	});

	queue.slow_full_sync();
	
	queue.submit(celerity::allow_by_ref, [=, &cfl_lin](celerity::handler& cgh) {
		celerity::accessor max_buf_acc{max_buf, cgh, celerity::access::all{}, celerity::read_only_host_task};
		cgh.host_task(celerity::on_master_node, [=, &cfl_lin]{
			cfl_lin = max_buf_acc[0];
		});
	});
	
	queue.slow_full_sync();

// ---------------------------------------------------------------	      
//      Former Kernel for Non-Parallel Execution
//----------------------------------------------------------------

	// auto computeStep = [](Reconstruction_Block reconst, const std::unique_ptr<Transformations>& Trafo, const std::unique_ptr<PhysFluxes>& PhysFlux, const std::vector<std::unique_ptr<RiemannSolver>>& Riemann, const ProblemType& Problem, const std::unique_ptr<EquationOfState>& eos, const Data& gdata, int ix, int iy, int iz, double& cfl_lin) {

	// 	// Reconstruction at given position
	// 	std::vector<phys_fields_0D> physVals;
	// 	for (int inum = 0; inum < 6; ++inum) {
	// 		physVals.push_back(phys_fields_0D(gdata, inum, gdata.fluid));
	// 	}
	// 	reconst.compute(gdata, physVals, ix, iy, iz);

	// 	if (ix >= -1 && iy >= -1 && iz >= -1) {

	// 		std::vector<phys_fields_0D> physValsOld;
	// 		for (int inum = 0; inum < 6; ++inum) {
	// 			physValsOld.push_back(phys_fields_0D(gdata, inum, gdata.fluid));
	// 		}

	// 		reconst.compute(gdata, physValsOld, ix - 1, iy, iz, DirX);
	// 		reconst.compute(gdata, physValsOld, ix, iy - 1, iz, DirY);
	// 		reconst.compute(gdata, physValsOld, ix, iy, iz - 1, DirZ);

	// 		physVals[gpu::FaceEast] = physValsOld[gpu::FaceEast];
	// 		physVals[gpu::FaceNorth] = physValsOld[gpu::FaceNorth];
	// 		physVals[gpu::FaceTop] = physValsOld[gpu::FaceTop];
	// 	}			
		
	// 	for (int dir = 0; dir < DirMax; ++dir) {

	// 		int ixOff = (dir == DirX) ? ix : ix - 1;
	// 		int iyOff = (dir == DirY) ? iy : iy - 1;
	// 		int izOff = (dir == DirZ) ? iz : iz - 1;

	// 		int face = dir * 2;

	// 		Trafo->get_Cons(gdata, Problem, *eos, physVals[face], ix, iy, iz, face);
	// 		PhysFlux->get_PhysFlux(gdata, Problem, physVals[face], ix, iy, iz, face);
	// 		// i-1,j,k - East
	// 		Trafo->get_Cons(gdata, Problem, *eos, physVals[face + 1], ix, iy, iz, face + 1);
	// 		PhysFlux->get_PhysFlux(gdata, Problem, physVals[face + 1], ixOff, iyOff, izOff, face + 1);

	// 	}

	// 	std::vector<num_fields_0D> numVals;
	// 	for (int idir = 0; idir < 3; ++idir) {
	// 		numVals.push_back(num_fields_0D(gdata, gdata.fluid));
	// 	}

	// 	if (ix >= -1 && iy >= -1 && iz >= -1) {
	// 		for (int dir = 0; dir < DirMax; ++dir) {
	// 			int face = dir * 2;
	// 			get_vChar2(gdata, Problem, physVals[face], physVals[face + 1], numVals[dir], dir, cfl_lin);
	// 			get_NumFlux2(gdata, physVals[face + 1], physVals[face], numVals[dir], dir);
	// 		}
	// 	}

	// 	return numVals;

	// };

	// for (int iz = -n_ghost[2]+1; iz <= gdata.mx[2]+n_ghost[2]-1; ++iz){

	// 	for (int iy = -n_ghost[1]+1; iy <= gdata.mx[1]+n_ghost[1]-1; ++iy){

	// 		for (int ix = ixStart; ix <= ixEnd; ++ix){
	// 			//if (ix == ixEnd && iy == gdata.mx[1]+n_ghost[1]-1 && iz == gdata.mx[2]+n_ghost[2]-1) cout << "reach limit: " << ix << "." <<  iy << "." << iz << "\n";
	// 			const int fluidType = Riemann[DirX]->get_Fluid_Type();

	// 			if(ix >= 0 && ix <= gdata.mx[0] && iy >= 0 && iy <= gdata.mx[1] && iz >= 0 && iz <= gdata.mx[2]) {
	// 				//cout << "seq kernel iteration " << ix << " " << iy << " " << iz <<"\n";
	// 				const auto numVals = computeStep(reconst, Trafo, PhysFlux, Riemann, Problem, eos, gdata, ix, iy, iz, cfl_lin);

	// 				//gdata.nom update n_OMINT
	// 				const auto numValsX = computeStep(reconst, Trafo, PhysFlux, Riemann, Problem, eos, gdata, ix + 1, iy, iz, cfl_lin);
	// 				get_Changes(gdata, numVals[DirX], numValsX[DirX], gdata.nom, ix, iy, iz, DirX, fluidType);

	// 				const auto numValsY = computeStep(reconst, Trafo, PhysFlux, Riemann, Problem, eos, gdata, ix, iy + 1, iz, cfl_lin);
	// 				get_Changes(gdata, numVals[DirY], numValsY[DirY], gdata.nom, ix, iy, iz, DirY, fluidType);

	// 				const auto numValsZ = computeStep(reconst, Trafo, PhysFlux, Riemann, Problem, eos, gdata, ix, iy, iz + 1, cfl_lin);
	// 				get_Changes(gdata, numVals[DirZ], numValsZ[DirZ], gdata.nom, ix, iy, iz, DirZ, fluidType);
		
	// 			}
	// 		}

	// 	}
	// }

// ----------------------------------------------------------------
//   Check for errors:
// ----------------------------------------------------------------

	

	for(int q = 0; q<n_omInt; ++q) {
		CheckNan(gdata.nom[q],q, 0, 1,"nom");
	}

// ----------------------------------------------------------------
//   Compute Courant number
// ----------------------------------------------------------------

	double cfl = compute_cfl(gdata, Problem, cfl_eta, cfl_lin, n);
	cout << "cfl: " << cfl << endl;
	//double cfl = 0.0;
// ----------------------------------------------------------------
//   Geometrical source terms:
// ----------------------------------------------------------------

#if (GEOM != CARTESIAN)
	sources->src_Geom(gdata, Problem, gdata.nom);
#endif

// ----------------------------------------------------------------
//   Geometrical source terms:
// ----------------------------------------------------------------

#if (USE_COROTATION == CRONOS_ON)
	Trafo->src_Corotating(gdata, Problem, gdata.nom);
#endif

// ----------------------------------------------------------------
//   User defined source terms:
// ----------------------------------------------------------------

#if (USE_COROTATION == CRONOS_ON)
	// Transform to co-rotating frame velocity for user src
	Trafo->TransInertToCorot(gdata, gfunc, Problem);
#endif

	Problem.src_User(gdata, gdata.nom, gdata.nom_user);

#if (USE_COROTATION == CRONOS_ON)
	// Transform back to inertial frame
	Trafo->TransCorotToInert(gdata, gfunc, Problem);
#endif

// ----------------------------------------------------------------
//   Check for errors again before applying the changes:
// ----------------------------------------------------------------

	for(int q = 0; q<n_omInt; ++q) {
		CheckNan(gdata.nom[q],q, 0, 2,"nom");
	}

// ----------------------------------------------------------------
//   Transform to conservative variables:
// ----------------------------------------------------------------

	Trafo->TransPrim2Cons(gdata, gfunc, Problem);

	Problem.TransPrim2Cons(gdata);

// ----------------------------------------------------------------
//   Determine domain to be integrated and apply changes:
// ----------------------------------------------------------------

	for (int q = 0; q < n_omInt; ++q) {
		
		TimeIntegratorGeneric[q]->Substep(queue, gdata, omRange, nomSYCL, gdata.om, n, nom_max);
		//TimeIntegratorGeneric[q]->Substep(gdata, Problem, gdata.nom[q], gdata.om, n);
	}

#if (OMS_USER == TRUE)

	for (int q=0; q<n_omIntUser; ++q) {

		TimeIntegratorUser[q]->Substep(gdata, Problem, gdata.nom_user[q],
		                              gdata.om_user, n);
		
	}

#endif

//	delete [] nom;
//#if (OMS_USER == TRUE)
//	delete [] nom_user;
//#endif

	gettimeofday(&tock, 0);
	cstep = clock();

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

	gettimeofday(&tock2, 0);
#if (USE_ANGULAR_MOMENTUM == TRUE)
	Trafo->TransAngMom2Vel(gdata, gfunc, Problem);
#else
	Trafo->TransMomen2Vel(gdata, gfunc, Problem);
#endif

#if (FLUID_TYPE == CRONOS_MHD)
	if(IntegrateA) {
		compute_B(gdata, gfunc, Problem);
	}
#endif

#if (USE_COROTATION == CRONOS_ON)
	// Transform to co-rotating frame velocity for BCs
	Trafo->TransInertToCorot(gdata, gfunc, Problem);
#endif

#if(ENERGETICS == FULL)
	int q_max = q_Eges;
#else
	int q_max = n_omInt;
#endif

	// Boundary conditions
	for(int q=0; q<q_max; ++q) {
		gfunc.boundary(gdata, Problem, gdata.om[q],B,q);
	}

	if(ENERGETICS == FULL) {

#if (USE_COROTATION == CRONOS_ON)
		// Transform back for energy trafo
		Trafo->TransCorotToInert(gdata, gfunc, Problem);
#endif

		if(thermal) {
#if(CRSWITCH_DUAL_ENERGY == CRONOS_ON)
			Trafo->TransE2Eth(gdata, gfunc, Problem, n, true);
#else
			Trafo->TransE2Eth(gdata, gfunc, Problem);
#endif
		} else {
			Trafo->TransE2T(gdata, gfunc, Problem);
		}

		for(int q=q_Eges; q<n_omInt; ++q) {
			gfunc.boundary(gdata, Problem, gdata.om[q],B,q);
		}

#if (USE_COROTATION == CRONOS_ON)
		// Transform to co-rotating frame velocity for BCs
		Trafo->TransInertToCorot(gdata, gfunc, Problem);
#endif

	}

	Problem.TransCons2Prim(gdata);

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

	double delt = ((tock.tv_sec + tock.tv_usec/1.e6) - 
	             (tick.tv_sec + tick.tv_usec/1.e6));

	double delt2 = ((tock2.tv_sec + tock2.tv_usec/1.e6) - 
	              (tock.tv_sec + tock.tv_usec/1.e6));

	double dtStep = ((tstep.tv_sec + tstep.tv_usec/1.e6) -
	               (tstepOld.tv_sec + tstepOld.tv_usec/1.e6));

	if(gdata.rank == 0) {
		if(n == RK_STEPS-1 && Problem.get_Info() && Problem.checkout(5)) {
			cout << "------------------------------------------------------" << endl;
			cout << " Time needed for substeps:  " << delt << " " << delt2 << endl;
			cout << "             for full step: " << dtStep << endl;
//       cout << " CPU Cycle times: " << gdata.rank << " ";
//       cout << (1.*(cstep-cstart))/CLOCKS_PER_SEC << " ";
//       cout << (1.*(cend-cstep))/CLOCKS_PER_SEC << endl;
		}

	}

// ----------------------------------------------------------------
//   Test physical state
// ----------------------------------------------------------------

#if (USE_COROTATION == CRONOS_ON)
	// Transform to inertial frame velocity for phystest
	Trafo->TransCorotToInert(gdata, gfunc, Problem);
#endif

	phystest(gdata, gfunc, Problem, n);

#if (USE_COROTATION == CRONOS_ON)
	// Transform back to co-rotating frame
	Trafo->TransInertToCorot(gdata, gfunc, Problem);
#endif

	// Computing div B only at end of timestep
	if(Problem.mag && n == RK_STEPS-1){
		compute_divB(gdata, gfunc, Problem);
	}

//	cout << " my cfl: " << cfl << endl;

	return cfl;

}
