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
#include <unistd.h>

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
//sleep(1);cout << "TransCons2Prim" << endl;
	// User defined transform
	Problem.TransCons2Prim(gdata);
//sleep(1);cout << "TransCons2Prim" << endl;
	// Check for nans in primitive variables
	gdata.CheckNan(1);
//sleep(1);cout << "CheckNan" << endl;
	// Compute the carbuncle-flag if necessary
	if(gdata.use_carbuncleFlag) {
		Riemann[DirX]->compute_carbuncleFlag(gdata);
	}

	double cfl_lin(0.), cfl_eta(0.);

	// range of ghost cells
	int n_ghost[3];
	for(int dir=0; dir<3; ++dir) {
		//if(gdata.mx[dir]<2) {
		//	n_ghost[dir] = 0;
		//} else {
			n_ghost[dir] = 3;
		//}
	}

//	phys_fields_0D physVals[9];

	//std::vector<std::unique_ptr<phys_fields_0D>> physVals;
	//for (int inum = 0; inum < 9; ++inum) {
	//	physVals.push_back(std::make_unique<phys_fields_0D>(gdata, inum, gdata.fluid));
	//}

	//std::vector<std::unique_ptr<num_fields_0D>> numVals;
	//for(int idir=0; idir<3; ++idir) {
	//	numVals.push_back(std::make_unique<num_fields_0D>(gdata, gdata.fluid));
	//}

	//phys_fields_0D physVals_oldE(gdata, 1, gdata.fluid);
	//phys_fields_0D physVals_oldN(gdata, 3, gdata.fluid);
	//phys_fields_0D physVals_oldT(gdata, 5, gdata.fluid);

	//num_fields_0D numVals_oldE(gdata, gdata.fluid);
	//num_fields_0D numVals_oldN(gdata, gdata.fluid);
	//num_fields_0D numVals_oldT(gdata, gdata.fluid);

	// TODO PHILGS: we need to specify a direction here, not sure whether it is used...
	//Reconstruction_Block reconst(gdata, 0, gdata.fluid);

	//cronos::vector<double> pos(0,0,0);
	//cronos::vector<int> ipos(0,0,0);

	int n_omInt = gdata.fluid.get_N_OMINT();
	//NumArray<double> uPriOld_E(n_omInt);
	//NumArray<double> flux_numOld_E(n_omInt);
	
	// TODO PHILGS: this variable was uninitialized
	//double ptotalOld_E = 0.0;

	//NumMatrix<double,2> uPriOld_N(Index::set(0,-n_ghost[0]), Index::set(n_omInt,gdata.mx[0]+n_ghost[0]));
	//NumMatrix<double,2> uPriCurr_N(Index::set(0,-n_ghost[0]), Index::set(n_omInt,gdata.mx[0]+n_ghost[0]));
	//NumMatrix<double,2> flux_numOld_N(Index::set(0,-n_ghost[0]), Index::set(n_omInt,gdata.mx[0]+n_ghost[0]));
	//NumMatrix<double,2> flux_numCurr_N(Index::set(0,-n_ghost[0]), Index::set(n_omInt,gdata.mx[0]+n_ghost[0]));
	//NumMatrix<double,1> ptotalOld_N(Index::set(-n_ghost[0]), Index::set(gdata.mx[0]+n_ghost[0]));
	//NumMatrix<double,1> ptotalCurr_N(Index::set(-n_ghost[0]), Index::set(gdata.mx[0]+n_ghost[0]));

	//NumMatrix<double,3> uPriOld_T(Index::set(0,-n_ghost[0],-n_ghost[1]), Index::set(n_omInt,gdata.mx[0]+n_ghost[0],gdata.mx[1]+n_ghost[1]));
	//NumMatrix<double,3> uPriCurr_T(Index::set(0,-n_ghost[0],-n_ghost[1]), Index::set(n_omInt,gdata.mx[0]+n_ghost[0],gdata.mx[1]+n_ghost[1]));
	//NumMatrix<double,3> flux_numOld_T(Index::set(0,-n_ghost[0],-n_ghost[1]), Index::set(n_omInt,gdata.mx[0]+n_ghost[0],gdata.mx[1]+n_ghost[1]));
	//NumMatrix<double,3> flux_numCurr_T(Index::set(0,-n_ghost[0],-n_ghost[1]), Index::set(n_omInt,gdata.mx[0]+n_ghost[0],gdata.mx[1]+n_ghost[1]));
	//NumMatrix<double,2> ptotalOld_T(Index::set(-n_ghost[0],-n_ghost[1]), Index::set(gdata.mx[0]+n_ghost[0],gdata.mx[1]+n_ghost[1]));
	//NumMatrix<double,2> ptotalCurr_T(Index::set(-n_ghost[0],-n_ghost[1]), Index::set(gdata.mx[0]+n_ghost[0],gdata.mx[1]+n_ghost[1]));

	// TODO: make this a perfectly nested loop nest, increment dimensions of NumMatrix (4D using std::vector) to always save the full state
	// TODO: make a second kernel that is 3D and leave the old code as-is
	//std::vector<Buffer<double, 3>> uPriOld(n_omInt);
	//std::vector<Buffer<double, 3>> uPriCur(n_omInt);
	//std::vector<Buffer<double, 3>> flux_numOld(n_omInt);
	//std::vector<Buffer<double, 3>> flux_numCur(n_omInt);
	//std::vector<Buffer<double, 3>> ptotalOld(n_omInt);
	//std::vector<Buffer<double, 3>> ptotalCur(n_omInt);
	//for (int q = 0; q < n_omInt; ++q) {
	//	const int xMax = gdata.mx[0] + n_ghost[0] * 2;
	//	const int yMax = gdata.mx[1] + n_ghost[1] * 2;
	//	const int zMax = gdata.mx[2] + n_ghost[2] * 2;
	//	//uPriOld[q] = Buffer<double, 3>(Range<3>(zMax, yMax, xMax));
	//	uPriCur[q] = Buffer<double, 3>(Range<3>(zMax, yMax, xMax));
	//	//flux_numOld[q] = Buffer<double, 3>(Range<3>(zMax, yMax, xMax));
	//	flux_numCur[q] = Buffer<double, 3>(Range<3>(zMax, yMax, xMax));
	//	//ptotalOld[q] = Buffer<double, 3>(Range<3>(zMax, yMax, xMax));
	//	ptotalCur[q] = Buffer<double, 3>(Range<3>(zMax, yMax, xMax));
	//}

	//std::vector<Buffer<double, 2>> physValsSYCL(gdata.omSYCL.size(), Buffer<double, 2>(Range<2>(gpu::FaceMax, gpu::TypeMax)));

	const int izStart = -n_ghost[2] + 1;
	const int izEnd = gdata.mx[2] + n_ghost[2] - 1;
	const int iyStart = -n_ghost[1] + 1;
	const int iyEnd = gdata.mx[1] + n_ghost[1] - 1;
	const int ixStart = -n_ghost[0] + 1;
	const int ixEnd = gdata.mx[0] + n_ghost[0] - 1;

	/*for (int q = 0; q < gdata.omSYCL.size(); ++q) {
		cout << gdata.omSYCL[q].get_range().size(); << "\n"
	}*/

	//auto range = Range<3>(izEnd-izStart, iyEnd-iyStart, ixEnd-ixStart);

	//for (int q = 0; q < gdata.omSYCL.size(); ++q) {
	//	queue.submit([&](sycl::handler& cgh) {
	//		auto r = gdata.omSYCL[q].get_range();
	//		//std::cout << "dims: " << r.get(0) << ", " << r.get(1) << ", " << r.get(2) << ", " << std::endl << std::flush;
	//		auto om_acc = gdata.omSYCL[q].get_access<cl::sycl::access::mode::read>(cgh);
	//		auto physVals_acc = physValsSYCL[q].get_access<cl::sycl::access::mode::read_write>(cgh);
	//		auto uPriCur_acc = uPriCur[q].get_access<cl::sycl::access::mode::read_write>(cgh);
	//		auto flux_numCur_acc = flux_numCur[q].get_access<cl::sycl::access::mode::read_write>(cgh);
	//		auto ptotalCur_acc = ptotalCur[q].get_access<cl::sycl::access::mode::read_write>(cgh);
	//		cgh.parallel_for<class ComputeReconst>(range, [=](Item<3> item) {
	//			size_t iz = item.get(0) - izStart;
	//			size_t iy = item.get(1) - iyStart;
	//			size_t ix = item.get(2) - ixStart;
	//			gpu::compute(om_acc, physVals_acc, ix, iy, iz);

	//			// TODO: uPriOld_E is missing
	//			physVals_acc[]
	//			physVals_oldN.uPri(q_ind) = uPriOld_N(q_ind, ix);
	//			physVals_oldT.uPri(q_ind) = uPriOld_T(q_ind, ix, iy);
	//		});
	//	});
	//}

	auto computeStep = [](Reconstruction_Block reconst, const std::unique_ptr<Transformations>& Trafo, const std::unique_ptr<PhysFluxes>& PhysFlux, const std::vector<std::unique_ptr<RiemannSolver>>& Riemann, const ProblemType& Problem, const std::unique_ptr<EquationOfState>& eos, const Data& gdata, int ix, int iy, int iz, auto& max_cfl_lin) {

		// Reconstruction at given position
		std::vector<phys_fields_0D> physVals;
		for (int inum = 0; inum < 6; ++inum) {
			physVals.push_back(phys_fields_0D(gdata, inum, gdata.fluid));
		}
		reconst.compute(gdata, physVals, ix, iy, iz);

		if (ix >= -1 && iy >= -1 && iz >= -1) {

			std::vector<phys_fields_0D> physValsOld;
			for (int inum = 0; inum < 6; ++inum) {
				physValsOld.push_back(phys_fields_0D(gdata, inum, gdata.fluid));
			}

			reconst.compute(gdata, physValsOld, ix - 1, iy, iz, DirX);
			reconst.compute(gdata, physValsOld, ix, iy - 1, iz, DirY);
			reconst.compute(gdata, physValsOld, ix, iy, iz - 1, DirZ);

			physVals[gpu::FaceEast] = physValsOld[gpu::FaceEast];
			physVals[gpu::FaceNorth] = physValsOld[gpu::FaceNorth];
			physVals[gpu::FaceTop] = physValsOld[gpu::FaceTop];
		}			
		
		for (int dir = 0; dir < DirMax; ++dir) {

			int ixOff = (dir == DirX) ? ix : ix - 1;
			int iyOff = (dir == DirY) ? iy : iy - 1;
			int izOff = (dir == DirZ) ? iz : iz - 1;

			int face = dir * 2;

			Trafo->get_Cons(gdata, Problem, *eos, physVals[face], ix, iy, iz, face);
			PhysFlux->get_PhysFlux(gdata, Problem, physVals[face], ix, iy, iz, face);
			// i-1,j,k - East
			Trafo->get_Cons(gdata, Problem, *eos, physVals[face + 1], ix, iy, iz, face + 1);
			PhysFlux->get_PhysFlux(gdata, Problem, physVals[face + 1], ixOff, iyOff, izOff, face + 1);

		}

		std::vector<num_fields_0D> numVals;
		for (int idir = 0; idir < 3; ++idir) {
			numVals.push_back(num_fields_0D(gdata, gdata.fluid));
		}

		if (ix >= -1 && iy >= -1 && iz >= -1) {
			for (int dir = 0; dir < DirMax; ++dir) {
				int face = dir * 2;
				double cfl_loc = get_vChar3(gdata, Problem, physVals[face], physVals[face + 1], numVals[dir], dir);
				max_cfl_lin.combine(cfl_loc);
				get_NumFlux2(gdata, physVals[face + 1], physVals[face], numVals[dir], dir);
			}
		}

		return numVals;

	};

	auto range = Range<3>(izEnd-izStart, iyEnd-iyStart, ixEnd-ixStart);
	
	//necessary due to a cuda memory freeing error for some reason
	sleep(2);
	//std::vector<CelerityBuffer<double, 2>> physValsSYCL;
	//for (int i = 0; i < gdata.omSYCL.size(); i++) physValsSYCL.push_back(CelerityBuffer<double, 2>(celerity::range<2>(2,2)));
	
	//buffers
	//std::vector<CelerityBuffer<double, 2>> physValsSYCL(gdata.omSYCL.size(), CelerityBuffer<double, 2>( /*reeval*/ celerity::range<2>(6,2)));
	//std::vector<CelerityBuffer<double, 2>> physValsSYCL_Old(gdata.omSYCL.size(), CelerityBuffer<double, 2>( /*reeval*/ celerity::range<2>(6,2)));
	//std::vector<CelerityBuffer<double, 2>> physPtotalPthermSYCL(gdata.omSYCL.size(), CelerityBuffer<double, 2>(celerity::range<2>(6,2)));
	//std::vector<CelerityBuffer<double, 2>> physPtotalPthermSYCL_Old(gdata.omSYCL.size(), CelerityBuffer<double, 2>(celerity::range<2>(6,2)));
	celerity::buffer<double, 1> max_buf{{1}};

	//looking for better solution later
	int fluidConst[] {gdata.fluid.get_q_rho(), gdata.fluid.get_q_sx(), gdata.fluid.get_q_sx(), gdata.fluid.get_q_sz(),
						gdata.fluid.get_q_Eges(), gdata.fluid.get_q_Eadd(), gdata.fluid.get_q_Bx(), gdata.fluid.get_q_By(),
						gdata.fluid.get_q_Bz()};

	const int problem_cs2 = Problem.get_cs2();
	const int fluidType = gdata.fluid.get_fluid_type();

	/*for (int q = 0; q < gdata.omSYCL.size()	; ++q) {
		queue.submit([=](celerity::handler& cgh) {
			celerity::accessor physVals_acc{physValsSYCL[q], cgh, celerity::access::one_to_one{}, celerity::write_only};	
			//celerity::accessor physValsOld_acc{physValsSYCL_Old[q], cgh, celerity::access::one_to_one{}, celerity::write_only};
			//celerity::accessor physPtotalPtherm_acc{physPtotalPthermSYCL[q], cgh, celerity::access::one_to_one{}, celerity::write_only};
			//celerity::accessor physPtotalPthermOld_acc{physPtotalPthermSYCL[q], cgh, celerity::access::one_to_one{}, celerity::write_only};
			cgh.parallel_for<class BufferInitializationKernel>(physValsSYCL[q].get_range(), [=](celerity::item<2> item) {
				physVals_acc[item.get_id(0)][item.get_id(1)] = 0;
			});
		});
	}*/

	for (int q = 0; q < gdata.omSYCL.size()	; ++q) {

		queue.submit(celerity::allow_by_ref, [=, &gdata](celerity::handler& cgh) {

			auto r = gdata.omSYCL[q].get_range();
			auto rd = celerity::reduction(max_buf, cgh, cl::sycl::maximum<double>{},
                                  cl::sycl::property::reduction::initialize_to_identity{});
			celerity::accessor om_acc{gdata.omSYCL[q], cgh, celerity::access::all{}, celerity::read_only};
			//celerity::accessor physVals_acc{physValsSYCL[q], cgh, celerity::access::all{}, celerity::read_write};	
			//celerity::accessor physValsOld_acc{physValsSYCL_Old[q], cgh, celerity::access::all{}, celerity::read_write};
			//celerity::accessor physPtotalPtherm_acc{physPtotalPthermSYCL[q], cgh, celerity::access::all{}, celerity::read_write};
			//celerity::accessor physPtotalPthermOld_acc{physPtotalPthermSYCL[q], cgh, celerity::access::all{}, celerity::read_write};

			cgh.parallel_for<class ReductionKernel>(range, rd, [=, &gdata](celerity::item<3> item, auto& max_cfl_lin) {

				size_t iz = item.get_id(0) - izStart;
				size_t iy = item.get_id(1) - iyStart;
				size_t ix = item.get_id(2) - ixStart;
				
				//no vector in device code
				double uPri[gpu::FaceMax][N_OMINT] = {{0.0}};
				double uCon[gpu::FaceMax][N_OMINT] = {{0.0}};
				double physFlux[gpu::FaceMax][N_OMINT] = {{0.0}};
				double physValPtherm[gpu::FaceMax] = {0};
				double physValPtotal[gpu::FaceMax] = {0};

				double uPriOld[gpu::FaceMax][N_OMINT] = {{0}};

				//if (ix == ixEnd && iy == gdata.mx[1]+n_ghost[1]-1 && iz == gdata.mx[2]+n_ghost[2]-1) cout << "reach limit: " << ix << "." <<  iy << "." << iz << "\n";
				//const int fluidType = Riemann[DirX]->get_Fluid_Type();

				//pointwise reconstruction (incomplete)
				gpu::compute(om_acc, uPri, ix, iy, iz, q);

				gpu::compute(om_acc, uPriOld, ix -1, iy, iz, q);
				gpu::compute(om_acc, uPriOld, ix, iy -1, iz, q);
				gpu::compute(om_acc, uPriOld, ix, iy, iz -1, q);

				uPri[gpu::FaceEast][0] = uPri[gpu::FaceEast][0];
				uPri[gpu::FaceEast][1] = uPri[gpu::FaceEast][1];
				uPri[gpu::FaceNorth][0] = uPri[gpu::FaceNorth][0];
				uPri[gpu::FaceNorth][1] = uPri[gpu::FaceNorth][1];
				uPri[gpu::FaceTop][0] = uPri[gpu::FaceTop][0];
				uPri[gpu::FaceTop][1] = uPri[gpu::FaceTop][1];
				
				for (int dir = 0; dir < DirMax; ++dir) {
					int ixOff = (dir == DirX) ? ix : ix - 1;
					int iyOff = (dir == DirY) ? iy : iy - 1;
					int izOff = (dir == DirZ) ? iz : iz - 1;

					int face = dir * 2;

					gpu::get_Cons(om_acc, uPri[face], uCon[face], physValPtherm[face], physValPtotal[face], ix, iy, iz, face,
							 		thermal, Problem.gamma, problem_cs2, *eos, fluidType, fluidConst);

					gpu::get_PhysFlux(om_acc, physFlux[face], uPri[face], uCon[face], physValPtherm[face], face, fluidConst);

					gpu::get_Cons(om_acc, uPri[face+1], uCon[face+1], physValPtherm[face+1], physValPtotal[face+1], ix, iy, iz, face+1,
							 		thermal, Problem.gamma, problem_cs2, *eos, fluidType, fluidConst);

					gpu::get_PhysFlux(om_acc, physFlux[face+1], uPri[face+1], uCon[face+1], physValPtherm[face+1], face+1, fluidConst);
				}

				/*if(ix >= 0 && ix <= gdata.mx[0] && iy >= 0 && iy <= gdata.mx[1] && iz >= 0 && iz <= gdata.mx[2]) {
					//std::vector<phys_fields_0D> physVals;
					for (int inum = 0; inum < 6; ++inum) {
						//physVals_acc[inum] = phys_fields_0D(gdata, inum, gdata.fluid);
					}
					//reconst_inKernel.compute(gdata, physVals, ix, iy, iz);
					//const auto numVals = computeStep(reconst, Trafo, PhysFlux, Riemann, Problem, eos, gdata, ix, iy, iz, max_cfl_lin);
				}*/
			});

			//cgh.parallel_for<class EmptyTestKernel>(celerity::range<2>(2,2), celerity::id<2>(2,2), [=](celerity::item<2> item) {int i = i + 1;});
		});
	}
/*cout << cfl_lin << " ";
	// Loop over grid-block
	for (int iz = -n_ghost[2]+1; iz <= gdata.mx[2]+n_ghost[2]-1; ++iz){

		for (int iy = -n_ghost[1]+1; iy <= gdata.mx[1]+n_ghost[1]-1; ++iy){

			for (int ix = ixStart; ix <= ixEnd; ++ix){
				if (ix == ixEnd && iy == gdata.mx[1]+n_ghost[1]-1 && iz == gdata.mx[2]+n_ghost[2]-1) cout << "reach limit: " << ix << "." <<  iy << "." << iz << "\n";
				const int fluidType = Riemann[DirX]->get_Fluid_Type();

				if(ix >= 0 && ix <= gdata.mx[0] && iy >= 0 && iy <= gdata.mx[1] && iz >= 0 && iz <= gdata.mx[2]) {

					const auto numVals = computeStep(reconst, Trafo, PhysFlux, Riemann, Problem, eos, gdata, ix, iy, iz, cfl_lin);
cout << cfl_lin << " ";
					//gdata.nom update n_OMINT
					const auto numValsX = computeStep(reconst, Trafo, PhysFlux, Riemann, Problem, eos, gdata, ix + 1, iy, iz, cfl_lin);
					get_Changes(gdata, numVals[DirX], numValsX[DirX], gdata.nom, ix, iy, iz, DirX, fluidType);
cout << cfl_lin << " ";
					const auto numValsY = computeStep(reconst, Trafo, PhysFlux, Riemann, Problem, eos, gdata, ix, iy + 1, iz, cfl_lin);
					get_Changes(gdata, numVals[DirY], numValsY[DirY], gdata.nom, ix, iy, iz, DirY, fluidType);
cout << cfl_lin << " ";
					const auto numValsZ = computeStep(reconst, Trafo, PhysFlux, Riemann, Problem, eos, gdata, ix, iy, iz + 1, cfl_lin);
					get_Changes(gdata, numVals[DirZ], numValsZ[DirZ], gdata.nom, ix, iy, iz, DirZ, fluidType);
cout << cfl_lin << endl;				}

			}

		}
	}*/

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
		
		TimeIntegratorGeneric[q]->Substep(gdata, Problem, gdata.nom[q], gdata.om, n);

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
