#include "gridgen.H"
#include <iomanip>
#include <cmath>
#include <vector>
#include "DissipationMHD.H"
#include "timewrapper.H"
#include "queue.H"

//void __assert_fail(const char * assertion, const char * file, unsigned int line, const char * function) {
//	abort();
//}

using namespace std;

#if(FLUID_TYPE==CRONOS_MULTIFLUID)
#include "singlestep_multifluid.C"
#else

REAL HyperbolicSolver::singlestep(Data &gdata, gridFunc &gfunc,
                                  ProblemType &Problem, int n, Queue& queue)
{
	if(eos == nullptr/*(EquationOfState*)0*/) {
		// TODO PHILGS: only initializes a few scalars
		eos = std::make_unique<EquationOfState>(Problem);
	}
	if(sources == NULL) {
		// TODO PHILGS: only initializes a few scalars
		sources = std::make_unique<SourceTerms>(gdata, Problem, IntegrateA, thermal);
	}

	// TODO PHILGS: can be ignored, does not apply in the HYDRO case
	Save = new Saves(gdata);

// ----------------------------------------------------------------
// Checking if there are negative Temperatures
// ----------------------------------------------------------------

	if(gdata.rank == 0) {
		if(ENERGETICS == FULL && TempNeg !=0) {
			cout << " Negative Temperatures: " << TempNeg << endl;
			TempNeg = 0;
		}
	}

// ----------------------------------------------------------------
// Starting the clock:
// ----------------------------------------------------------------

	gettimeofday(&tick, 0);
	gettimeofday(&tin1, 0);
	cstart = clock();  

// ----------------------------------------------------------------
// Relevant declarations
// ----------------------------------------------------------------

//	nom = new NumMatrix<REAL,3> [n_omInt];
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
//#endif
//		nom[q].clear();
//	}

	for (int q = 0; q < n_omInt; ++q) {
		gdata.nom[q].clear();
	}


#if (OMS_USER == TRUE)

	for (int q = 0; q < n_omIntUser; ++q) {
		gdata.nom_user[q].clear();
	}

#endif

// ---------------------------------------------------------------	      
//      Aux transormations
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


//---------------------------------------------------------------
//		Fill mass flux array with zeros
//---------------------------------------------------------------

#if(GEOM == SPHERICAL)
	gdata.massFlux.clear();
#endif

	// TODO PHILGS: Accesses: gdata.om, AllReduce and throws exception if NaN found anywhere
	gdata.CheckNan(1);

	// Compute the carbuncle-flag if necessary
	// TODO PHILGS: Does not seem to be used for the simple test case SodToro5HllcRef, ignore for now
	//if(gdata.use_carbuncleFlag) {
		//RiemannX->compute_carbuncleFlag(gdata);
	//}

	REAL cfl_lin(0.), cfl_eta(0.);

	cronos::vector<REAL> pos(0,0,0);
	cronos::vector<REAL> ipos(0, 0, 0);

	// x-direction
	int dir=0;

	//Buffer<double, 2> local_omLoc(cl::sycl::range<2>(n_omInt, gdataSize));
	
	for(int k=-1; k<=gdata.mx[2]+1; ++k) {
		// Update of perp part of position:
		//pos.set(2, gdata.getCen_z(k));
		ipos.set(2, k);
		for (int j=-1; j<=gdata.mx[1]+1; ++j){
			// Update of perp part of position:
			//pos.set(1, gdata.getCen_y(j));
			ipos.set(1, j);

			// ORIG
			// Assigning 1D array data
			for (int q = 0; q < n_omInt; ++q) {
				if ((FLUID_TYPE == CRONOS_MHD && (q<q_By || q>q_Bz)) ||
					FLUID_TYPE != CRONOS_MHD) {
					for (int i = -3; i <= gdata.mx[dir]+2; ++i) {
						//fields[DIR_X]->omLocORIG[q](i - 3) = gdata.om[q](i - 3, j, k);
						fieldsX->omLocORIG[q](i) = gdata.om[q](i, j, k);
					}
#ifdef USE_SYCL
					fieldsX->omLocSYCL[q] = extractPoleToSYCL(gdata.om[q], gdata.mx[dir] + 6, dir, -3, j, k);
#endif
				}
			}

			//for (int q = 0; q < n_omInt; ++q) {
			//	//assert_sycl_eq(fields[DIR_X]->omLocORIG[q], fields[DIR_X]->omLocSYCL[q]) << "HyperbolicSolver::singlestep:" << __LINE__ << " checking equality for q=" << q << std::endl;
			//	assert_sycl_eq(fieldsX->omLocORIG[q], fieldsX->omLocSYCL[q]) << "HyperbolicSolver::singlestep:" << __LINE__ << " checking equality for q=" << q << std::endl;
			//}

			// Assigning shifted fields for magnetic induction
#if (FLUID_TYPE == CRONOS_MHD)
			for (int i = -3; i <= gdata.mx[0]+2; ++i){
				fieldsX->omLocP[q_By](i) = gdata.om[q_By](i,j,k);
			}
			for (int i = -3; i <= gdata.mx[0]+2; ++i){
				fieldsX->omLocM[q_By](i) = gdata.om[q_By](i,j-1,k);
			}
			for (int i = -3; i <= gdata.mx[0]+2; ++i){
				fieldsX->omLocP[q_Bz](i) = gdata.om[q_Bz](i,j,k);
			}
			for (int i = -3; i <= gdata.mx[0]+2; ++i){
				fieldsX->omLocM[q_Bz](i) = gdata.om[q_Bz](i,j,k-1);
			}
#endif
			// Assing 1D array data for user fields:
#if (OMS_USER == TRUE)
			for (int q = 0; q < n_omIntUser; ++q){
				for (int i = -3; i <= gdata.mx[0]+2; ++i){
					fieldsXUser->omLoc[q](i) = gdata.om_user[q](i,j,k);
				}
			}
#endif

			// Assigning carbuncle flag values if necessary
			//if(gdata.use_carbuncleFlag) {
			//	for (int i = -2; i <= gdata.mx[0]+2; ++i){
			//		fieldsX->carbuncle_flag(i) = gdata.carbuncleFlag(i,j,k);
			//	}
			//}

			// Doing the reconstruction:
			ReconstX[n]->compute(queue, gdata, *fieldsX, *physValxL, *physValxR);
			//ReconstX[n]->compute(queue, gdata, *fields[DIR_X], *physValxL, *physValxR);

#if (OMS_USER == TRUE)
			// Reconstruction for User-data:
			ReconstXUser[n]->compute(gdata, *fieldsXUser,
			                      *physValxLUser, *physValxRUser);
#endif

			// Get conservative variables:
			// on LHS
			Trafo->get_Cons(queue, gdata, Problem, *eos, ipos, *physValxL, dir, -0.5);
			// on RHS
			Trafo->get_Cons(queue, gdata, Problem, *eos, ipos, *physValxR, dir,  0.5);

#if (OMS_USER == TRUE) 
			// Get conservative variables (needed for Riemann solver)
			// -- default is just a copy
			// LHS
			Problem.get_Cons(gdata, ipos, *physValxL, *physValxLUser, dir, -0.5);
			// Trafo->get_ConsUser(gdata, *physValxLUser, n_omIntUser);
			// RHS
			Problem.get_Cons(gdata, ipos, *physValxR, *physValxRUser, dir,  0.5);
			// Trafo->get_ConsUser(gdata, *physValxRUser, n_omIntUser);
#endif

			// Compute physical fluxes:
			// on LHS
			PhysFlux->get_PhysFlux(queue, gdata, Problem, ipos, *physValxL, dir, -0.5);
			// on RHS
			PhysFlux->get_PhysFlux(queue, gdata, Problem, ipos, *physValxR, dir,  0.5);


			// Possible modification by user
			// TODO PHILGS: Not used by current integration tests, disable for now
			//Problem.get_PhysFlux(gdata, ipos, *physValxL, dir, -0.5);

			// TODO PHILGS: Not used by current integration tests, disable for now
			//Problem.get_PhysFlux(gdata, ipos, *physValxR, dir, 0.5);


#if (OMS_USER == TRUE)
			// Compute physical fluxes for user fields:
			// on LHS
// 			PhysFluxUser->get_PhysFlux(gdata, Problem, ipos, 
// 			                           *physValxL, *physValxLUser,
// 			                           dir, -0.5);
			Problem.get_PhysFluxUser(gdata, ipos, *physValxL, *physValxLUser,
			                         dir, -0.5);
			// on RHS
// 			PhysFluxUser->get_PhysFlux(gdata, Problem, ipos, 
// 			                           *physValxR, *physValxRUser,
// 			                           dir,  0.5);
			Problem.get_PhysFluxUser(gdata, ipos, *physValxR, *physValxRUser,
			                         dir,  0.5);

#endif

			// Compute characteristic velocities:
			RiemannX->get_vChar(queue, gdata, Problem, ipos, *physValxL, *physValxR,
					fieldsX->v_ch_mORIG, fieldsX->v_ch_pORIG, fieldsX->v_ch_mSYCL, fieldsX->v_ch_pSYCL, dir, cfl_lin);
			//RiemannX->get_vChar(gdata, Problem, ipos, *physValxL, *physValxR,
			//	fields[DIR_X]->v_ch_m, fields[DIR_X]->v_ch_p, dir, cfl_lin);
//			RiemannX->get_vChar(gdata, Problem, ipos, *physValxL
//			                    *physValxR, *fieldsX, dir, cfl_lin);

#if (OMS_USER == TRUE)
			fieldsXUser->v_ch_p = fieldsX->v_ch_p;
			fieldsXUser->v_ch_m = fieldsX->v_ch_m;
#endif


			// Compute numerical fluxes:
			RiemannX->get_NumFlux(queue, gdata, *physValxL, *physValxR, *fieldsX, dir);
			//RiemannX->get_NumFlux(gdata, *physValxL, *physValxR, *fields[DIR_X], dir);

#if (OMS_USER == TRUE)
			RiemannXUser->get_NumFlux(gdata, *physValxLUser, *physValxRUser,
			                          *fieldsXUser, dir);

#endif


			if((j>=0 && j<=gdata.mx[1]) && (k>=0 && k<=gdata.mx[2])) {

				get_Changes(gdata, fieldsX->fluxORIG, gdata.nom, fieldsX->ptotalORIG,
				            0, j, k, n_omInt, RiemannX->get_Fluid_Type());
				//get_Changes(gdata, fields[DIR_X]->flux, gdata.nom, fields[DIR_X]->ptotal,
				//	0, j, k, n_omInt, RiemannX->get_Fluid_Type());

#if (OMS_USER == TRUE)

				get_Changes(gdata, fieldsXUser->flux, gdata.nom_user,
				            fieldsXUser->ptotal,
				            0, j, k, n_omIntUser,
				            RiemannXUser->get_Fluid_Type());
#endif

			}

			Save->save_vars(gdata, *fieldsX, 0, j, k);
			//Save->save_vars(gdata, *fields[DIR_X], 0, j, k);
      
		}
	}

	//for (int i = 0; i < 5; ++i) {
	//	dump(gdata.nom[i]);
	//}
	//exit(0);


	// Compute sum of mass flux along phi and theta for MPI parallel case
#if (GEOM == SPHERICAL)
#ifdef parallel
	NumMatrix<double,1> massFlux_global;
	int len = gdata.massFlux.getHigh(0) - gdata.massFlux.getLow(0) + 1;
	massFlux_global.resize(gdata.massFlux);
	MPI_Reduce(gdata.massFlux, massFlux_global, len, MPI_DOUBLE,
			MPI_SUM, 0, gdata.comm_constx);
	MPI_Bcast(massFlux_global, len, MPI_DOUBLE, 0, gdata.comm_constx);
	MPI_Barrier(gdata.comm_constx);
	gdata.massFlux = massFlux_global;
#endif
#endif


	// y-direction
	dir = 1;

	for (int k = -1; k <= gdata.mx[2]+1; ++k){
		// Update of perp part of position:
		pos.set(2, gdata.getCen_z(k));
		ipos.set(2, k);
		for (int i = -1; i <= gdata.mx[0]+1; ++i){
			// Update of perp part of position:
			pos.set(0, gdata.getCen_x(i));
			ipos.set(0, i);
      
			// Assigning 1D array data
			for (int q = 0; q < n_omInt; ++q){
				if((FLUID_TYPE == CRONOS_MHD && q!=q_Bx && q!=q_Bz) ||
				   FLUID_TYPE != CRONOS_MHD) {
					for (int j = -3; j <= gdata.mx[1]+2; ++j){
						fieldsY->omLocORIG[q](j) = gdata.om[q](i,j,k);
					}
				}
				fieldsY->omLocSYCL[q] = extractPoleToSYCL(gdata.om[q], gdata.mx[dir] + 6, dir, i, -3, k);
			}


			// Assigning shifted fields for magnetic induction
#if (FLUID_TYPE == CRONOS_MHD)			
			for (int j = -3; j <= gdata.mx[1]+2; ++j){
				fieldsY->omLocP[q_Bx](j) = gdata.om[q_Bx](i,j,k);
//				fieldsY->omLocP[q_Bx](j) = gdata.om[q_Bx](i,j,k) - Problem.get_FieldAnalytical(4,gdata.getEdgL_x(i+1),gdata.getCen_y(j),gdata.getCen_z(k));
			}
			for (int j = -3; j <= gdata.mx[1]+2; ++j){
				fieldsY->omLocM[q_Bx](j) = gdata.om[q_Bx](i-1,j,k);
//				fieldsY->omLocM[q_Bx](j) = gdata.om[q_Bx](i-1,j,k) - Problem.get_FieldAnalytical(4,gdata.getEdgL_x(i),gdata.getCen_y(j),gdata.getCen_z(k));
			}
			for (int j = -3; j <= gdata.mx[1]+2; ++j){
				fieldsY->omLocP[q_Bz](j) = gdata.om[q_Bz](i,j,k);
			}
			for (int j = -3; j <= gdata.mx[1]+2; ++j){
				fieldsY->omLocM[q_Bz](j) = gdata.om[q_Bz](i,j,k-1);
			}
#endif

#if (OMS_USER == TRUE)
			for (int q = 0; q < n_omIntUser; ++q){
				for (int j = -3; j <= gdata.mx[dir]+2; ++j){
					fieldsYUser->omLoc[q](j) = gdata.om_user[q](i,j,k);
				}
			}
#endif


			// Assigning carbuncle flag values if necessary
			if(gdata.use_carbuncleFlag) {
				for (int j = -2; j <= gdata.mx[1]+2; ++j){
					fieldsY->carbuncle_flag(j) = gdata.carbuncleFlag(i,j,k);
				}
			}

			// Doing the reconstruction:
			ReconstY[n]->compute(queue, gdata, *fieldsY, *physValyL, *physValyR);

#if (OMS_USER == TRUE)
			// Reconstruction for User-data:
			ReconstYUser[n]->compute(gdata, *fieldsYUser,
			                      *physValyLUser, *physValyRUser);
#endif



//			// Reset prim values:
//			for (int j = -2; j <= gdata.mx[1]+1; ++j){
//				if(false) {
//				physValyL->uPri[4](j) += Problem.get_FieldAnalytical(4,gdata.getCen_x(i),gdata.getEdgL_y(j),gdata.getCen_z(k));
//				physValyR->uPri[4](j) += Problem.get_FieldAnalytical(4,gdata.getCen_x(i),gdata.getEdgL_y(j+1),gdata.getCen_z(k));
//				}
//			}

      
			// Get conservative variables:
			// on LHS
			Trafo->get_Cons(queue, gdata, Problem, *eos, ipos, *physValyL, dir, -0.5);
			// on RHS
			Trafo->get_Cons(queue, gdata, Problem, *eos, ipos, *physValyR, dir,  0.5);

#if (OMS_USER == TRUE)
			// Get conservative variables (needed for Riemann solver)
			// -- just by copying for the time being
			// LHS
			Problem.get_Cons(gdata, ipos, *physValyL, *physValyLUser, dir, -0.5);
			// Trafo->get_ConsUser(gdata, *physValyLUser, n_omIntUser);
			// RHS
			Problem.get_Cons(gdata, ipos, *physValyR, *physValyRUser, dir, 0.5);
			// Trafo->get_ConsUser(gdata, *physValyRUser, n_omIntUser);
#endif


			// Compute physical fluxes:
			// on LHS
			PhysFlux->get_PhysFlux(queue, gdata, Problem, ipos, *physValyL, dir, -0.5);
			// on RHS
			PhysFlux->get_PhysFlux(queue, gdata, Problem, ipos, *physValyR, dir,  0.5);


			// Possible modification by user
			// TODO PHILGS: Not used by current integration tests, disable for now
			//Problem.get_PhysFlux(gdata, ipos, *physValyL, dir, -0.5);

			// TODO PHILGS: Not used by current integration tests, disable for now
			//Problem.get_PhysFlux(gdata, ipos, *physValyR, dir,  0.5);

#if (OMS_USER == TRUE)
			// Compute physical fluxes for user fields:
			// on LHS
// 			PhysFluxUser->get_PhysFlux(gdata, Problem, ipos, 
// 			                           *physValyL, *physValyLUser,
// 			                           dir, -0.5);
			Problem.get_PhysFluxUser(gdata, ipos, *physValyL, *physValyLUser,
			                         dir, -0.5);

			// on RHS
// 			PhysFluxUser->get_PhysFlux(gdata, Problem, ipos, 
// 			                           *physValyR, *physValyRUser,
// 			                           dir,  0.5);
			Problem.get_PhysFluxUser(gdata, ipos, *physValyR, *physValyRUser,
			                         dir,  0.5);
#endif
			
			// Compute characteristic velocities:
			RiemannY->get_vChar(queue, gdata, Problem, ipos, *physValyL, *physValyR,
					fieldsY->v_ch_mORIG, fieldsY->v_ch_pORIG, fieldsY->v_ch_mSYCL, fieldsY->v_ch_pSYCL, dir, cfl_lin);
//			RiemannY->get_vChar(gdata, Problem, ipos, *physValyL, *physValyR, *fieldsY,
//			                   dir, cfl_lin);

#if (OMS_USER == TRUE)
			fieldsYUser->v_ch_p = fieldsY->v_ch_p;
			fieldsYUser->v_ch_m = fieldsY->v_ch_m;
#endif

			// Compute numerical flux:
			RiemannY->get_NumFlux(queue, gdata, *physValyL, *physValyR, *fieldsY, dir);


#if (OMS_USER == TRUE)
			RiemannYUser->get_NumFlux(gdata, *physValyLUser, *physValyRUser,
						                          *fieldsYUser, dir);
#endif

			// Compute changes (noms)
			if((i>=0 && i<=gdata.mx[0]) && (k>=0 && k<=gdata.mx[2])) {
				get_Changes(gdata, fieldsY->fluxORIG, gdata.nom, fieldsY->ptotalORIG,
				            dir, i, k, n_omInt, RiemannY->get_Fluid_Type());
			}

#if (OMS_USER == TRUE)
			if((i>=0 && i<=gdata.mx[0]) && (k>=0 && k<=gdata.mx[2])) {

				get_Changes(gdata, fieldsYUser->flux, gdata.nom_user,
				            fieldsYUser->ptotal,
				            dir, i, k, n_omIntUser,
				            RiemannYUser->get_Fluid_Type());
			}
#endif


			// Save 1D data if necessary
			Save->save_vars(gdata, *fieldsY, dir, i, k);




			

		}
	}



	// z -direction
	dir = 2;

	for (int j = -1; j <= gdata.mx[1]+1; ++j){
		// Update of perp part of position:
		pos.set(1, gdata.getCen_y(j));
		ipos.set(1, j);
		for (int i = -1; i <= gdata.mx[0]+1; ++i){
			// Update of perp part of position:
			pos.set(0, gdata.getCen_x(i));
			ipos.set(0, i);

			// Assigning 1D array data
			for (int q = 0; q < n_omInt; ++q){
				if((FLUID_TYPE == CRONOS_MHD && (q<q_Bx || q>q_By)) ||
					FLUID_TYPE != CRONOS_MHD) {
					for (int k = -3; k <= gdata.mx[2]+2; ++k){
						fieldsZ->omLocORIG[q](k) = gdata.om[q](i,j,k);
					}
				}
				fieldsZ->omLocSYCL[q] = extractPoleToSYCL(gdata.om[q], gdata.mx[dir] + 6, dir, i, j, -3);
			}
      
			// Assigning shifted fields for magnetic induction
#if (FLUID_TYPE == CRONOS_MHD)
			for (int k = -3; k <= gdata.mx[dir]+2; ++k){
				fieldsZ->omLocP[q_Bx](k) = gdata.om[q_Bx](i,j,k);
			}
			for (int k = -3; k <= gdata.mx[dir]+2; ++k){
				fieldsZ->omLocM[q_Bx](k) = gdata.om[q_Bx](i-1,j,k);
			}
			for (int k = -3; k <= gdata.mx[dir]+2; ++k){
				fieldsZ->omLocP[q_By](k) = gdata.om[q_By](i,j,k);
			}
			for (int k = -3; k <= gdata.mx[dir]+2; ++k){
				fieldsZ->omLocM[q_By](k) = gdata.om[q_By](i,j-1,k);
			}
#endif

#if (OMS_USER == TRUE)
			for (int q = 0; q < n_omIntUser; ++q){
				for (int k = -3; k <= gdata.mx[dir]+2; ++k){
					fieldsZUser->omLoc[q](k) = gdata.om_user[q](i,j,k);
				}
			}
#endif

			// Assigning carbuncle flag values if necessary
			if(gdata.use_carbuncleFlag) {
				for (int k = -2; k <= gdata.mx[2]+2; ++k){
					fieldsZ->carbuncle_flag(k) = gdata.carbuncleFlag(i,j,k);
				}
			}

			// Doing the reconstruction:
			ReconstZ[n]->compute(queue, gdata, *fieldsZ, *physValzL, *physValzR);

#if (OMS_USER == TRUE)
			// Reconstruction for User-data:
			ReconstZUser[n]->compute(gdata, *fieldsZUser,
			                      *physValzLUser, *physValzRUser);
#endif
      
			// Get conservative variables:
			// on LHS
			Trafo->get_Cons(queue, gdata, Problem, *eos, ipos, *physValzL, dir, -0.5);
			// on RHS
			Trafo->get_Cons(queue, gdata, Problem, *eos, ipos, *physValzR, dir,  0.5);

#if (OMS_USER == TRUE)
			// Get conservative variables (needed for Riemann solver)
			// -- just by copying for the time being
			// LHS
			Problem.get_Cons(gdata, ipos, *physValzL, *physValzLUser, dir, -0.5);
			// Trafo->get_ConsUser(gdata, *physValzLUser, n_omIntUser);
			// RHS
			Problem.get_Cons(gdata, ipos, *physValzR, *physValzRUser, dir,  0.5);
			// Trafo->get_ConsUser(gdata, *physValzRUser, n_omIntUser);
#endif

			// Compute physical fluxes:
			// on LHS
			PhysFlux->get_PhysFlux(queue, gdata, Problem, ipos, *physValzL, dir, -0.5);
			// on RHS
			PhysFlux->get_PhysFlux(queue, gdata, Problem, ipos, *physValzR, dir,  0.5);


			// Possible modification by user
			// TODO PHILGS: Not used by current integration tests, disable for now
			//Problem.get_PhysFlux(gdata, ipos, *physValzL, dir, -0.5);

			// TODO PHILGS: Not used by current integration tests, disable for now
			//Problem.get_PhysFlux(gdata, ipos, *physValzR, dir,  0.5);

#if (OMS_USER == TRUE)
// 			PhysFluxUser->get_PhysFlux(gdata, Problem, ipos,
// 			                           *physValzL, *physValzLUser,
// 			                           dir, -0.5);
			Problem.get_PhysFluxUser(gdata, ipos, *physValzL, *physValzLUser,
			                         dir, -0.5);
			// on RHS
// 			PhysFluxUser->get_PhysFlux(gdata, Problem, ipos,
// 			                           *physValzR, *physValzRUser,
// 			                           dir,  0.5);
			Problem.get_PhysFluxUser(gdata, ipos, *physValzR, *physValzRUser,
			                         dir,  0.5);
#endif

			// Compute characteristic velocities:
			RiemannZ->get_vChar(queue, gdata, Problem, ipos, *physValzL, *physValzR,
						fieldsZ->v_ch_mORIG, fieldsZ->v_ch_pORIG, fieldsZ->v_ch_mSYCL, fieldsZ->v_ch_pSYCL, dir, cfl_lin);
//			RiemannZ->get_vChar(gdata, Problem, ipos, *physValzL, *physValzR, *fieldsZ,
//			                   dir, cfl_lin);

#if (OMS_USER == TRUE)
			fieldsZUser->v_ch_p = fieldsZ->v_ch_p;
			fieldsZUser->v_ch_m = fieldsZ->v_ch_m;
#endif

			// Compute numerical flux:
			RiemannZ->get_NumFlux(queue, gdata, *physValzL, *physValzR, *fieldsZ, dir);

#if (OMS_USER == TRUE)
			RiemannZUser->get_NumFlux(gdata, *physValzLUser, *physValzRUser,
			                          *fieldsZUser, dir);
#endif

			// Compute changes (noms)
			if((i>=0 && i<=gdata.mx[0]) && (j>=0 && j<=gdata.mx[1])) {
				get_Changes(gdata, fieldsZ->fluxORIG, gdata.nom, fieldsZ->ptotalORIG,
				            dir, i, j, n_omInt, RiemannZ->get_Fluid_Type());
			}

#if (OMS_USER == TRUE)
			if((i>=0 && i<=gdata.mx[0]) && (j>=0 && j<=gdata.mx[1])) {
				get_Changes(gdata, fieldsZUser->flux, gdata.nom_user,
				            fieldsZUser->ptotal,
				            dir, i, j, n_omIntUser,
				            RiemannZUser->get_Fluid_Type());
			}
			
#endif

			// Save 1D data if necessary
			Save->save_vars(gdata, *fieldsZ, dir, i, j);





		}
	}

// ---------------------------------------------------------------
//   Evolution for EMF via Constrained Transport
// ----------------------------------------------------------------

	gettimeofday(&tin5, 0);

#if (FLUID_TYPE == CRONOS_MHD)

#if (USE_COROTATION == CRONOS_ON)
	// Transform to co-rotating frame velocity for user src
	Trafo->TransInertToCorot(gdata, gfunc, Problem);
#endif

	if(gdata.tstep==0 && n==0) {
		if(IntegrateA) {
			gdata.nom[q_Bx].resize(Index::set(0,0,0),
					Index::set(gdata.mx[0],gdata.mx[1]+1,gdata.mx[2]+1));
			gdata.nom[q_By].resize(Index::set(0,0,0),
					Index::set(gdata.mx[0]+1,gdata.mx[1],gdata.mx[2]+1));
			gdata.nom[q_Bz].resize(Index::set(0,0,0),
					Index::set(gdata.mx[0]+1,gdata.mx[1]+1,gdata.mx[2]));
		} else {
			gdata.nom[q_Bx].resize(Index::set(-1,0,0),
					Index::set(gdata.mx[0],gdata.mx[1],gdata.mx[2]));
			gdata.nom[q_By].resize(Index::set(0,-1,0),
					Index::set(gdata.mx[0],gdata.mx[1],gdata.mx[2]));
			gdata.nom[q_Bz].resize(Index::set(0,0,-1),
					Index::set(gdata.mx[0],gdata.mx[1],gdata.mx[2]));
		}
	}
	gdata.nom[q_Bx].clear();
	gdata.nom[q_By].clear();
	gdata.nom[q_Bz].clear();

//#if (USE_COROTATION == CRONOS_ON)
//	// E_x:
//	CTSolveX->get_NumEMFCT(gdata, Problem, *Save, *Trafo, nom);
//
//
//	// Ey:
//	CTSolveY->get_NumEMFCT(gdata, Problem, *Save, *Trafo, nom);
//
//
//	// Ez:
//	CTSolveZ->get_NumEMFCT(gdata, Problem, *Save, *Trafo, nom);
//#else
	// E_x:
	CTSolveX->get_NumEMFCT(gdata, Problem, *Save, gdata.nom);


	// Ey:
	CTSolveY->get_NumEMFCT(gdata, Problem, *Save, gdata.nom);
	
	
	// Ez:
	CTSolveZ->get_NumEMFCT(gdata, Problem, *Save, gdata.nom);
//#endif

#if (USE_COROTATION == CRONOS_ON)
		// Transform back for energy trafo
		Trafo->TransCorotToInert(gdata, gfunc, Problem);
#endif

#endif


	delete Save;


	gettimeofday(&tin6, 0);

// ----------------------------------------------------------------
//   Check for errors:
// ----------------------------------------------------------------

	for(int q = 0; q<n_omInt; ++q) {
		CheckNan(gdata.nom[q],q, 0, 1,"nom");
	}

// ----------------------------------------------------------------
//   Compute Courant number
// ----------------------------------------------------------------

	REAL cfl = compute_cfl(gdata, Problem, cfl_eta, cfl_lin, n);
	
// ----------------------------------------------------------------
//   Geometrical source terms:
// ----------------------------------------------------------------

#if (GEOM != CARTESIAN)
	sources->src_Geom(gdata, Problem, gdata.nom);
#endif

// ----------------------------------------------------------------
//   Corotation source terms:
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

// 	if(ENERGETICS == FULL) {
// 		if(thermal) {
// 			Trafo->TransEth2E(gdata, gfunc, Problem);
// 		} else {
// 			Trafo->TransT2E(gdata, gfunc, Problem);
// 		}
// 	}  


// #if (USE_ANGULAR_MOMENTUM == TRUE)
// 	Trafo->TransVel2AngMom(gdata, gfunc, Problem);
// #else
// 	Trafo->TransVel2Momen(gdata, gfunc, Problem);
// #endif

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

	gettimeofday(&tock, 0);
	cstep = clock();

// ----------------------------------------------------------------
//   Determine time at intermediate steps:
// ----------------------------------------------------------------

	// Time must come before boundary for time dependent boundaries
	gdata.time += TimeIntegratorGeneric[0]->get_dt(gdata, n);


// ----------------------------------------------------------------
//   Transformation to base variables and bcs
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
//	for(int q=0; q<q_Eges; ++q) {
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
//   Performance estimate
// ----------------------------------------------------------------

	gettimeofday(&tock2, 0);
	cend = clock();

	timeval tstepOld = tstep;
	if(n == TIME_SUBSTEPS-1) {
		gettimeofday(&tstep, 0);
	}

	REAL delt = ((tock.tv_sec + tock.tv_usec/1.e6) - 
	             (tick.tv_sec + tick.tv_usec/1.e6));

	REAL delt2 = ((tock2.tv_sec + tock2.tv_usec/1.e6) - 
	              (tock.tv_sec + tock.tv_usec/1.e6));

	REAL dtStep = ((tstep.tv_sec + tstep.tv_usec/1.e6) -
	               (tstepOld.tv_sec + tstepOld.tv_usec/1.e6));


	if(gdata.rank == 0) {
		if(n == TIME_SUBSTEPS-1 && Problem.get_Info() && Problem.checkout(5)) {
			cout << "------------------------------------------------------" << endl;
			cout << " Time needed for substeps:  " << delt << " " << delt2 << endl;
			cout << "             for full step: " << dtStep << endl;
//       cout << " CPU Cycle times: " << gdata.rank << " ";
//       cout << (1.*(cstep-cstart))/CLOCKS_PER_SEC << " ";
//       cout << (1.*(cend-cstep))/CLOCKS_PER_SEC << endl;
		}

		// REAL deltat[15];

		// deltat[4] = ((tin6.tv_sec + tin6.tv_usec/1.e6) - 
		//              (tin5.tv_sec + tin5.tv_usec/1.e6));
		// deltat[3] = ((tin5.tv_sec + tin5.tv_usec/1.e6) - 
		//              (tin4.tv_sec + tin4.tv_usec/1.e6));
		// deltat[2] = ((tin4.tv_sec + tin4.tv_usec/1.e6) - 
		//              (tin3.tv_sec + tin3.tv_usec/1.e6));
		// deltat[1] = ((tin3.tv_sec + tin3.tv_usec/1.e6) - 
		//              (tin2.tv_sec + tin2.tv_usec/1.e6));
		// deltat[0] = ((tin2.tv_sec + tin2.tv_usec/1.e6) - 
		//              (tin1.tv_sec + tin1.tv_usec/1.e6));

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
	if(Problem.mag && n == TIME_SUBSTEPS-1){
		compute_divB(gdata, gfunc, Problem);
	}


	return cfl;



}
#endif
