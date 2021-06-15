#include "gridgen.H"
#include <iomanip>
#include <cmath>
#include <vector>
#include "DissipationMHD.H"

using namespace std;

// This file is only applicable for the multifluid case

#if(FLUID_TYPE == CRONOS_MULTIFLUID)

REAL HyperbolicSolver::singlestep(Data &gdata, gridFunc &gfunc,
		ProblemType &Problem, int n)
{

	if(eos == (EquationOfState*)0) {
		eos = new EquationOfState(Problem);
	}
	if(sources == NULL) {
		sources = new SourceTerms(gdata, Problem, IntegrateA, thermal);
	}


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

	// get total number of fields
	n_omInt = gdata.fluids->get_N_OMINT();
	n_Omega = gdata.fluids->get_N_OMEGA();

//	nom = new NumMatrix<REAL,3> [n_omInt];
//	for (int q = 0; q < n_omInt; ++q) {
//		if((gdata.om[q].getName() == "B_x" ||
//				gdata.om[q].getName() == "B_y" ||
//				gdata.om[q].getName() == "B_z")) {
//			nom[q].resize(Index::set(0,0,0),
//					Index::set(0,0,0));
//		} else {
//			nom[q].resize(Index::set(0,0,0),
//					Index::set(gdata.mx[0],gdata.mx[1],gdata.mx[2]));
//		}

	for (int q = 0; q < n_omInt; ++q) {
		gdata.nom[q].clear();
	}

#if (OMS_USER == TRUE)

	// get total number of user fields
	n_omIntUser = gdata.fluids->get_N_OMINT_USER();
	n_OmegaUser = gdata.fluids->get_N_OMEGA_USER();

//	nom_user = new NumMatrix<REAL,3> [n_omIntUser];
//	for (int q = 0; q < n_omIntUser; ++q) {
//		nom_user[q].resize(Index::set(0,0,0),
//				Index::set(gdata.mx[0],gdata.mx[1],gdata.mx[2]));
//		nom_user[q].clear();
//	}

	for (int q = 0; q < n_omIntUser; ++q) {
		gdata.nom_user[q].clear();
	}

#endif

	//	cout << " Boz " << gdata.om[6](0,399,0) << " " << gdata.om[6](0,400,0) << endl;

	//	cout << endl << " Start " << endl;
	//	cout << gdata.om[8](399,0,0) << " ";
	//	cout << gdata.om[9](399,0,0) << " ";
	//	cout << gdata.om[10](399,0,0) << " ";
	//	cout << gdata.om[11](399,0,0) << " ";
	//	cout << gdata.om[12](399,0,0) << " ";
	//	cout << endl;

	// ---------------------------------------------------------------
	//      Trafo of variables to conservative form
	//----------------------------------------------------------------

	// Loop over all fluids:
	for(int iFluid=0; iFluid<gdata.fluids->get_numFluids(); ++iFluid) {
		// Set indices for each fluid
		//		Trafo->reset_Indices(gdata.fluids->fluids[iFluid]);
#if(CRSWITCH_DUAL_ENERGY == CRONOS_ON)
		if(n_omInt > 7) {
			TrafoMulti[iFluid]->computeEntropyFromE(gdata, gfunc, Problem);
		}
#endif

		TrafoMulti[iFluid]->TransPrim2Cons(gdata, gfunc, Problem);

		// User defined transform
		Problem.TransPrim2Cons(gdata);
	}


#if(FLUID_TYPE == CRONOS_MULTIFLUID)
	bool has_MagField = gdata.fluids->with_magField();
#elif(FLUID_TYPE == CRONOS_MHD)
	bool has_MagField = true;
#else
	bool has_MagField = false;
#endif

	// ---------------------------------------------------------------
	//      Using primitive variables for reconstruction
	//----------------------------------------------------------------

	// Loop over all fluids:
	for(int iFluid=0; iFluid<gdata.fluids->get_numFluids(); ++iFluid) {
		// Set indices for each fluid
		//		Trafo->reset_Indices(gdata.fluids->fluids[iFluid]);

		TrafoMulti[iFluid]->TransCons2Prim(gdata, gfunc, Problem);

		// User defined transform
		Problem.TransCons2Prim(gdata);
	}



	gdata.CheckNan(1);


	REAL cfl_lin(0.), cfl_eta(0.);

	cronos::vector<REAL> pos(0,0,0);
	cronos::vector<REAL> ipos(0,0,0);

	//	cout << gdata.om[8](399,0,0) << " ";
	//	cout << gdata.om[9](399,0,0) << " ";
	//	cout << gdata.om[10](399,0,0) << " ";
	//	cout << gdata.om[11](399,0,0) << " ";
	//	cout << gdata.om[12](399,0,0) << " ";
	//	cout << endl;

	int dir=0;

	for(int k=-1; k<=gdata.mx[2]+1; ++k) {
		// Update of perp part of position:
		pos.set(2, gdata.getCen_z(k));
		ipos.set(2, k);
		for (int j=-1; j<=gdata.mx[1]+1; ++j){
			// Update of perp part of position:
			pos.set(1, gdata.getCen_y(j));
			ipos.set(1, j);

			// Assigning 1D array data
			// Loop over all data for multifluid case - fill up all local fields

			// Loop over all fluids:
			for(int iFluid=0; iFluid<numFluids; ++iFluid) { // Loop over all fluids
				// Store arrays for each fluid
				n_omInt = gdata.fluids->fluids[iFluid].get_N_OMINT();
				int has_MagFieldLoc = gdata.fluids->fluids[iFluid].has_MagField();

				for (int q = 0; q < n_omInt; ++q){

					if((has_MagFieldLoc && (q<q_By || q>q_Bz)) || !has_MagFieldLoc) {
						// Get global index:
						int q_global = gdata.fluids->fluids[iFluid].get_IndexGlobal(q);
						for (int i = -3; i <= gdata.mx[0]+2; ++i){
							fieldsMultiX[iFluid]->omLoc[q](i) = gdata.om[q_global](i,j,k);
						}
					}
				}
				//				if(j==0 && k==0 && iFluid==0) {
				//					cout << " Werte " << iFluid << "  - ";
				//					cout << fieldsMultiX[iFluid]->omLoc[0](401) << " ";
				//					cout << fieldsMultiX[iFluid]->omLoc[0](400) << " ";
				//					cout << fieldsMultiX[iFluid]->omLoc[0](399) << " ";
				//					cout << endl;
				//					cout << fieldsMultiX[iFluid]->omLoc[7](401) << " ";
				//					cout << fieldsMultiX[iFluid]->omLoc[7](400) << " ";
				//					cout << fieldsMultiX[iFluid]->omLoc[7](399) << " ";
				//					cout << endl;
				//				}

				// Assigning shifted fields for magnetic induction
				if(has_MagFieldLoc) {
					int q_By_global = gdata.fluids->fluids[iFluid].get_q_By_global();
					int q_Bz_global = gdata.fluids->fluids[iFluid].get_q_Bz_global();
					for (int i = -3; i <= gdata.mx[0]+2; ++i){
						fieldsMultiX[iFluid]->omLocP[q_By](i) = gdata.om[q_By_global](i,j,k);
					}
					for (int i = -3; i <= gdata.mx[0]+2; ++i){
						fieldsMultiX[iFluid]->omLocM[q_By](i) = gdata.om[q_By_global](i,j-1,k);
					}
					for (int i = -3; i <= gdata.mx[0]+2; ++i){
						fieldsMultiX[iFluid]->omLocP[q_Bz](i) = gdata.om[q_Bz_global](i,j,k);
					}
					for (int i = -3; i <= gdata.mx[0]+2; ++i){
						fieldsMultiX[iFluid]->omLocM[q_Bz](i) = gdata.om[q_Bz_global](i,j,k-1);
					}
				}

				// Doing the reconstruction:
				ReconstMultiX[iFluid][n]->compute(gdata, *fieldsMultiX[iFluid],
						*physValMultixL[iFluid], *physValMultixR[iFluid]);

				//				if(j==0 && k==0 && iFluid==1) {
				//					cout << " uPri " << iFluid << "  - ";
				//					cout << physValMultixL[iFluid]->uPri[0](401) << " ";
				//					cout << physValMultixL[iFluid]->uPri[0](400) << " ";
				//					cout << physValMultixL[iFluid]->uPri[0](399) << " ";
				//					cout << endl;
				//				}

				// Get conservative variables:
				// on LHS
				TrafoMulti[iFluid]->get_Cons(gdata, Problem, *eos, ipos, *physValMultixL[iFluid], dir, -0.5);
				// on RHS
				TrafoMulti[iFluid]->get_Cons(gdata, Problem, *eos, ipos, *physValMultixR[iFluid], dir,  0.5);

				//				if(j==0 && k==0 && iFluid==1) {
				//					cout << " uCon" << iFluid << "  - ";
				////					cout << "     ";
				////					cout << physValMultixL[iFluid]->uCon[0](401) << " ";
				////					cout << physValMultixL[iFluid]->uCon[0](400) << " ";
				////					cout << physValMultixL[iFluid]->uCon[0](399) << " ";
				////					cout << endl;
				////					cout << "     ";
				////					cout << physValMultixL[iFluid]->uCon[1](401) << " ";
				////					cout << physValMultixL[iFluid]->uCon[1](400) << " ";
				////					cout << physValMultixL[iFluid]->uCon[1](399) << " ";
				////					cout << endl;
				////					cout << "     ";
				////					cout << physValMultixL[iFluid]->uCon[2](401) << " ";
				////					cout << physValMultixL[iFluid]->uCon[2](400) << " ";
				////					cout << physValMultixL[iFluid]->uCon[2](399) << " ";
				////					cout << endl;
				////					cout << "     ";
				////					cout << physValMultixL[iFluid]->uCon[7](401) << " ";
				////					cout << physValMultixL[iFluid]->uCon[7](400) << " ";
				////					cout << physValMultixL[iFluid]->uCon[7](399) << " ";
				////					cout << endl;
				////					cout << "     ";
				//					cout << physValMultixL[iFluid]->uCon[4](401) << " ";
				//					cout << physValMultixL[iFluid]->uCon[4](400) << " ";
				//					cout << physValMultixL[iFluid]->uCon[4](399) << " ";
				//					cout << endl;
				////					cout << "     ";
				////					cout << physValMultixL[iFluid]->uCon[5](401) << " ";
				////					cout << physValMultixL[iFluid]->uCon[5](400) << " ";
				////					cout << physValMultixL[iFluid]->uCon[5](399) << " ";
				////					cout << endl;
				////					cout << "     ";
				////					cout << physValMultixL[iFluid]->uCon[6](401) << " ";
				////					cout << physValMultixL[iFluid]->uCon[6](400) << " ";
				////					cout << physValMultixL[iFluid]->uCon[6](399) << " ";
				////					cout << endl;
				//
				//				}

				// Compute physical fluxes:
				// on LHS
				PhysFluxMulti[iFluid]->get_PhysFlux(gdata, Problem, ipos, *physValMultixL[iFluid], dir, -0.5);
				// on RHS
				PhysFluxMulti[iFluid]->get_PhysFlux(gdata, Problem, ipos, *physValMultixR[iFluid], dir,  0.5);

				//				if(j==0 && k==0 && iFluid==1) {
				//					cout << " pFlux: " << endl;
				//					cout << physValMultixL[iFluid]->flux[4](399) << " ";
				//					cout << physValMultixL[iFluid]->flux[4](400) << " ";
				//					cout << physValMultixL[iFluid]->flux[4](401) << " ";
				//					cout << endl;
				//				}

				// Compute characteristic velocities:
				RiemannSolversX[iFluid]->get_vChar(gdata, Problem, ipos, *physValMultixL[iFluid], *physValMultixR[iFluid],
						fieldsMultiX[iFluid]->v_ch_m, fieldsMultiX[iFluid]->v_ch_p, dir, cfl_lin);

				//				if(iFluid==1 && j==0 && k==0) {
				//					cout << " vChar " << endl;
				//					cout << fieldsMultiX[iFluid]->v_ch_m(399) << " ";
				//					cout << fieldsMultiX[iFluid]->v_ch_m(400) << " ";
				//					cout << fieldsMultiX[iFluid]->v_ch_m(401) << " ";
				//					cout << endl;
				//				}

				// Compute numerical fluxes:
				RiemannSolversX[iFluid]->get_NumFlux(gdata,
						*physValMultixL[iFluid], *physValMultixR[iFluid], *fieldsMultiX[iFluid], dir);

				// Changes are computed for each field (not fluid) individually
				if((j>=0 && j<=gdata.mx[1]) && (k>=0 && k<=gdata.mx[2])) {

					//					if(iFluid==1) {
					//					cout << " flu ";
					//					cout << fieldsMultiX[iFluid]->flux[4](402) << " ";
					//					cout << fieldsMultiX[iFluid]->flux[4](401) << " ";
					//					cout << fieldsMultiX[iFluid]->flux[4](400) << " ";
					//					cout << fieldsMultiX[iFluid]->flux[4](399) << " ";
					//					cout << fieldsMultiX[iFluid]->flux[4](398) << " ";
					//					cout << endl;
					//					}
					get_Changes(gdata, fieldsMultiX[iFluid]->flux, gdata.nom, fieldsMultiX[iFluid]->ptotal,
							0, j, k, n_omInt, RiemannSolversX[iFluid]->get_Fluid_Type(), iFluid);

				}
				if(has_MagFieldLoc) {
					Save->save_vars(gdata, *fieldsMultiX[iFluid], 0, j, k);
				}
			}
#if (OMS_USER == TRUE)
			// Assing 1D array data for user fields:

			// Loop over all fluids:
			for(int iFluid=0; iFluid<numFluids; ++iFluid) { // Loop over all fluids

				n_omIntUser = gdata.fluids->fluids[iFluid].get_N_OMINT_USER();

				for (int q = 0; q < n_omIntUser; ++q){
					// Get global index:
					int q_global = gdata.fluids->fluids[iFluid].get_IndexGlobalUser(q);

					for (int i = -3; i <= gdata.mx[0]+2; ++i){
						fieldsXUser->omLoc[q](i) = gdata.om_user[q_global](i,j,k);
					}
				}

				// Reconstruction for User-data:
				ReconstXUser[n]->compute(gdata, *fieldsXUser,
						*physValxLUser, *physValxRUser);

				// Get conservative variables (needed for Riemann solver)
				// -- default is just a copy
				// LHS
				Problem.get_Cons(gdata, ipos, *physValxL, *physValxLUser, dir, -0.5);
				// Trafo->get_ConsUser(gdata, *physValxLUser, n_omIntUser);
				// RHS
				Problem.get_Cons(gdata, ipos, *physValxR, *physValxRUser, dir,  0.5);
				// Trafo->get_ConsUser(gdata, *physValxRUser, n_omIntUser);

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

				fieldsXUser->v_ch_p = fieldsX->v_ch_p;
				fieldsXUser->v_ch_m = fieldsX->v_ch_m;

				RiemannXUser->get_NumFlux(gdata, *physValxLUser, *physValxRUser,
						*fieldsXUser, dir);


				if((j>=0 && j<=gdata.mx[1]) && (k>=0 && k<=gdata.mx[2])) {
					get_Changes(gdata, fieldsXUser->flux, gdata.nom_user,
							fieldsXUser->ptotal,
							0, j, k, n_omIntUser,
							RiemannXUser->get_Fluid_Type());
				}
			}
#endif



		}
	}

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
			for(int iFluid=0; iFluid<gdata.fluids->get_numFluids(); ++iFluid) { // Loop over all fluids

				// Store arrays for each fluid
				n_omInt = gdata.fluids->fluids[iFluid].get_N_OMINT();
				int has_MagFieldLoc = gdata.fluids->fluids[iFluid].has_MagField();

				for (int q = 0; q < n_omInt; ++q){
					if((has_MagFieldLoc && q!=q_Bx && q!=q_Bz) || !has_MagFieldLoc) {
						// Get global index:
						int q_global = gdata.fluids->fluids[iFluid].get_IndexGlobal(q);
						for (int j = -3; j <= gdata.mx[1]+2; ++j){
							fieldsMultiY[iFluid]->omLoc[q](j) = gdata.om[q_global](i,j,k);
						}
					}
				}


				// Assigning shifted fields for magnetic induction
				if(has_MagFieldLoc) {
					int q_Bx_global = gdata.fluids->fluids[iFluid].get_q_Bx_global();
					int q_Bz_global = gdata.fluids->fluids[iFluid].get_q_Bz_global();
					for (int j = -3; j <= gdata.mx[1]+2; ++j){
						fieldsMultiY[iFluid]->omLocP[q_Bx](j) = gdata.om[q_Bx_global](i,j,k);
					}
					for (int j = -3; j <= gdata.mx[1]+2; ++j){
						fieldsMultiY[iFluid]->omLocM[q_Bx](j) = gdata.om[q_Bx_global](i-1,j,k);
					}
					for (int j = -3; j <= gdata.mx[1]+2; ++j){
						fieldsMultiY[iFluid]->omLocP[q_Bz](j) = gdata.om[q_Bz_global](i,j,k);
					}
					for (int j = -3; j <= gdata.mx[1]+2; ++j){
						fieldsMultiY[iFluid]->omLocM[q_Bz](j) = gdata.om[q_Bz_global](i,j,k-1);
					}
				}

				//				if(i==0 && k==0 && iFluid==0) {
				//					cout << " Werte: " << q_Bx << " " << q_By << " " << q_Bz << endl;
				//					cout << fieldsMultiY[iFluid]->omLoc[0](399) << " ";
				//					cout << fieldsMultiY[iFluid]->omLoc[0](400) << " ";
				//					cout << fieldsMultiY[iFluid]->omLoc[0](401) << " ";
				//					cout << endl;
				//					cout << fieldsMultiY[iFluid]->omLocP[4](399) << " ";
				//					cout << fieldsMultiY[iFluid]->omLocP[4](400) << " ";
				//					cout << fieldsMultiY[iFluid]->omLocP[4](401) << " ";
				//					cout << endl;
				//					cout << fieldsMultiY[iFluid]->omLoc[5](399) << " ";
				//					cout << fieldsMultiY[iFluid]->omLoc[5](400) << " ";
				//					cout << fieldsMultiY[iFluid]->omLoc[5](401) << " ";
				//					cout << endl;
				//					cout << fieldsMultiY[iFluid]->omLocP[6](399) << " ";
				//					cout << fieldsMultiY[iFluid]->omLocP[6](400) << " ";
				//					cout << fieldsMultiY[iFluid]->omLocP[6](401) << " ";
				//					cout << endl;
				//					cout << fieldsMultiY[iFluid]->omLoc[7](399) << " ";
				//					cout << fieldsMultiY[iFluid]->omLoc[7](400) << " ";
				//					cout << fieldsMultiY[iFluid]->omLoc[7](401) << " ";
				//					cout << endl;
				//					cout << endl << " Bz: " << endl;
				//					cout << gdata.om[6](i,390,k) << " ";
				//					cout << gdata.om[6](i,399,k) << " ";
				//					cout << gdata.om[6](i,400,k) << " ";
				//					cout << gdata.om[6](i,401,k) << " ";
				//					cout << endl;
				//
				//				}

				// Doing the reconstruction:
				ReconstMultiY[iFluid][n]->compute(gdata, *fieldsMultiY[iFluid],
						*physValMultiyL[iFluid], *physValMultiyR[iFluid]);


				// Get conservative variables:
				// on LHS
				TrafoMulti[iFluid]->get_Cons(gdata, Problem, *eos, ipos, *physValMultiyL[iFluid], dir, -0.5);
				// on RHS
				TrafoMulti[iFluid]->get_Cons(gdata, Problem, *eos, ipos, *physValMultiyR[iFluid], dir,  0.5);


				// Compute physical fluxes:
				// on LHS
				PhysFluxMulti[iFluid]->get_PhysFlux(gdata, Problem, ipos, *physValMultiyL[iFluid], dir, -0.5);
				// on RHS
				PhysFluxMulti[iFluid]->get_PhysFlux(gdata, Problem, ipos, *physValMultiyR[iFluid], dir,  0.5);


				// Compute characteristic velocities:
				RiemannSolversY[iFluid]->get_vChar(gdata, Problem, ipos,
						*physValMultiyL[iFluid], *physValMultiyR[iFluid],
						fieldsMultiY[iFluid]->v_ch_m, fieldsMultiY[iFluid]->v_ch_p, dir, cfl_lin);

				// Compute numerical flux:
				RiemannSolversY[iFluid]->get_NumFlux(gdata, *physValMultiyL[iFluid],
						*physValMultiyR[iFluid], *fieldsMultiY[iFluid], dir);


				// Compute changes (noms)
				if((i>=0 && i<=gdata.mx[0]) && (k>=0 && k<=gdata.mx[2])) {
					get_Changes(gdata, fieldsMultiY[iFluid]->flux, gdata.nom, fieldsMultiY[iFluid]->ptotal,
							dir, i, k, n_omInt, RiemannSolversY[iFluid]->get_Fluid_Type(), iFluid);
				}


				if(has_MagFieldLoc) {
					// Save 1D data if necessary
					Save->save_vars(gdata, *fieldsMultiY[iFluid], dir, i, k);
				}
			}



#if (OMS_USER == TRUE)
			for (int q = 0; q < n_omIntUser; ++q){
				for (int j = -3; j <= gdata.mx[dir]+2; ++j){
					fieldsYUser->omLoc[q](j) = gdata.om_user[q](i,j,k);
				}
			}


			// Reconstruction for User-data:
			ReconstYUser[n]->compute(gdata, *fieldsYUser,
					*physValyLUser, *physValyRUser);


			// Get conservative variables (needed for Riemann solver)
			// -- just by copying for the time being
			// LHS
			Problem.get_Cons(gdata, ipos, *physValyL, *physValyLUser, dir, -0.5);
			// Trafo->get_ConsUser(gdata, *physValyLUser, n_omIntUser);
			// RHS
			Problem.get_Cons(gdata, ipos, *physValyR, *physValyRUser, dir, 0.5);
			// Trafo->get_ConsUser(gdata, *physValyRUser, n_omIntUser);


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


			fieldsYUser->v_ch_p = fieldsY->v_ch_p;
			fieldsYUser->v_ch_m = fieldsY->v_ch_m;


			RiemannYUser->get_NumFlux(gdata, *physValyLUser, *physValyRUser,
					*fieldsYUser, dir);


			if((i>=0 && i<=gdata.mx[0]) && (k>=0 && k<=gdata.mx[2])) {

				get_Changes(gdata, fieldsYUser->flux, gdata.nom_user,
						fieldsYUser->ptotal,
						dir, i, k, n_omIntUser,
						RiemannYUser->get_Fluid_Type());
			}

#endif


		}
	}



	// z-direction

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

			// Each of the following is to be done for each fluid individually
			for(int iFluid=0; iFluid<gdata.fluids->get_numFluids(); ++iFluid) { // Loop over all fluids

				n_omInt = gdata.fluids->fluids[iFluid].get_N_OMINT();
				int has_MagFieldLoc = gdata.fluids->fluids[iFluid].has_MagField();

				for (int q = 0; q < n_omInt; ++q){
					if((has_MagFieldLoc && (q<q_Bx || q>q_By)) || !has_MagFieldLoc) {
						// Get global index:
						int q_global = gdata.fluids->fluids[iFluid].get_IndexGlobal(q);
						for (int k = -3; k <= gdata.mx[2]+2; ++k){
							fieldsMultiZ[iFluid]->omLoc[q](k) = gdata.om[q_global](i,j,k);
						}
					}
				}

				// Assigning shifted fields for magnetic induction
				if(has_MagFieldLoc) {
					int q_Bx_global = gdata.fluids->fluids[iFluid].get_q_Bx_global();
					int q_By_global = gdata.fluids->fluids[iFluid].get_q_By_global();
					for (int k = -3; k <= gdata.mx[dir]+2; ++k){
						fieldsMultiZ[iFluid]->omLocP[q_Bx](k) = gdata.om[q_Bx_global](i,j,k);
					}
					for (int k = -3; k <= gdata.mx[dir]+2; ++k){
						fieldsMultiZ[iFluid]->omLocM[q_Bx](k) = gdata.om[q_Bx_global](i-1,j,k);
					}
					for (int k = -3; k <= gdata.mx[dir]+2; ++k){
						fieldsMultiZ[iFluid]->omLocP[q_By](k) = gdata.om[q_By_global](i,j,k);
					}
					for (int k = -3; k <= gdata.mx[dir]+2; ++k){
						fieldsMultiZ[iFluid]->omLocM[q_By](k) = gdata.om[q_By_global](i,j-1,k);
					}
				}

				// Doing the reconstruction:
				ReconstMultiZ[iFluid][n]->compute(gdata, *fieldsMultiZ[iFluid],
						*physValMultizL[iFluid], *physValMultizR[iFluid]);


				// Get conservative variables:
				// on LHS
				TrafoMulti[iFluid]->get_Cons(gdata, Problem, *eos, ipos, *physValMultizL[iFluid], dir, -0.5);
				// on RHS
				TrafoMulti[iFluid]->get_Cons(gdata, Problem, *eos, ipos, *physValMultizR[iFluid], dir,  0.5);


				// Compute physical fluxes:
				// on LHS
				PhysFluxMulti[iFluid]->get_PhysFlux(gdata, Problem, ipos, *physValMultizL[iFluid], dir, -0.5);
				// on RHS
				PhysFluxMulti[iFluid]->get_PhysFlux(gdata, Problem, ipos, *physValMultizR[iFluid], dir,  0.5);



				// Compute characteristic velocities:
				RiemannSolversZ[iFluid]->get_vChar(gdata, Problem, ipos,
						*physValMultizL[iFluid], *physValMultizR[iFluid],
						fieldsMultiZ[iFluid]->v_ch_m, fieldsMultiZ[iFluid]->v_ch_p, dir, cfl_lin);
				//			RiemannZ->get_vChar(gdata, Problem, ipos, *physValzL, *physValzR, *fieldsZ,
				//			                   dir, cfl_lin);


				// Compute numerical flux:
				RiemannSolversZ[iFluid]->get_NumFlux(gdata, *physValMultizL[iFluid],
						*physValMultizR[iFluid], *fieldsMultiZ[iFluid], dir);

				// Compute changes (noms)
				if((i>=0 && i<=gdata.mx[0]) && (j>=0 && j<=gdata.mx[1])) {
					get_Changes(gdata, fieldsMultiZ[iFluid]->flux, gdata.nom, fieldsMultiZ[iFluid]->ptotal,
							dir, i, j, n_omInt, RiemannSolversZ[iFluid]->get_Fluid_Type(), iFluid);
				}

				// Save 1D data if necessary
				if(has_MagFieldLoc) {
					Save->save_vars(gdata, *fieldsMultiZ[iFluid], dir, i, j);
				}

			}


#if (OMS_USER == TRUE)
			for (int q = 0; q < n_omIntUser; ++q){
				for (int k = -3; k <= gdata.mx[dir]+2; ++k){
					fieldsZUser->omLoc[q](k) = gdata.om_user[q](i,j,k);
				}
			}


			// Reconstruction for User-data:
			ReconstZUser[n]->compute(gdata, *fieldsZUser,
					*physValzLUser, *physValzRUser);


			// Get conservative variables (needed for Riemann solver)
			// -- just by copying for the time being
			// LHS
			Problem.get_Cons(gdata, ipos, *physValzL, *physValzLUser, dir, -0.5);
			// Trafo->get_ConsUser(gdata, *physValzLUser, n_omIntUser);
			// RHS
			Problem.get_Cons(gdata, ipos, *physValzR, *physValzRUser, dir,  0.5);
			// Trafo->get_ConsUser(gdata, *physValzRUser, n_omIntUser);




			// Compute physical fluxes for user fields:
			// on LHS
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


			fieldsZUser->v_ch_p = fieldsZ->v_ch_p;
			fieldsZUser->v_ch_m = fieldsZ->v_ch_m;


			RiemannZUser->get_NumFlux(gdata, *physValzLUser, *physValzRUser,
					*fieldsZUser, dir);


			if((i>=0 && i<=gdata.mx[0]) && (j>=0 && j<=gdata.mx[1])) {
				get_Changes(gdata, fieldsZUser->flux, gdata.nom_user,
						fieldsZUser->ptotal,
						dir, i, j, n_omIntUser,
						RiemannZUser->get_Fluid_Type());
			}

#endif


		}
	}

	// ---------------------------------------------------------------
	//   Evolution for EMF via Constrained Transport
	// ----------------------------------------------------------------

	gettimeofday(&tin5, 0);


	if(has_MagField) {

#if(FLUID_TYPE == CRONOS_MULTIFLUID)
		int iMagFluid = gdata.fluids->get_i_magFluid();
		int q_Bx = gdata.fluids->fluids[iMagFluid].get_q_Bx_global();
		int q_By = gdata.fluids->fluids[iMagFluid].get_q_By_global();
		int q_Bz = gdata.fluids->fluids[iMagFluid].get_q_Bz_global();
		int q_sx = gdata.fluids->fluids[iMagFluid].get_q_sx_global();
		int q_sy = gdata.fluids->fluids[iMagFluid].get_q_sy_global();
		int q_sz = gdata.fluids->fluids[iMagFluid].get_q_sz_global();
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


		// E_x:
		CTSolveX->get_NumEMFCT(gdata, Problem, *Save, gdata.nom);


		// Ey:
		CTSolveY->get_NumEMFCT(gdata, Problem, *Save, gdata.nom);


		// Ez:
		CTSolveZ->get_NumEMFCT(gdata, Problem, *Save, gdata.nom);

	}

	delete Save;


	gettimeofday(&tin6, 0);

	// Compute total number of fields for multifluid case
#if(FLUID_TYPE==CRONOS_MULTIFLUID)
	n_omInt = gdata.fluids->get_N_OMINT();
#endif

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
	//   User defined source terms:
	// ----------------------------------------------------------------

	Problem.src_User(gdata, gdata.nom, gdata.nom_user);

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

	// Loop over all fluids:
	for(int iFluid=0; iFluid<gdata.fluids->get_numFluids(); ++iFluid) {
		// Set indices for each fluid
		//		Trafo->reset_Indices(gdata.fluids->fluids[iFluid]);

		TrafoMulti[iFluid]->TransPrim2Cons(gdata, gfunc, Problem);

		// User defined transform
		Problem.TransPrim2Cons(gdata);
	}

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



	Pot & Az = gdata.om[n_omIntAll+N_ADD+2];
	//	cout << " Az: " << Az(400,1,0) << " " << nom[6](400,1,0) << " " << n_omIntAll << endl;

	//	if(gdata.time > 0.01) exit(3);
	//	if(n==1 && gdata.time>gdata.dt)
	//	exit(3);

//	delete [] nom;
//#if (OMS_USER == TRUE)
//	delete [] nom_user;
//#endif

	gettimeofday(&tock, 0);
	cstep = clock();

	// ----------------------------------------------------------------
	//   Determine time at intermediate steps:
	// ----------------------------------------------------------------
	gdata.time += TimeIntegratorGeneric[0]->get_dt(gdata, n);

	// ----------------------------------------------------------------
	//   Transformation to base variables and bcs
	// ----------------------------------------------------------------

	gettimeofday(&tock2, 0);

	// Loop over all fluids:
	for(int iFluid=0; iFluid<gdata.fluids->get_numFluids(); ++iFluid) {
		// Set indices for each fluid
		//		Trafo->reset_Indices(gdata.fluids->fluids[iFluid]);
#if (USE_ANGULAR_MOMENTUM == TRUE)
		TrafoMulti[iFluid]->TransAngMom2Vel(gdata, gfunc, Problem);
#else
		TrafoMulti[iFluid]->TransMomen2Vel(gdata, gfunc, Problem);
#endif
	}

	//	cout << " Machfeld " << endl;
	//	cout << gdata.om[0](399,0,0) << " ";
	//	cout << gdata.om[1](399,0,0) << " ";
	//	cout << gdata.om[2](399,0,0) << " ";
	//	cout << gdata.om[3](399,0,0) << " ";
	//	cout << gdata.om[4](399,0,0) << " ";
	//	cout << gdata.om[5](399,0,0) << " ";
	//	cout << gdata.om[6](399,0,0) << " ";
	//	cout << gdata.om[7](399,0,0) << " ";
	//	cout << endl;

#if (FLUID_TYPE != CRONOS_HYDRO)
	if(has_MagField) {
		if(IntegrateA) {
			compute_B(gdata, gfunc, Problem);
		}
	}
#endif
	//	cout << gdata.om[0](399,0,0) << " ";
	//	cout << gdata.om[1](399,0,0) << " ";
	//	cout << gdata.om[2](399,0,0) << " ";
	//	cout << gdata.om[3](399,0,0) << " ";
	//	cout << gdata.om[4](399,0,0) << " ";
	//	cout << gdata.om[5](399,0,0) << " ";
	//	cout << gdata.om[6](399,0,0) << " ";
	//	cout << gdata.om[7](399,0,0) << " ";
	//	cout << endl;

	//	cout << gdata.om[8](399,0,0) << " ";
	//	cout << gdata.om[9](399,0,0) << " ";
	//	cout << gdata.om[10](399,0,0) << " ";
	//	cout << gdata.om[11](399,0,0) << " ";
	//	cout << gdata.om[12](399,0,0) << " ";
	//	cout << endl;

	// Loop over all fluids:
	for(int iFluid=0; iFluid<gdata.fluids->get_numFluids(); ++iFluid) {
		// Get range of indices
		int q_min = gdata.fluids->fluids[iFluid].get_IndexGlobal(0);
#if(ENERGETICS == FULL)
		int q_Eges = gdata.fluids->fluids[iFluid].get_q_Eges_global();
#else
		int q_Eges = gdata.fluids->fluids[iFluid].get_N_OMINT() + q_min;
#endif
		n_omInt = gdata.fluids->get_N_OMINT(iFluid);

		for(int q=q_min; q<q_Eges; ++q) {
			gfunc.boundary(gdata, Problem, gdata.om[q],B,q, iFluid);
		}


		if(ENERGETICS == FULL) {
			// Set indices for each fluid
			//			Trafo->reset_Indices(gdata.fluids->fluids[iFluid]);
			if(thermal) {
#if(CRSWITCH_DUAL_ENERGY == CRONOS_ON)
				TrafoMulti[iFluid]->TransE2Eth(gdata, gfunc, Problem, n, true);
#else
				TrafoMulti[iFluid]->TransE2Eth(gdata, gfunc, Problem);
#endif
			} else {
				TrafoMulti[iFluid]->TransE2T(gdata, gfunc, Problem);
			}

			for(int q=q_Eges; q<n_omInt; ++q) {
				gfunc.boundary(gdata, Problem, gdata.om[q],B,q, iFluid);
			}
		}

		Problem.TransCons2Prim(gdata);

		//		cout << " So ist es: " << iFluid;
		//		cout << gdata.om[1].getName() << " ";
		//		cout << gdata.om[2].getName() << " ";
		//		cout << gdata.om[3].getName() << " ";
		//		cout << gdata.om[9].getName() << " ";
		//		cout << gdata.om[10].getName() << " ";
		//		cout << gdata.om[11].getName() << " ";
		//		cout << endl << endl;
	}

	//	cout << " So ist es: ";
	//	cout << gdata.om[1].getName() << " ";
	//	cout << gdata.om[2].getName() << " ";
	//	cout << gdata.om[3].getName() << " ";
	//	cout << gdata.om[9].getName() << " ";
	//	cout << gdata.om[10].getName() << " ";
	//	cout << gdata.om[11].getName() << " ";
	//	cout << endl << endl;

	//	cout << gdata.om[8](399,0,0) << " ";
	//	cout << gdata.om[9](399,0,0) << " ";
	//	cout << gdata.om[10](399,0,0) << " ";
	//	cout << gdata.om[11](399,0,0) << " ";
	//	cout << gdata.om[12](399,0,0) << " ";
	//	cout << endl;
	//	cout << gdata.om[8](400,0,0) << " ";
	//	cout << gdata.om[9](400,0,0) << " ";
	//	cout << gdata.om[10](400,0,0) << " ";
	//	cout << gdata.om[11](400,0,0) << " ";
	//	cout << gdata.om[12](400,0,0) << " ";
	//	cout << endl << endl << endl;

#if (OMS_USER == TRUE)

	for(int iFluid=0; iFluid<gdata.fluids->get_numFluids(); ++iFluid) {
		for(int q=0; q<gdata.fluids->fluids[iFluid].get_N_OMINT_USER(); ++q) {
			int q_global = gdata.fluids->fluids[iFluid].get_IndexGlobal(q);
			gfunc.boundary(gdata, Problem, gdata.om_user[q_global],B,q_global, iFluid);
		}
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

	phystest(gdata, gfunc, Problem, n);

	//	cout << " So ist es nach phystest: ";
	//	cout << gdata.om[1].getName() << " ";
	//	cout << gdata.om[2].getName() << " ";
	//	cout << gdata.om[3].getName() << " ";
	//	cout << gdata.om[9].getName() << " ";
	//	cout << gdata.om[10].getName() << " ";
	//	cout << gdata.om[11].getName() << " ";
	//	cout << endl << endl;

	// Computing div B only at end of timestep
	if(Problem.mag && n == TIME_SUBSTEPS-1){
		compute_divB(gdata, gfunc, Problem);
	}


	//	cout << " So ist es: ";
	//	cout << gdata.om[1].getName() << " ";
	//	cout << gdata.om[2].getName() << " ";
	//	cout << gdata.om[3].getName() << " ";
	//	cout << gdata.om[9].getName() << " ";
	//	cout << gdata.om[10].getName() << " ";
	//	cout << gdata.om[11].getName() << " ";
	//	cout << endl << endl;
	//
	//	exit(3);

	return cfl;



}

#endif // Check wether we have multifluid case
