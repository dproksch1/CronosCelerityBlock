#include "gridgen.H"
#include "specific.H"
#include <iomanip>

#include "RiemannSolverHD.H"
#include "RiemannSolverMHD.H"
#include "PhysFluxesHD.H"
#include "PhysFluxesMHD.H"
#include "CTLondrilloDelZanna.H"
//#include "CTStone.H"
#include "timewrapper.H"

HyperbolicSolver::HyperbolicSolver(Data &gdata, ProblemType &Problem)
	: Ekin_tave(0.),
	  Emag_tave(0.),
	  Etherm_tave(0.),
	  Eges_tave(0.),
	  Ekfluc_tave(0.),
	  Ebfluc_tave(0.),
	  Vorticity_tave(0.),
	  Enstrophy_tave(0.),
	  cs_tave(0.),
	  v2ave_tave(0.),
	  time_ave(0.)
{

	// Copying some constants from Problem-class to use locally
//	this->q_rho = Problem.q_rho;
//	this->q_sx = Problem.q_sx;
//	this->q_sy = Problem.q_sy;
//	this->q_sz = Problem.q_sz;
//	this->q_Bx = Problem.q_Bx;
//	this->q_By = Problem.q_By;
//	this->q_Bz = Problem.q_Bz;
//	this->q_Eges = Problem.q_Eges;
//	this->q_Eadd = Problem.q_Eadd;

#if (FLUID_TYPE == CRONOS_MULTIFLUID)
	with_mag = gdata.fluids->with_magField();
	q_rho = gdata.fluids->get_q_rho(0);
	q_sx = gdata.fluids->get_q_sx(0);
	q_sy = gdata.fluids->get_q_sy(0);
	q_sz = gdata.fluids->get_q_sz(0);
	q_Bx = gdata.fluids->get_q_Bx();
	q_By = gdata.fluids->get_q_By();
	q_Bz = gdata.fluids->get_q_Bz();
	q_Eges = gdata.fluids->get_q_Eges(0);
	q_Eadd = gdata.fluids->get_q_Eadd(0);
#else
	with_mag = (FLUID_TYPE == CRONOS_MHD);
	q_rho = gdata.fluid.get_q_rho();
	q_sx = gdata.fluid.get_q_sx();
	q_sy = gdata.fluid.get_q_sy();
	q_sz = gdata.fluid.get_q_sz();
	q_Bx = gdata.fluid.get_q_Bx();
	q_By = gdata.fluid.get_q_By();
	q_Bz = gdata.fluid.get_q_Bz();
	q_Eges = gdata.fluid.get_q_Eges();
	q_Eadd = gdata.fluid.get_q_Eadd();
#endif


#if (FLUID_TYPE == CRONOS_MULTIFLUID)
	numFluids = gdata.fluids->get_numFluids();
	n_om = gdata.fluids->get_N_OM();
	n_omInt = gdata.fluids->get_N_OMINT();
	n_Omega = gdata.fluids->get_N_OMEGA();
	n_omUser = gdata.fluids->get_N_OM_USER();
	n_omIntUser = gdata.fluids->get_N_OMINT_USER();
	n_OmegaUser = gdata.fluids->get_N_OMEGA_USER();
	n_omIntAll = gdata.fluids->get_N_OMINT_ALL();
#else
	n_om = N_OM;
	n_omInt = N_OMINT;
	n_Omega = N_OMEGA;
#if(OMS_USER == TRUE)
	n_omUser = N_OM_USER;
	n_omIntUser = N_OMINT_USER;
	n_OmegaUser = N_OMEGA_USER;
#else
	n_omUser = 0;
	n_omIntUser = 0;
	n_OmegaUser = 0;
#endif
	n_omIntAll = N_OMINT_ALL;
#endif

	init_constants(gdata);

#if(FLUID_TYPE != CRONOS_MULTIFLUID)
	//Trafo = new Transformations(Problem, false);
	Trafo = std::make_unique<Transformations>(gdata.fluid, Problem, false);
#else
	// Make trafo for each fluid
	for(int iFluid=0; iFluid<numFluids; ++iFluid) {
		TrafoMulti.push_back(new Transformations(gdata.fluids->fluids[iFluid], Problem,false,iFluid));
	}
#endif

	IntegrateA = true;
	bcVecPotResized = false;

	debug = static_cast<int>(value((char*)"debug"));

	// string Fluid_Type_string = svalue((char*)"Fluid_Type");
	// if (strcmp(svalue((char*)"Fluid_Type"),"MHD") == 0) {
// #if (FLUID_TYPE == CRONOS_MHD)
// 		Fluid_Type = CRONOS_MHD;
// 	// } else if (strcmp(svalue((char*)"Fluid_Type"),"Hydro") == 0){
// #elif (FLUID_TYPE == CRONOS_HYDRO)
// 		Fluid_Type = CRONOS_HYDRO;
// #endif
	// } else {
	// 	if(gdata.rank == 0) {
	// 		cerr << " No such Fluid_Type available " << endl;
	// 	}
	// 	exit(-23);
	// }

	
	string RiemannSolver = svalue((char*)"RiemannSolver");

	Riemann.resize(DirMax);
	
#if(FLUID_TYPE != CRONOS_MULTIFLUID)

	if(RiemannSolver == "hll") {
		if(gdata.rank==0) {
			cout << " Using the HLL solver " << endl;
		}
		if(FLUID_TYPE == CRONOS_MHD) {
			Riemann[DirX] = std::make_unique<HLLSolver>(gdata, 0, CRONOS_MHD);
			Riemann[DirY] = std::make_unique<HLLSolver>(gdata, 1, CRONOS_MHD);
			Riemann[DirZ] = std::make_unique<HLLSolver>(gdata, 2, CRONOS_MHD);
		} else {
			assert_fail() << "not implemented";
			//RiemannX = std::make_unique<HLLSolver_gen>(gdata, 0, CRONOS_HYDRO);
			//RiemannY = std::make_unique<HLLSolver_gen>(gdata, 1, CRONOS_HYDRO);
			//RiemannZ = std::make_unique<HLLSolver_gen>(gdata, 2, CRONOS_HYDRO);
		}
	} else if(RiemannSolver == "hlld") {
		if(FLUID_TYPE != CRONOS_MHD) {
			if(gdata.rank==0) {
				cerr << " hlld Only available for Fluid_Solver = MHD " << endl;
				exit(2);
			}
		}
		if(gdata.rank==0) {
			cout << " Using the HLLD solver " << endl;
		}
		assert_fail() << "not implemented";
		//RiemannX = std::make_unique<HLLDSolver>(gdata, 0, CRONOS_MHD);
		//RiemannY = std::make_unique<HLLDSolver>(gdata, 1, CRONOS_MHD);
		//RiemannZ = std::make_unique<HLLDSolver>(gdata, 2, CRONOS_MHD);
	} else if(RiemannSolver == "hllc") {
		if(FLUID_TYPE != CRONOS_HYDRO) {
			if(gdata.rank==0) {
				cerr << " hllc Only available for Fluid_Solver = Hydro " << endl;
				exit(2);
			}
		}
#if (ENERGETICS == ISOTHERMAL)
		if(gdata.rank==0) {
			cout << " Using the HLL solver for isothermal case " << endl;
		}
		RiemannX = std::make_unique<HLLSolver_gen>(gdata, 0, CRONOS_HYDRO);
		RiemannY = std::make_unique<HLLSolver_gen>(gdata, 1, CRONOS_HYDRO);
		RiemannZ = std::make_unique<HLLSolver_gen>(gdata, 2, CRONOS_HYDRO);
#else
		if(gdata.rank==0) {
			cout << " Using the HLLC solver " << endl;
		}
		Riemann[DirX] = std::make_unique<HLLCSolver_Hydro>(gdata, 0, CRONOS_HYDRO);
		Riemann[DirY] = std::make_unique<HLLCSolver_Hydro>(gdata, 1, CRONOS_HYDRO);
		Riemann[DirZ] = std::make_unique<HLLCSolver_Hydro>(gdata, 2, CRONOS_HYDRO);
#endif
    } else {
		if(gdata.rank==0) {
			cerr << " No such Riemann solver " << endl;
		}
		cout << " Choices are: " << endl;
		cout << " (hll, hllc, hlld) " << endl;
		exit(-201);
	}

//	// Set indices for Riemann solver
//	for(int iFluid=0; iFluid<numFluids; ++iFluid) {
//		RiemannX->reset_Indices(gdata.fluid);
//		RiemannY->reset_Indices(gdata.fluid);
//		RiemannZ->reset_Indices(gdata.fluid);
//	}
	fieldsX = std::make_unique<fields_1D>(gdata, 0, gdata.fluid);
	fieldsY = std::make_unique<fields_1D>(gdata, 1, gdata.fluid);
	fieldsZ = std::make_unique<fields_1D>(gdata, 2, gdata.fluid);

	//for (int i = 0; i < (int)Direction::Num; ++i) {
	//	fields[i] = std::make_unique<fields_1D>(gdata, i, gdata.fluid);
	//}


	physValxL = std::make_unique<phys_fields_1D>(gdata, 0);
	physValxR = std::make_unique<phys_fields_1D>(gdata, 0);
	
	physValyL = std::make_unique<phys_fields_1D>(gdata, 1);
	physValyR = std::make_unique<phys_fields_1D>(gdata, 1);

	physValzL = std::make_unique<phys_fields_1D>(gdata, 2);
	physValzR = std::make_unique<phys_fields_1D>(gdata, 2);

	ReconstX.resize(TIME_SUBSTEPS);// = new Reconstruction * [TIME_SUBSTEPS];
	ReconstY.resize(TIME_SUBSTEPS);// = new Reconstruction * [TIME_SUBSTEPS];
	ReconstZ.resize(TIME_SUBSTEPS);// = new Reconstruction * [TIME_SUBSTEPS];

	for (int i=0; i < TIME_SUBSTEPS; i++) {
		ReconstX[i] = std::make_unique<Reconstruction>(gdata, 0, gdata.fluid, i);
		ReconstY[i] = std::make_unique<Reconstruction>(gdata, 1, gdata.fluid, i);
		ReconstZ[i] = std::make_unique<Reconstruction>(gdata, 2, gdata.fluid, i);
	}


#if (FLUID_TYPE == CRONOS_MHD)
	PhysFlux = std::make_unique<PhysFluxesMHD>(gdata, gdata.fluid);
#elif (FLUID_TYPE == CRONOS_HYDRO)
	PhysFlux = std::make_unique<PhysFluxesHD>(gdata, gdata.fluid);
#endif

#else // Definitions for multifluid type
	// For multifluid we destinguish between hll and hllX Riemann solver

	// Obtain number of fluids
	numFluids = gdata.fluids->get_numFluids();

	// Make a bunch of Riemann solvers:
	if(RiemannSolver == "hll") {
		if(gdata.rank==0) {
			cout << " Using the HLL solver for all fluids " << endl;
		}

		for(int iFluid=0; iFluid<numFluids; ++iFluid) {
			cout << " fluid: " << iFluid << " " << numFluids << endl;
			// Now prepare each Riemann solvers according to type
			if(gdata.fluids->get_fluidType(iFluid) == CRONOS_MHD) {
				cout << " mag " << endl;
				RiemannSolversX.push_back(new HLLSolver(gdata, 0, CRONOS_MHD));
				RiemannSolversY.push_back(new HLLSolver(gdata, 1, CRONOS_MHD));
				RiemannSolversZ.push_back(new HLLSolver(gdata, 2, CRONOS_MHD));
			} else {
				cout << " hyd " << iFluid << endl;
				RiemannSolversX.push_back(new HLLSolver_gen(gdata, 0, CRONOS_HYDRO));
				RiemannSolversY.push_back(new HLLSolver_gen(gdata, 1, CRONOS_HYDRO));
				RiemannSolversZ.push_back(new HLLSolver_gen(gdata, 2, CRONOS_HYDRO));
			}
		}
	} else if(RiemannSolver == "hllX") {
		if(gdata.rank==0) {
			cout << " Using the improved HLL solver where possible " << endl;
		}
		for(int iFluid=0; iFluid<numFluids; ++iFluid) {
			// Now prepare each Riemann solvers according to type
			if(gdata.fluids->get_fluidType(iFluid) == CRONOS_MHD) {
				RiemannSolversX.push_back(new HLLDSolver(gdata, 0, CRONOS_MHD));
				RiemannSolversY.push_back(new HLLDSolver(gdata, 1, CRONOS_MHD));
				RiemannSolversZ.push_back(new HLLDSolver(gdata, 2, CRONOS_MHD));
			} else {
#if (ENERGETICS == ISOTHERMAL)
//				RiemannSolversX[iFluid] = HLLSolver_gen(gdata, 0, CRONOS_HYDRO);
//				RiemannSolversX[iFluid] = HLLSolver_gen(gdata, 1, CRONOS_HYDRO);
//				RiemannSolversX[iFluid] = HLLSolver_gen(gdata, 2, CRONOS_HYDRO);
				RiemannSolversX.push_back(new HLLSolver_gen(gdata, 0, CRONOS_HYDRO));
				RiemannSolversY.push_back(new HLLSolver_gen(gdata, 1, CRONOS_HYDRO));
				RiemannSolversZ.push_back(new HLLSolver_gen(gdata, 2, CRONOS_HYDRO));
#else
				RiemannSolversX.push_back(new HLLCSolver_Hydro(gdata, 0, CRONOS_HYDRO));
				RiemannSolversY.push_back(new HLLCSolver_Hydro(gdata, 1, CRONOS_HYDRO));
				RiemannSolversZ.push_back(new HLLCSolver_Hydro(gdata, 2, CRONOS_HYDRO));
#endif
			}
		}
    } else {
		if(gdata.rank==0) {
			cerr << " No such Riemann solver " << endl;
		}
		cout << " Choices are: " << endl;
		cout << " (hll, hllX) " << endl;
		exit(-201);
	}
	// Set indices for Riemann solver
	for(int iFluid=0; iFluid<numFluids; ++iFluid) {
		RiemannSolversX[iFluid]->reset_Indices(gdata.fluids->fluids[iFluid]);
		RiemannSolversY[iFluid]->reset_Indices(gdata.fluids->fluids[iFluid]);
		RiemannSolversZ[iFluid]->reset_Indices(gdata.fluids->fluids[iFluid]);
	}

	// Make a bunch of things:
//	fieldsMultiX = new fields_1D[numFluids];
//	fieldsMultiY = new fields_1D[numFluids];
//	fieldsMultiZ = new fields_1D[numFluids];
//	fieldsX = new fields_1D(gdata.fluids, 0);
//	fieldsY = new fields_1D(gdata.fluids, 0);
//	fieldsZ = new fields_1D(gdata.fluids, 0);
//	fieldsX = new fields_1D(gdata, 0);
//	fieldsY = new fields_1D(gdata, 0);
//	fieldsZ = new fields_1D(gdata, 0);

//	// physical fields
//	physValMultixL = new phys_fields_1D[numFluids];
//	physValMultiyL = new phys_fields_1D[numFluids];
//	physValMultizL = new phys_fields_1D[numFluids];
//	physValMultixR = new phys_fields_1D[numFluids];
//	physValMultiyR = new phys_fields_1D[numFluids];
//	physValMultizR = new phys_fields_1D[numFluids];

//	// Reconstruction
//	ReconstMultiX = new Reconstruction[numFluids];
//	ReconstMultiY = new Reconstruction[numFluids];
//	ReconstMultiZ = new Reconstruction[numFluids];

//	// Physical fluxes
//	PhysFluxMulti = new PhysFluxesMHD[numFluids];
	for(int iFluid=0; iFluid<numFluids; ++iFluid) {
		// get number of fields
		int n_omIntLocal = gdata.fluids->get_N_OMINT(iFluid);
		cout << " mache was " << iFluid << endl;
		fieldsMultiX.push_back(new fields_1D(gdata, 0, gdata.fluids->fluids[iFluid]));
		fieldsMultiY.push_back(new fields_1D(gdata, 1, gdata.fluids->fluids[iFluid]));
		fieldsMultiZ.push_back(new fields_1D(gdata, 2, gdata.fluids->fluids[iFluid]));
//		fieldsMultiX[iFluid] = fields_1D(gdata, 0, gdata.fluids->fluids[iFluid]);
//		fieldsMultiY[iFluid] = fields_1D(gdata, 1, gdata.fluids->fluids[iFluid]);
//		fieldsMultiZ[iFluid] = fields_1D(gdata, 2, gdata.fluids->fluids[iFluid]);

//		cout << endl << endl << " x " << fieldsMultiX[iFluid]->omLoc[0].getHigh(0) << endl;
//		cout << endl << " y " << fieldsMultiY[iFluid]->omLoc[0].getHigh(0) <<endl << endl;
//		cout << endl << " z " << fieldsMultiZ[iFluid]->omLoc[0].getHigh(0) <<endl << endl;

//		physValMultixL[iFluid]  = phys_fields_1D(gdata, 0);
//		physValMultiyL[iFluid]  = phys_fields_1D(gdata, 1);
//		physValMultizL[iFluid]  = phys_fields_1D(gdata, 2);
//		physValMultixR[iFluid]  = phys_fields_1D(gdata, 0);
//		physValMultiyR[iFluid]  = phys_fields_1D(gdata, 1);
//		physValMultizR[iFluid]  = phys_fields_1D(gdata, 2);

		physValMultixL.push_back(new phys_fields_1D(gdata, 0));
		physValMultiyL.push_back(new phys_fields_1D(gdata, 1));
		physValMultizL.push_back(new phys_fields_1D(gdata, 2));
		physValMultixR.push_back(new phys_fields_1D(gdata, 0));
		physValMultiyR.push_back(new phys_fields_1D(gdata, 1));
		physValMultizR.push_back(new phys_fields_1D(gdata, 2));

//		ReconstMultiX[iFluid]  = Reconstruction(gdata, 0);
//		ReconstMultiY[iFluid]  = Reconstruction(gdata, 1);
//		ReconstMultiZ[iFluid]  = Reconstruction(gdata, 2);

		ReconstMultiX.push_back(new Reconstruction*[TIME_SUBSTEPS]);
		ReconstMultiY.push_back(new Reconstruction*[TIME_SUBSTEPS]);
		ReconstMultiZ.push_back(new Reconstruction*[TIME_SUBSTEPS]);

		for (int iStep=0; iStep < TIME_SUBSTEPS; ++iStep) {
			ReconstMultiX[iFluid][iStep] = new Reconstruction(gdata, 0, gdata.fluids->fluids[iFluid], iStep);
			ReconstMultiY[iFluid][iStep] = new Reconstruction(gdata, 1, gdata.fluids->fluids[iFluid], iStep);
			ReconstMultiZ[iFluid][iStep] = new Reconstruction(gdata, 2, gdata.fluids->fluids[iFluid], iStep);
		}

		if(gdata.fluids->get_fluidType(iFluid) == CRONOS_MHD) {
			PhysFluxMulti.push_back(new PhysFluxesMHD(gdata, gdata.fluids->fluids[iFluid]));
		} else {
			PhysFluxMulti.push_back(new PhysFluxesHD(gdata, gdata.fluids->fluids[iFluid]));
		}
	}
#endif // End of Multifluid block

#if (OMS_USER == TRUE)
	UserEquations(gdata);
	set_UserPde(gdata);
#endif

	for(int dir=0; dir < DIM; ++dir) {
		gflux[dir].resize(Index::set(-2),Index::set(gdata.mx[dir]+1));
		gflux[dir].clear();
	}

	// Get starting time
	gettimeofday(&tstep, 0);

	// Initialising classes with NULL.
	eos = nullptr;
	sources = NULL;
}


HyperbolicSolver::~HyperbolicSolver() {

#if(FLUID_TYPE != CRONOS_MULTIFLUID)
	//delete Trafo;
#else
	for(int iFluid=0; iFluid<numFluids; ++iFluid) {
		delete TrafoMulti[iFluid];
	}
#endif

	//if(eos != NULL) {
	//	delete eos;
	//}
	//if(sources != NULL) {
	//	delete sources;
	//}

#if(FLUID_TYPE != CRONOS_MULTIFLUID)
	//if(RiemannX != NULL) {
	//	delete RiemannX;
	//}
	//if(RiemannY != NULL) {
	//	delete RiemannY;
	//}
	//if(RiemannZ != NULL) {
	//	delete RiemannZ;
	//}
	//if(fieldsX != NULL) {
	//	delete fieldsX;
	//}
	//if(fieldsY != NULL) {
	//	delete fieldsY;
	//}
	//if(fieldsZ != NULL) {
	//	delete fieldsZ;
	//}

	//if(physValxL != NULL) {
	//	delete physValxL;
	//}
	//if(physValxR != NULL) {
	//	delete physValxR;
	//}
	//if(physValyL != NULL) {
	//	delete physValyL;
	//}
	//if(physValyR != NULL) {
	//	delete physValyR;
	//}
	//if(physValzL != NULL) {
	//	delete physValzL;
	//}
	//if(physValzR != NULL) {
	//	delete physValzR;
	//}

	//if(ReconstX != NULL) {
	//	delete [] ReconstX;
	//}
	//if(ReconstY != NULL) {
	//	delete [] ReconstY;
	//}
	//if(ReconstZ != NULL) {
	//	delete [] ReconstZ;
	//}

	//if(PhysFlux != NULL) {
	//	delete PhysFlux;
	//}
#else
//	delete [] RiemannSolversX;
//	delete [] RiemannSolversY;
//	delete [] RiemannSolversZ;
	for(int iFluid=0; iFluid<numFluids; ++iFluid) {
		delete RiemannSolversX[iFluid];
		delete RiemannSolversY[iFluid];
		delete RiemannSolversZ[iFluid];
	}

//	delete [] fieldsMultiX;
//	delete [] fieldsMultiY;
//	delete [] fieldsMultiZ;
	for(int iFluid=0; iFluid<numFluids; ++iFluid) {
		delete fieldsMultiX[iFluid];
		delete fieldsMultiY[iFluid];
		delete fieldsMultiZ[iFluid];
	}

	for(int iFluid=0; iFluid<numFluids; ++iFluid) {
		delete physValMultixL[iFluid];
		delete physValMultiyL[iFluid];
		delete physValMultizL[iFluid];
		delete physValMultixR[iFluid];
		delete physValMultiyR[iFluid];
		delete physValMultizR[iFluid];
	}

//	delete [] ReconstMultiX;
//	delete [] ReconstMultiY;
//	delete [] ReconstMultiZ;
	for(int iFluid=0; iFluid<numFluids; ++iFluid) {
		delete [] ReconstMultiX[iFluid];
		delete [] ReconstMultiY[iFluid];
		delete [] ReconstMultiZ[iFluid];
	}

//	delete [] PhysFluxMulti;
	for(int iFluid=0; iFluid<numFluids; ++iFluid) {
		delete PhysFluxMulti[iFluid];
	}

#endif
  
#ifdef SAVEMEM
	if(Save != NULL) {
		delete Save;
	}
#endif

#if (FLUID_TYPE == CRONOS_MHD)
	if(CTSolveX != NULL) {
		delete CTSolveX;
	}

	if(CTSolveY != NULL) {
		delete CTSolveY;
	}

	if(CTSolveZ != NULL) {
		delete CTSolveZ;
	}
#endif

	//if(TimeIntegratorGeneric != NULL) {
		//delete [] TimeIntegratorGeneric;
	//}

#if (OMS_USER == TRUE)
	if(TimeIntegratorUser != NULL) {
		delete [] TimeIntegratorUser;
	}

	if(PhysFluxUser != NULL) delete PhysFluxUser;

	delete fieldsXUser;
	delete fieldsYUser;
	delete fieldsZUser;

	delete physValxLUser;
	delete physValxRUser;
	delete physValyLUser;
	delete physValyRUser;
	delete physValzLUser;
	delete physValzRUser;

	delete [] ReconstXUser;
	delete [] ReconstYUser;
	delete [] ReconstZUser;
	
	delete RiemannXUser;
	delete RiemannYUser;
	delete RiemannZUser;

#endif

}


void HyperbolicSolver::init_constants(Data &gdata)
{
	/*
	  Here all (!) relevant quantities have to be defined.
	*/

	yama = value((char*)"Adiabatic_exponent");

	c2_iso = sqr(value((char*)"Isothermal_Soundspeed")); 
	rho0 = value((char*)"Initial_density");


	denomPres = pow(rho0,1-yama);

	thermal = static_cast<int>(value((char*)"thermal"));

	// Values for phystest:
	rhomin = value((char*)"rhomin");
	Ethmin = value((char*)"pmin");
	vmax = value((char*)"vmax");

}


void HyperbolicSolver::addConstants_toH5(Data &gdata, Hdf5Stream &h5out) {
	//! Add constants to each hdf5 output file

	// Write adiabatic exponent
	h5out.AddGlobalAttr("HyperbolicSolver_Gamma", yama);

	// Write isothermal soundspeed
	h5out.AddGlobalAttr("HyperbolicSolver_Isothermal_Soundspeed", sqrt(c2_iso));

	// Write isothermal soundspeed
	h5out.AddGlobalAttr("HyperbolicSolver_SQR_Soundspeed", c2_iso);

	// Write initial density
	h5out.AddGlobalAttr("HyperbolicSolver_Initial_density", rho0);

	// Write allowed extreme values:
	h5out.AddGlobalAttr("HyperbolicSolver_rhomin",rhomin);
	h5out.AddGlobalAttr("HyperbolicSolver_pmin",Ethmin);
	h5out.AddGlobalAttr("HyperbolicSolver_vmax",vmax);

	// Write chosen Riemann solver:
	// string RiemannSolver = svalue((char*)"RiemannSolver");
	// h5out.AddGlobalAttr("HyperbolicSolver_RiemannSolver",RiemannSolver);

	// Indicate if vector potential is integrated:
	h5out.AddGlobalAttr("HyperbolicSolver_IntegrateA", IntegrateA);
	
}

void HyperbolicSolver::getConstants_fromH5(Data &gdata, Hdf5iStream &h5in) {
	//! Get any necessary constant from hdf5 input file:
	
	// So far we only use this to make sure the correct treatment of
	// the magnetic field is used - can be extended in the future.
	int use_vecpot; // Hack to avoid having to write additional bool
					// input-routine.

	// Need to check if attribute is available:
	if(h5in.doesAttrExist("HyperbolicSolver_IntegrateA")) {
		h5in.ReadGlobalAttr("HyperbolicSolver_IntegrateA", use_vecpot);
		is_OldRestart = false;
	} else {
		is_OldRestart = true;
	}
	IntegrateA = use_vecpot;

}


void HyperbolicSolver::init(Data &gdata, gridFunc &gfunc,
                            ProblemType &Problem)
{

	if(gdata.time < 0.1*gdata.dt) {

		// initial conditions
     
		for (int l = 0; l < n_Omega; l++) {
			gdata.om[l].clear();
		}

		if(gdata.rank==0){
			cout << "calling init (type=" << Problem.get_Type() << ")..." << endl;
			cout << Problem.get_Name() << endl;
			cout << gdata.mx[0] << " " << gdata.mx[1] << " " << gdata.mx[2] << endl;
		}
    
		Problem.init_fields(gdata);

#if(CRSWITCH_DUAL_ENERGY == CRONOS_ON)

#if (FLUID_TYPE == CRONOS_MULTIFLUID)
		for(int iFluid=0; iFluid<numFluids; ++iFluid) {
		// Obtain indices:
		q_Eadd = gdata.fluids->get_q_Eadd(iFluid);
		q_rho  = gdata.fluids->get_q_rho(iFluid);
#endif

		gdata.om[q_Eadd] = gdata.om[q_Eges];
#if(AUX_ENERGY == ENTROPY)
		for (int k = -B; k <= gdata.mx[2]+B; ++k) {
			for (int j = -B; j <= gdata.mx[1]+B; ++j) {
				for (int i = -B; i <= gdata.mx[0]+B; ++i) {
					gdata.om[q_Eadd](i,j,k) *= pow(gdata.om[q_rho](i,j,k), 1.-Problem.gamma);
				}
			}
		}
#endif

#if (FLUID_TYPE == CRONOS_MULTIFLUID)
		} // Close loop from above
#endif

#endif


		if(gdata.rank==0) {
			cout << "all ICs set." << endl;
		}

	}

	// Store values at boundaries for fixed boundary conditions
	gfunc.prep_boundaries(gdata, Problem);

	// if Vector-Potential is given -> Transform to magnetic field
	if (FLUID_TYPE == CRONOS_MHD || with_mag) {
		if(gdata.om[q_Bx].getName() == "A_x" ||
				gdata.om[q_By].getName() == "A_y" ||
				gdata.om[q_Bz].getName() == "A_z") {
			TransPot2Mag(gdata, gfunc, Problem);
			IntegrateA = true;
			if(gdata.rank == 0) {
				cout << " Integrating the vector potential " << endl;
			}
		} else {
			IntegrateA = false;
			if(gdata.rank == 0) {
				cout << " Integrating the magnetic induction " << endl;
			}
		}
	}


	init_general(gdata, gfunc, Problem);
}


void HyperbolicSolver::init_general(Data &gdata, gridFunc &gfunc,
                                    ProblemType &Problem) {

	// General initialisation procedure
	// This routine is executed for normal start AND for a restart
	
#if(FLUID_TYPE==CRONOS_MULTIFLUID)
	int i_magFluid = gdata.fluids->get_i_magFluid();
#else
	int i_magFluid = 0;
#endif

	if(with_mag) {
#if CT_TYPE == CONSISTENT
		CTSolveX = new CTLondrilloDelZanna(gdata, 0, IntegrateA, i_magFluid);
		CTSolveY = new CTLondrilloDelZanna(gdata, 1, IntegrateA, i_magFluid);
		CTSolveZ = new CTLondrilloDelZanna(gdata, 2, IntegrateA, i_magFluid);
#elif CT_TYPE == STONE
		CTSolveX = new CTStone(gdata, 0, IntegrateA, i_magFluid);
		CTSolveY = new CTStone(gdata, 1, IntegrateA, i_magFluid);
		CTSolveZ = new CTStone(gdata, 2, IntegrateA, i_magFluid);
#endif
	}

	// Introducing unique ids to all fields om, om_user - can be used
	// for identification (like, e.g., in bcs)
	gdata.set_fieldIds();

	set_TimeIntegrator(gdata, gfunc);


	// Check if using thermal energy or temperature as base variable:
#if (ENERGETICS == FULL)
	if(gdata.om[q_Eges].getName() == "Etherm") {
		thermal = true;
	} else {
		thermal = false;
	}
#if(FLUID_TYPE==CRONOS_MULTIFLUID)
	for(int iFluid=0; iFluid<numFluids; ++iFluid) {
		TrafoMulti[iFluid]->set_thermal(thermal);
	}
#else
	Trafo->set_thermal(thermal);
#endif
#endif

	gfunc.boundary(gdata, Problem); 

#if (USE_COROTATION == CRONOS_ON)
		// Transform inertial frame for phystest
		Trafo->TransCorotToInert(gdata, gfunc, Problem);
#endif

#if(FLUID_TYPE==CRONOS_MULTIFLUID)
	// Loop over all fluids:
	for(int iFluid=0; iFluid<gdata.fluids->get_numFluids(); ++iFluid) {
		CronosFluid fluid = gdata.fluids->fluids[iFluid];
		q_rho = fluid.get_q_rho_global();
		q_sx = fluid.get_q_sx_global();
		q_sy = fluid.get_q_sy_global();
		q_sz = fluid.get_q_sz_global();
		q_Eges = fluid.get_q_Eges_global();
		q_Eadd= fluid.get_q_Eadd_global();
//		Trafo->reset_Indices(fluid);

		/* Initial status checks */
		phystest(gdata, gfunc, Problem, -1, iFluid);
	}

	// Check wheter we have a magnetised fluid
	bool with_MagField = gdata.fluids->with_magField();
	if(with_MagField) {
		compute_divB(gdata, gfunc, Problem);
	}

#else // single-fluid case

	phystest(gdata, gfunc, Problem, -1);

#if(FLUID_TYPE == CRONOS_MHD)
	compute_divB(gdata, gfunc, Problem);
#endif

#endif // end of single-fluid case

#if (USE_COROTATION == CRONOS_ON)
	// Transform back to co-rotating frame
	Trafo->TransInertToCorot(gdata, gfunc, Problem);
#endif

	// Prepare fields to hold changes for magnetic field
#if (FLUID_TYPE == CRONOS_MHD)
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
#endif

}


void HyperbolicSolver::set_TimeIntegrator(const Data &gdata,
                                          gridFunc &gfunc) {

#if(FLUID_TYPE == CRONOS_MULTIFLUID)
	//TimeIntegratorGeneric = new TimeIntegrator*[n_omIntAll];
	TimeIntegratorGeneric.resize(n_omIntAll);
	for(int iom=0; iom<n_omIntAll; iom++) {
#else
	//TimeIntegratorGeneric = new TimeIntegrator*[n_omInt];
	TimeIntegratorGeneric.resize(n_omInt);
	for(int iom=0; iom<n_omInt; iom++) {
#endif

#if (TIME_INTEGRATOR == RUNGEKUTTA)
	//TimeIntegratorGeneric[iom] = new RKSteps;
	TimeIntegratorGeneric[iom] = std::make_unique<RKSteps>();
#elif (TIME_INTEGRATOR == VANLEER)
	//TimeIntegratorGeneric[iom] = new VanLeerIntegrator;
	TimeIntegratorGeneric[iom] = std::make_unique<VanLeerIntegrator>();
#endif
	}

	// Normal settings for 
	for(int q=0; q<n_omInt; ++q) {

		int qchange(q);

		int ibeg[3] = {0,0,0};
		int iend[3] = {gdata.mx[0], gdata.mx[1], gdata.mx[2]};
		
		if(with_mag) {
			if(q>=q_Bx && q<=q_Bz) {
				if(IntegrateA) {
#if(FLUID_TYPE == CRONOS_MULTIFLUID)
					qchange = n_omIntAll+N_ADD+q-q_Bx;
#else
					qchange = n_omInt+N_ADD+q-q_Bx;
#endif
					if(q == q_Bx) {
						iend[1] = gdata.mx[1]+1;
						iend[2] = gdata.mx[2]+1;
					} else if (q == q_By) {
						iend[0] = gdata.mx[0]+1;
						iend[2] = gdata.mx[2]+1;
					} else if (q == q_Bz) {
						iend[0] = gdata.mx[0]+1;
						iend[1] = gdata.mx[1]+1;
					}
				} else {
					if(q == q_Bx) {
						ibeg[0] = -1;
					} else if (q == q_By) {
						ibeg[1] = -1;
					} else if (q == q_Bz) {
						ibeg[2] = -1;
					}
				}
			}
		}

		if(gfunc.get_bc_Type(1) == 1) {
			iend[0] -= 1;
		}

		TimeIntegratorGeneric[q]->set_corrField(qchange);
		TimeIntegratorGeneric[q]->set_IntRange(ibeg, iend);

	}

		
#if (OMS_USER == TRUE)
	// Setting user-oms just to normalgrid-size - maybe this should be
	// made accessible for the user?

	TimeIntegratorUser = new TimeIntegrator*[n_omIntUser];
	for(int iom=0; iom<n_omIntUser; iom++) {

#if (TIME_INTEGRATOR == RUNGEKUTTA)
		TimeIntegratorUser[iom] = new RKSteps;
#elif (TIME_INTEGRATOR == VANLEER)
		TimeIntegratorUser[iom] = new VanLeerIntegrator;
#endif
	}

	for(int q=0; q<n_omIntUser; ++q) {

		int qchange(q), qold(q);

		int ibeg[3] = {0,0,0};
		int iend[3] = {gdata.mx[0], gdata.mx[1], gdata.mx[2]};

		TimeIntegratorUser[q]->set_corrField(qchange);
		TimeIntegratorUser[q]->set_IntRange(ibeg, iend);

	}
#endif		

	

}



void HyperbolicSolver::gen_MagBoundValues(Data &gdata, gridFunc &gfunc,
                                          ProblemType &Problem) {
	//! Compute the values at the boundary for the magnetic field


	// Now for the special case that we do not have the boundary cells
	// for the magnetic induction we need to include them via the
	// divergence constraint:

	int rim = B;

	Pot & Bx = gdata.om[q_Bx];
	Pot & By = gdata.om[q_By];
	Pot & Bz = gdata.om[q_Bz];


	gfunc.boundary(gdata, Problem, gdata.om[q_Bx], B, q_Bx);
	gfunc.boundary(gdata, Problem, gdata.om[q_By], B, q_By);
	gfunc.boundary(gdata, Problem, gdata.om[q_Bz], B, q_Bz);
//	for (int q = 4; q < N_OMINT; q++) {
//		gfunc.boundary(gdata, Problem, gdata.om[q],B,q);
//	}

#if (NON_LINEAR_GRID == CRONOS_OFF)
	REAL idx[3]= {(1./gdata.dx[0]), (1./gdata.dx[1]), (1./gdata.dx[2])};
	REAL dx[3]= {(gdata.dx[0]), (gdata.dx[1]), (gdata.dx[2])};
#endif

	// Lower x-boundary:


#ifdef parallel
	if(std::abs(gdata.global_xb[0] - gdata.xb[0]) < gdata.dx[0])
#endif
	for (int iz = -rim+1; iz <= gdata.mx[2]+rim; ++iz){
		for (int iy = -rim+1; iy <= gdata.mx[1]+rim; ++iy){

#if (NON_LINEAR_GRID == CRONOS_ON)
			REAL idx[3] = {(gdata.getCen_idx(0,0)), (gdata.getCen_idx(1,iy)),
			               (gdata.getCen_idx(2,iz))};
			REAL dx[3]= {(gdata.getCen_dx(0,0)), (gdata.getCen_dx(1,iy)),
			             (gdata.getCen_dx(2,iz))};
#endif
				
#if (GEOM==1) // Cartesian version
			Bx(-1,iy,iz) = Bx(0,iy,iz) + dx[0]*((By(0,iy  ,iz) -
			                                     By(0,iy-1,iz))*idx[1] +
			                                    (Bz(0,iy,iz  ) -
			                                     Bz(0,iy,iz-1))*idx[2]);
#else // Non-Cartesian version
			double inv_const = 1./(gdata.h1(0,iy,iz,-1,0,0)*gdata.h2(0,iy,iz,-1,0,0));
			Bx(-1,iy,iz) = (gdata.h1(0,iy,iz, 1,0,0)*gdata.h2(0,iy,iz, 1,0,0)*Bx(0,iy,iz) +
			                dx[0]*((gdata.h0(0,iy,iz,0, 1,0)*gdata.h2(0,iy,iz,0, 1,0)*By(0,iy  ,iz) - 
			                        gdata.h0(0,iy,iz,0,-1,0)*gdata.h2(0,iy,iz,0,-1,0)*By(0,iy-1,iz))*idx[1] +
			                       (gdata.h0(0,iy,iz,0,0, 1)*gdata.h1(0,iy,iz,0,0, 1)*Bz(0,iy,iz  ) -
			                        gdata.h0(0,iy,iz,0,0,-1)*gdata.h1(0,iy,iz,0,0,-1)*Bz(0,iy,iz-1))*idx[2]))*inv_const;
			
#endif
		}
	}

	// Lower y-boundary:
#ifdef parallel
	if(std::abs(gdata.global_xb[1] - gdata.xb[1]) < gdata.dx[1])
#endif
	for (int iz = -rim+1; iz <= gdata.mx[2]+rim; ++iz){
		for (int ix = -rim+1; ix <= gdata.mx[0]+rim; ++ix){

#if (NON_LINEAR_GRID == CRONOS_ON)
			REAL idx[3] = {(gdata.getCen_idx(0,ix)), (gdata.getCen_idx(1,0)),
			               (gdata.getCen_idx(2,iz))};
			REAL dx[3]= {(gdata.getCen_dx(0,ix)), (gdata.getCen_dx(1,0)),
			             (gdata.getCen_dx(2,iz))};
#endif
			
#if (GEOM==1) // Cartesian version
			By(ix,-1,iz) = By(ix,0,iz) + dx[1]*((Bx(ix  ,0,iz  ) -
			                                     Bx(ix-1,0,iz  ))*idx[0] +
			                                    (Bz(ix  ,0,iz  ) -
			                                     Bz(ix  ,0,iz-1))*idx[2]);
#else // Non-Cartesian version
			double inv_const = 1./(gdata.h0(ix,0,iz,0,-1,0)*gdata.h2(ix,0,iz,0,-1,0));
			By(ix,-1,iz) = (gdata.h0(ix,0,iz,0, 1,0)*gdata.h2(ix,0,iz,0, 1,0)*By(ix,0,iz) +
			                dx[1]*((gdata.h1(ix,0,iz, 1,0,0)*gdata.h2(ix,0,iz, 1,0,0)*Bx(ix  ,0,iz  ) - 
			                        gdata.h1(ix,0,iz,-1,0,0)*gdata.h2(ix,0,iz,-1,0,0)*Bx(ix-1,0,iz  ))*idx[0] +
			                       (gdata.h0(ix,0,iz,0,0, 1)*gdata.h1(ix,0,iz,0,0, 1)*Bz(ix  ,0,iz  ) -
			                        gdata.h0(ix,0,iz,0,0,-1)*gdata.h1(ix,0,iz,0,0,-1)*Bz(ix  ,0,iz-1))*idx[2]))*inv_const;
#endif

		}
	}

	// Lower z-boundary
#ifdef parallel
	if(std::abs(gdata.global_xb[2] - gdata.xb[2]) < gdata.dx[2])
#endif
	for (int iy = -rim+1; iy <= gdata.mx[1]+rim; ++iy){
		for (int ix = -rim+1; ix <= gdata.mx[0]+rim; ++ix){

#if (NON_LINEAR_GRID == CRONOS_ON)
			REAL idx[3] = {(gdata.getCen_idx(0,ix)), (gdata.getCen_idx(1,iy)),
			               (gdata.getCen_idx(2,0))};
			REAL dx[3]= {(gdata.getCen_dx(0,ix)), (gdata.getCen_dx(1,iy)),
			             (gdata.getCen_dx(2,0))};
#endif

#if (GEOM==1) // Cartesian version
			Bz(ix,iy,-1) = Bz(ix,iy,0) + dx[2]*((Bx(ix  ,iy  ,0) -
			                                     Bx(ix-1,iy  ,0))*idx[0] +
			                                    (By(ix  ,iy  ,0) -
			                                     By(ix  ,iy-1,0))*idx[1]);
#else // Non-Cartesian version
			double inv_const = 1./(gdata.h0(ix,iy,0,0,0,-1)*gdata.h1(ix,iy,0,0,0,-1));
			Bz(ix,iy,-1) = (gdata.h0(ix,iy,0,0,0, 1)*gdata.h1(ix,iy,0,0,0, 1)*Bz(ix,iy, 0) +
			                dx[2]*((gdata.h1(ix,iy,0, 1,0,0)*gdata.h2(ix,iy,0, 1,0,0)*Bx(ix  ,iy  ,0) - 
			                        gdata.h1(ix,iy,0,-1,0,0)*gdata.h2(ix,iy,0,-1,0,0)*Bx(ix-1,iy  ,0))*idx[0] +
			                       (gdata.h0(ix,iy,0,0, 1,0)*gdata.h2(ix,iy,0,0, 1,0)*By(ix  ,iy  ,0) - 
			                        gdata.h0(ix,iy,0,0,-1,0)*gdata.h2(ix,iy,0,0,-1,0)*By(ix  ,iy-1,0))*idx[1]))*inv_const;
#endif
		}
	}

#if (NON_LINEAR_GRID == CRONOS_ON)
	cerr << " div B recovery not implemented yet " << endl;
	exit(3);
#endif
}



void HyperbolicSolver::restart(Data &gdata, gridFunc &gfunc,
                               ProblemType &Problem) {
	// If a vector potential is supplied in the restart process -
	// compute the magnetic induction
//#if (FLUID_TYPE == CRONOS_MHD)
	if(with_mag) {
		// Only needed for old restart files -- otherwise the option should
		// be fixed already.

		if(is_OldRestart) {
			if(N_SUBS > 2 && gdata.om[n_omInt+N_ADD].getName() == "A_x" &&
					gdata.om[n_omInt+N_ADD+1].getName() == "A_y" &&
					gdata.om[n_omInt+N_ADD+2].getName() == "A_z") {
				IntegrateA = true;
				compute_B(gdata, gfunc, Problem, true);
				if(gdata.rank == 0) {
					cout << " Integrating the vector potential " << endl;
				}
			} else {
				if(gdata.rank == 0) {
					cout << " Integrating the magnetic induction " << endl;
				}
			}
		} else {
			if(IntegrateA) {
				if(gdata.rank == 0) {
					cout << " Integrating the vector potential " << endl;
				}
				// Better recompute magnetic field - to avoid problems in case of grid resize
				compute_B(gdata, gfunc, Problem, true);
			} else {
				if(gdata.rank == 0) {
					cout << " Integrating the magnetic induction " << endl;
				}
			}
		}

		if(BOUT_DBL < 1) {
			gen_MagBoundValues(gdata, gfunc, Problem);
			gfunc.prep_boundaries(gdata, Problem);
		}
	} else {
		// Store boundary values for fixed boundary conditions
		gfunc.prep_boundaries(gdata, Problem);
	}
	init_general(gdata, gfunc, Problem);
}



REAL HyperbolicSolver::compute_divB(Data &gdata, gridFunc &gfunc,
                                    ProblemType &Problem)
{
	if(Problem.get_Info() && Problem.checkout(4)) {

		REAL divVal;
		if(N_ADD >= 2) {
#if(FLUID_TYPE == CRONOS_MULTIFLUID)
			Pot & divB = gdata.om[n_omIntAll+1];
#else
			Pot & divB = gdata.om[n_omInt+1];
#endif

			divVal = compute_divB(gdata, gfunc, Problem, divB);
		} else {
			NumMatrix<double,3> divB(Index::set(-3,-3,-3),
			                         Index::set(gdata.mx[0]+3,gdata.mx[1]+3,gdata.mx[2]+3));
			divVal = compute_divB(gdata, gfunc, Problem, divB);
		}
		return divVal;
	} else {
		return 42.;
	}
}



REAL HyperbolicSolver::compute_divB(Data &gdata, gridFunc &gfunc,
                                    ProblemType &Problem,
                                    NumMatrix<double,3> &divB)
{
	// Calculation of divB and saving into additional array

	double divBSum  = 0;
	double BSum     = 0;
	double JSum     = 0;
	double max_divB = -10;
	double min_divB = 10;
	int minloc[3] = {0, 0, 0};
	int maxloc[3] = {0, 0, 0};
	int min[3] = {0, 0, 0};
	int max[3] = {gdata.mx[0], gdata.mx[1], gdata.mx[2]};
  
	if(gfunc.get_bc_Type(1) == 1) {
		max[0] = gdata.mx[0]-1;
	}
	if(gfunc.get_bc_Type(3) == 1) {
		max[1] = gdata.mx[1]-1;
	}
	if(gfunc.get_bc_Type(5) == 1) {
		max[2] = gdata.mx[2]-1;
	}
  
#if (NON_LINEAR_GRID == CRONOS_OFF)
	REAL idx[3]= {(1./gdata.dx[0]), (1./gdata.dx[1]), (1./gdata.dx[2])};
#endif
	for (int k = min[2]; k <= max[2]; k++) {   // compute J = rot(B) and 
		for (int j = min[1]; j <= max[1]; j++) {   // div(B) as om fields
			for (int i = min[0]; i <= max[0]; i++) {

#if (NON_LINEAR_GRID == CRONOS_ON)
				REAL idx[3] = {(gdata.getCen_idx(0,i)), (gdata.getCen_idx(1,j)),
				               (gdata.getCen_idx(2,k))};
#endif


#if GEOM > 1
				REAL f_geom = 1./gdata.get_CellGeomTrafo(i,j,k);
#endif
#if GEOM > 1
				double dxBx = (gdata.h1(i,j,k, 1,0,0)*gdata.h2(i,j,k, 1,0,0)*gdata.om[q_Bx](i,j,k) -
				               gdata.h1(i,j,k,-1,0,0)*gdata.h2(i,j,k,-1,0,0)*gdata.om[q_Bx](i-1,j,k))*idx[0]*f_geom;
				double dxBy = ((gdata.h1(i+1,j,k,0, 1,0)*gdata.om[q_By](i+1,j  ,k) -
				                gdata.h1(i-1,j,k,0, 1,0)*gdata.om[q_By](i-1,j  ,k)) +
				               (gdata.h1(i+1,j,k,0,-1,0)*gdata.om[q_By](i+1,j-1,k) -
				                gdata.h1(i-1,j,k,0,-1,0)*gdata.om[q_By](i-1,j-1,k)))*idx[0]*0.25;
				double dxBz = ((gdata.h2(i+1,j,k,0,0, 1)*gdata.om[q_Bz](i+1,j,k  ) -
				                gdata.h2(i-1,j,k,0,0, 1)*gdata.om[q_Bz](i-1,j,k  )) +
				               (gdata.h2(i+1,j,k,0,0,-1)*gdata.om[q_Bz](i+1,j,k-1) -
				                gdata.h2(i-1,j,k,0,0,-1)*gdata.om[q_Bz](i-1,j,k-1)))*idx[0]*0.25;
#else
				double dxBx = (gdata.om[q_Bx](i  ,j,k)-gdata.om[q_Bx](i-1,j,k))*idx[0];
				double dxBy = ((gdata.om[q_By](i+1,j  ,k) -
				                gdata.om[q_By](i-1,j  ,k)) +
				               (gdata.om[q_By](i+1,j-1,k) -
				                gdata.om[q_By](i-1,j-1,k)))*idx[0]*0.25;
				double dxBz = ((gdata.om[q_Bz](i+1,j,k  ) -
				                gdata.om[q_Bz](i-1,j,k  )) +
				               (gdata.om[q_Bz](i+1,j,k-1) - 
				                gdata.om[q_Bz](i-1,j,k-1)))*idx[0]*0.25;
#endif
	
#if (GEOM > 1)
				double dyBx = ((gdata.h0(i,j+1,k, 1,0,0)*gdata.om[q_Bx](i  ,j+1,k) - 
				                gdata.h0(i,j-1,k, 1,0,0)*gdata.om[q_Bx](i  ,j-1,k)) +
				               (gdata.h0(i,j+1,k,-1,0,0)*gdata.om[q_Bx](i-1,j+1,k) -
				                gdata.h0(i,j-1,k,-1,0,0)*gdata.om[q_Bx](i-1,j-1,k)))*idx[1]*0.25;
				double dyBy = (gdata.h0(i,j,k,0, 1,0)*gdata.h2(i,j,k,0, 1,0)*gdata.om[q_By](i,j,k) - 
				               gdata.h0(i,j,k,0,-1,0)*gdata.h2(i,j,k,0,-1,0)*gdata.om[q_By](i,j-1,k))*idx[1]*f_geom;
				double dyBz = ((gdata.h2(i,j+1,k,0,0, 1)*gdata.om[q_Bz](i,j+1,k  ) - 
				                gdata.h2(i,j-1,k,0,0, 1)*gdata.om[q_Bz](i,j-1,k  )) +
				               (gdata.h2(i,j+1,k,0,0,-1)*gdata.om[q_Bz](i,j+1,k-1) -
				                gdata.h2(i,j-1,k,0,0,-1)*gdata.om[q_Bz](i,j-1,k-1)))*idx[1]*0.25;
#else
				double dyBx = ((gdata.om[q_Bx](i  ,j+1,k) -
				                gdata.om[q_Bx](i  ,j-1,k)) +
				               (gdata.om[q_Bx](i-1,j+1,k) -
				                gdata.om[q_Bx](i-1,j-1,k)))*idx[1]*0.25;
				double dyBy = (gdata.om[q_By](i,j,k) - gdata.om[q_By](i,j-1,k))*idx[1];
				double dyBz = ((gdata.om[q_Bz](i,j+1,k  ) - 
				                gdata.om[q_Bz](i,j-1,k  )) +
				               (gdata.om[q_Bz](i,j+1,k-1) - 
				                gdata.om[q_Bz](i,j-1,k-1)))*idx[1]*0.25;
#endif
	
#if (GEOM > 1)
				double dzBx = ((gdata.h0(i,j,k+1, 1,0,0)*gdata.om[q_Bx](i  ,j,k+1) - 
				                gdata.h0(i,j,k-1, 1,0,0)*gdata.om[q_Bx](i  ,j,k-1)) +
				               (gdata.h0(i,j,k+1,-1,0,0)*gdata.om[q_Bx](i-1,j,k+1) -
				                gdata.h0(i,j,k-1,-1,0,0)*gdata.om[q_Bx](i-1,j,k-1)))*idx[2]*0.25;
				double dzBy = ((gdata.h1(i,j,k+1,0, 1,0)*gdata.om[q_By](i,j  ,k+1) - 
				                gdata.h1(i,j,k-1,0, 1,0)*gdata.om[q_By](i,j  ,k-1)) +
				               (gdata.h1(i,j,k+1,0,-1,0)*gdata.om[q_By](i,j-1,k+1) -
				                gdata.h1(i,j,k-1,0,-1,0)*gdata.om[q_By](i,j-1,k-1)))*idx[2]*0.25;
				double dzBz = (gdata.h0(i,j,k,0,0, 1)*gdata.h1(i,j,k,0,0, 1)*gdata.om[q_Bz](i,j,k) -
				               gdata.h0(i,j,k,0,0,-1)*gdata.h1(i,j,k,0,0,-1)*gdata.om[q_Bz](i,j,k-1))*idx[2]*f_geom;
#else
				double dzBx = ((gdata.om[q_Bx](i  ,j,k+1) - 
				                gdata.om[q_Bx](i  ,j,k-1)) +
				               (gdata.om[q_Bx](i-1,j,k+1) - 
				                gdata.om[q_Bx](i-1,j,k-1)))*idx[2]*0.25;
				double dzBy = ((gdata.om[q_By](i,j  ,k+1) - 
				                gdata.om[q_By](i,j  ,k-1)) +
				               (gdata.om[q_By](i,j-1,k+1) - 
				                gdata.om[q_By](i,j-1,k-1)))*idx[2]*0.25;

				double dzBz = (gdata.om[q_Bz](i,j,k)-gdata.om[q_Bz](i,j,k-1))*idx[2];	
#endif
	
				double BAbs = sqrt(sqr(0.5*(gdata.om[q_Bx](i,j,k) + gdata.om[q_Bx](i-1,j,k))) +
				                   sqr(0.5*(gdata.om[q_By](i,j,k) + gdata.om[q_By](i,j-1,k))) +
				                   sqr(0.5*(gdata.om[q_Bz](i,j,k) + gdata.om[q_Bz](i,j,k-1))));
	
#if (GEOM > 1)
				double Jx   = (dyBz-dzBy)/(gdata.getCen_h1(i,j,k)*gdata.getCen_h2(i,j,k));
				double Jy   = (dzBx-dxBz)/(gdata.getCen_h0(i,j,k)*gdata.getCen_h2(i,j,k));  // J = rot(B)
				double Jz   = (dxBy-dyBx)/(gdata.getCen_h0(i,j,k)*gdata.getCen_h1(i,j,k));
#else
				double Jx   = dyBz-dzBy;
				double Jy   = dzBx-dxBz;  // J = rot(B)
				double Jz   = dxBy-dyBx;
#endif
				double JAbs = sqrt(sqr(Jx)+sqr(Jy)+sqr(Jz));
	
				divB(i,j,k) = dxBx+dyBy+dzBz;
	
				double divBval = divB(i,j,k);
	
				if(divBval > max_divB){
					max_divB = std::max(max_divB,divBval);
					maxloc[0] = i;
					maxloc[1] = j;
					maxloc[2] = k;
				}
				if(divBval < min_divB){
					minloc[0] = i;
					minloc[1] = j;
					minloc[2] = k;
				}
				min_divB = std::min(min_divB,divBval);
				divBSum += abs(divB(i,j,k));
				BSum    += BAbs;
				JSum    += JAbs;
			}
		}
	}

	if(gdata.rank==0) {
		cout << "-------div B:-----------------------------------------" << endl;
		//    cout << "........................div B:................." << endl;
		cout << " Min / Max: " <<  min_divB << " " << max_divB << endl;
		cout << " Location: " << minloc[0] << " " << minloc[1];
		cout << " " << minloc[2] << "    ";
		cout << maxloc[0] << " " << maxloc[1] << " " << maxloc[2] << endl;
	}
	// cout << gdata.rank << " Min / Max: " <<  min_divB << " " << max_divB << " Location: " << minloc[0] << " " << minloc[1] << " " << minloc[2] << endl;
	//   cout << "-------div B:-----------------------------------------" << endl;
	//     //    cout << "........................div B:................." << endl;
	//   cout << " Min / Max: " <<  min_divB << " " << max_divB << endl;
	//   cout << " Location: " << minloc[0] << " " << minloc[1];
	//   cout << " " << minloc[2] << "    ";
	//   cout << maxloc[0] << " " << maxloc[1] << " " << maxloc[2] << " ";
	//   cout << gdata.rank << endl;
#ifdef parallel  
	MPI_Barrier(gdata.comm3d);
#endif

	gfunc.boundary(gdata, Problem, divB,3);
  
#ifdef parallel
	gdata.MpiSum(divBSum);
	gdata.MpiSum(BSum);
	gdata.MpiSum(JSum);
#endif
	if(gdata.rank==0){
		divBSum /= BSum;
		cout << " Ave: " << divBSum << " " << divBSum*BSum/JSum << endl;
		//   cout << "Sum of B: " << BSum << endl;
		//     cout << "---------------------------------------------" << endl;
	}
#ifdef parallel
	MPI_Barrier(gdata.comm3d);
#endif
	REAL divVal =  divBSum;
	return divVal;
}



void HyperbolicSolver::phystest(Data &gdata, gridFunc &gfunc,
                                ProblemType &Problem)
{
	phystest(gdata, gfunc, Problem, -1);
}



void HyperbolicSolver::phystest(Data &gdata, gridFunc &gfunc,
                                ProblemType &Problem, int rkstep, int iFluid)
{

	bool init(false);
	if(rkstep < 0) {
		init = true;
		rkstep = TIME_SUBSTEPS-1;
	}
  
	REAL dV(gdata.get_CellVolume(0,0,0));
	double Vol(gdata.get_Volume());
	double Mass(0.), Eges(0.);
	double mv[3] = {0., 0., 0.};
	double Ekin(0.), Emag(0.), Etherm(0.);
	double cs(0.);
	double Ekfluc(0.), Ebfluc(0.);
	double Vorticity(0.), Enstrophy(0.);
	double sqrdivv(0.), sqrdissv(0.);
	double v2ave(0.);
	REAL mvrms[3] = {0., 0., 0.};
#if (FLUID_TYPE == CRONOS_MHD)
	REAL BAve[3] = {0., 0., 0.};
	REAL Brms[3] = {0., 0., 0.};
#endif
    

#ifdef parallel
	gdata.MpiSum(Vol);
#endif

	// For multifluid case use only relevant range of fields
#if(FLUID_TYPE==CRONOS_MULTIFLUID)
//	cout << " Indizes: " << q_sx << " " << q_sy << endl;
	int q_min = gdata.fluids->fluids[iFluid].get_IndexGlobal(0);
	int q_Eges = gdata.fluids->fluids[iFluid].get_q_Eges_global();
	n_omInt = gdata.fluids->get_N_OMINT(iFluid);
	int q_rho = gdata.fluids->fluids[iFluid].get_q_rho_global();
	int q_sx = gdata.fluids->fluids[iFluid].get_q_sx_global();
	int q_sy = gdata.fluids->fluids[iFluid].get_q_sy_global();
	int q_sz = gdata.fluids->fluids[iFluid].get_q_sz_global();
	int q_Bx = gdata.fluids->fluids[iFluid].get_q_Bx_global();
	int q_By = gdata.fluids->fluids[iFluid].get_q_By_global();
	int q_Bz = gdata.fluids->fluids[iFluid].get_q_Bz_global();
#else
	int q_min = 0;
#endif

	// Enforce min and max values if desired by the user -- this is
	// done just at the beginning in order to compute volume
	// integrated quantities for the corrected values.
	for(int q=q_min; q<n_omInt; ++q) {
		if(Problem.force_min(q)) {
			gdata.om[q].set_min(Problem.min_Val(q));
		}
		if(Problem.force_max(q)) {
			gdata.om[q].set_max(Problem.max_Val(q));
		}
	}

	// Conservative quantities only to be evaluated, when output is
	// desired -- not anymore...
	//  if(Problem.get_AsciiOut() || Problem.get_Info()) {

	// Use conservative variables:
#if (ENERGETICS == FULL)
	if(gdata.om[q_Eges].getName() == "Etherm" || 
	   gdata.om[q_Eges].getName() == "Temp") {
		// Making certain to use velocity for energy trafo
		if(gdata.om[q_sx].getName() == "s_x" && gdata.om[q_sy].getName() == "s_y" &&
		   gdata.om[q_sz].getName() == "s_z") {
#if(FLUID_TYPE==CRONOS_MULTIFLUID)
			TrafoMulti[iFluid]->TransMomen2Vel(gdata, gfunc, Problem);
#else
			Trafo->TransMomen2Vel(gdata, gfunc, Problem);
#endif
		}
		if(thermal) {
#if(FLUID_TYPE==CRONOS_MULTIFLUID)
			TrafoMulti[iFluid]->TransEth2E(gdata, gfunc, Problem);
#else
			Trafo->TransEth2E(gdata, gfunc, Problem);
#endif
		} else {
#if(FLUID_TYPE==CRONOS_MULTIFLUID)
			TrafoMulti[iFluid]->TransT2E(gdata, gfunc, Problem);
#else
			Trafo->TransT2E(gdata, gfunc, Problem);
#endif
		}

	}
#endif
  
	if(gdata.om[q_sx].getName() == "v_x" && gdata.om[q_sy].getName() == "v_y" &&
	   gdata.om[q_sz].getName() == "v_z") {
#if(FLUID_TYPE==CRONOS_MULTIFLUID)
		TrafoMulti[iFluid]->TransVel2Momen(gdata, gfunc, Problem);
#else
		Trafo->TransVel2Momen(gdata, gfunc, Problem);
#endif
	}

	// Compute conservative quantities:
	Mass = gdata.computeInt(q_rho);

	mv[0] = gdata.computeInt(q_sx);
	mv[1] = gdata.computeInt(q_sy);
	mv[2] = gdata.computeInt(q_sz);
  
#if (FLUID_TYPE == CRONOS_MHD)
	if(Problem.mag) {
		BAve[0] = gdata.computeInt(q_Bx);
		BAve[1] = gdata.computeInt(q_By);
		BAve[2] = gdata.computeInt(q_Bz);
		BAve[0] /= Vol;
		BAve[1] /= Vol;
		BAve[2] /= Vol;
	}
#endif  

	if(ENERGETICS == FULL) {
		Eges = gdata.computeInt(q_Eges);
	}
  
	mvrms[0] = gdata.computeRMS(q_sx);
	mvrms[1] = gdata.computeRMS(q_sy);
	mvrms[2] = gdata.computeRMS(q_sz);
#if (FLUID_TYPE == CRONOS_MHD)
	if(Problem.mag) {
		Brms[0] = gdata.computeRMS(q_Bx);
		Brms[1] = gdata.computeRMS(q_By);
		Brms[2] = gdata.computeRMS(q_Bz);
	}
#endif


	// Making certain to use primitve Variables for the rest of the
	// computations

	if(gdata.om[q_sx].getName() == "s_x" && gdata.om[q_sy].getName() == "s_y" &&
	   gdata.om[q_sz].getName() == "s_z") {
#if(FLUID_TYPE==CRONOS_MULTIFLUID)
		TrafoMulti[iFluid]->TransMomen2Vel(gdata, gfunc, Problem);
#else
		Trafo->TransMomen2Vel(gdata, gfunc, Problem);
#endif
	}
// if(ENERGETICS == FULL) {
	// 	if(thermal) {
	// 		Trafo->TransE2Eth(gdata, gfunc, Problem);
	// 	} else {
	// 		Trafo->TransE2T(gdata, gfunc, Problem);
	// 	}
	// }
	if(ENERGETICS == FULL) {
#if(FLUID_TYPE==CRONOS_MULTIFLUID)
		TrafoMulti[iFluid]->TransE2Eth(gdata, gfunc, Problem);
#else
		Trafo->TransE2Eth(gdata, gfunc, Problem);
#endif
	}
	//    exit(2);


	if(rkstep == TIME_SUBSTEPS-1) {
  
		// Energetics only to be evaluated, when output is desired
		//    if((Problem.get_AsciiOut() || Problem.get_Info()) && !init) {
		// Now energetics is computed every time -- values are needed for averaging

		// Compute Energetics & speed of sound
		Vol = 0.;
		for (int k = 0; k < gdata.mx[2]; k++) {
			for (int j = 0; j < gdata.mx[1]; j++) {
				for (int i = 0; i < gdata.mx[0]; i++) {
#if ((GEOM > 1) || (NON_LINEAR_GRID == CRONOS_ON))
					dV = gdata.get_CellVolume(i,j,k);
#endif
	  
					Ekin += (sqr(gdata.om[q_sx](i,j,k)) + 
					         sqr(gdata.om[q_sy](i,j,k)) +
					         sqr(gdata.om[q_sz](i,j,k)))*0.5*dV*gdata.om[q_rho](i,j,k);
#if (FLUID_TYPE == CRONOS_MHD)
					Emag += (sqr(0.5*(gdata.om[q_Bx](i,j,k) + 
					                  gdata.om[q_Bx](i-1,j,k))) +
					         sqr(0.5*(gdata.om[q_By](i,j,k) + 
					                  gdata.om[q_By](i,j-1,k))) +
					         sqr(0.5*(gdata.om[q_Bz](i,j,k) + 
					                  gdata.om[q_Bz](i,j,k-1))))*0.5*dV;
#endif

					if(ENERGETICS == FULL) {
						Etherm += gdata.om[q_Eges](i,j,k)*dV;
						cs += sqrt((yama - 1.)*gdata.om[q_Eges](i,j,k)/gdata.om[q_rho](i,j,k))*dV;
					} else {
						cs += Problem.c2_iso(gdata,i,j,k)*dV;
					}
	  
					Vol += dV;
				}
			}
		}
    

#ifdef parallel
		gdata.MpiSum(Ekin);
#if (FLUID_TYPE == CRONOS_MHD)
		gdata.MpiSum(Emag);
#endif
		gdata.MpiSum(Etherm);
		gdata.MpiSum(cs);
		gdata.MpiSum(Vol);
#endif
    
		cs = cs/Vol;
		//     if(N_OMINT < 8) {
		//       Etherm = Mass*sqr(cs);
		//       if(yama != 1.) {
		// 	Etherm /= (yama - 1.);
		//       }
		//       Eges = Etherm + Ekin + Emag;
		//     }
		if(ENERGETICS != FULL) {
			Eges = Ekin;
#if (FLUID_TYPE == CRONOS_MHD)
			Eges += Emag;
#endif
		}
    
    
		// Set initial values at first step:
		if(init) {
			Mass0 = Mass;
			mv0[0] = mv[0];
			mv0[1] = mv[1];
			mv0[2] = mv[2];
#if (FLUID_TYPE == CRONOS_MHD)
			BAve0[0] = BAve[0];
			BAve0[1] = BAve[1];
			BAve0[2] = BAve[2];
#endif
			Eges0 = Eges;
		}
    
		// Saving time-averaged data:
		Ekin_tave += gdata.dt*Ekin;
#if (FLUID_TYPE == CRONOS_MHD)
		Emag_tave += gdata.dt*Emag;
#endif
		Etherm_tave += gdata.dt*Etherm;
		Eges_tave += gdata.dt*Eges;
		cs_tave += gdata.dt*cs;
		time_ave += gdata.dt;


		// Compute Velocity Stats:
		Vol = 0.;
#if (NON_LINEAR_GRID == CRONOS_OFF)
		REAL hx[3]= {(gdata.hx[0]), (gdata.hx[1]), (gdata.hx[2])};
#endif
		for (int k = 0; k < gdata.mx[2]; k++) {
			for (int j = 0; j < gdata.mx[1]; j++) {
				for (int i = 0; i < gdata.mx[0]; i++) {

#if ((GEOM > 1) || (NON_LINEAR_GRID == CRONOS_ON))
					dV = gdata.get_CellVolume(i,j,k);
#endif
#if (NON_LINEAR_GRID == CRONOS_ON)
					REAL hx[3]= {(gdata.getCen_hx(0,i)), (gdata.getCen_hx(1,j)),
					             (gdata.getCen_hx(2,k))};	
#endif

	  
#if (GEOM > 1)
					REAL f_geom = 1./gdata.get_CellGeomTrafo(i,j,k);

					REAL dxVx = (gdata.getCen_h1(i+1,j,k)*gdata.getCen_h2(i+1,j,k)*gdata.om[q_sx](i+1,j,k) -
					             gdata.getCen_h1(i-1,j,k)*gdata.getCen_h2(i-1,j,k)*gdata.om[q_sx](i-1,j,k))*hx[0]*f_geom;
					REAL dxVy = (gdata.getCen_h1(i+1,j,k)*gdata.om[q_sy](i+1,j,k) -
					             gdata.getCen_h1(i-1,j,k)*gdata.om[q_sy](i-1,j,k))*hx[0];
					REAL dxVz = (gdata.getCen_h2(i+1,j,k)*gdata.om[q_sz](i+1,j,k) -
					             gdata.getCen_h2(i-1,j,k)*gdata.om[q_sz](i-1,j,k))*hx[0];
					REAL dyVx = (gdata.getCen_h0(i,j+1,k)*gdata.om[q_sx](i,j+1,k) -
					             gdata.getCen_h0(i,j-1,k)*gdata.om[q_sx](i,j-1,k))*hx[1];
					REAL  dyVy = (gdata.getCen_h0(i,j+1,k)*gdata.getCen_h2(i,j+1,k)*gdata.om[q_sy](i,j+1,k) -
					              gdata.getCen_h0(i,j-1,k)*gdata.getCen_h2(i,j-1,k)*gdata.om[q_sy](i,j-1,k))*hx[1]*f_geom;
					REAL dyVz = (gdata.getCen_h2(i,j+1,k)*gdata.om[q_sz](i,j+1,k) -
					             gdata.getCen_h2(i,j-1,k)*gdata.om[q_sz](i,j-1,k))*hx[1];
					REAL dzVx = (gdata.getCen_h0(i,j,k+1)*gdata.om[q_sx](i,j,k+1) -
					             gdata.getCen_h0(i,j,k-1)*gdata.om[q_sx](i,j,k-1))*hx[2];
					REAL dzVy = (gdata.getCen_h1(i,j,k+1)*gdata.om[q_sy](i,j,k+1) -
					             gdata.getCen_h1(i,j,k-1)*gdata.om[q_sy](i,j,k-1))*hx[2];
					REAL dzVz = (gdata.getCen_h0(i,j,k+1)*gdata.getCen_h1(i,j,k+1)*gdata.om[q_sz](i,j,k+1) -
					             gdata.getCen_h0(i,j,k-1)*gdata.getCen_h1(i,j,k-1)*gdata.om[q_sz](i,j,k-1))*hx[2]*f_geom;
#else
					double dxVx = (gdata.om[q_sx](i+1,j,k) -
					               gdata.om[q_sx](i-1,j,k))*hx[0];
					double dxVy = (gdata.om[q_sy](i+1,j,k) -
					               gdata.om[q_sy](i-1,j,k))*hx[0];
					double dxVz = (gdata.om[q_sz](i+1,j,k) -
					               gdata.om[q_sz](i-1,j,k))*hx[0];
					double dyVx = (gdata.om[q_sx](i,j+1,k) -
					               gdata.om[q_sx](i,j-1,k))*hx[1];
					double dyVy = (gdata.om[q_sy](i,j+1,k) -
					               gdata.om[q_sy](i,j-1,k))*hx[1];
					double dyVz = (gdata.om[q_sz](i,j+1,k) -
					               gdata.om[q_sz](i,j-1,k))*hx[1];
					double dzVx = (gdata.om[q_sx](i,j,k+1) -
					               gdata.om[q_sx](i,j,k-1))*hx[2];
					double dzVy = (gdata.om[q_sy](i,j,k+1) -
					               gdata.om[q_sy](i,j,k-1))*hx[2];
					double dzVz = (gdata.om[q_sz](i,j,k+1) -
					               gdata.om[q_sz](i,j,k-1))*hx[2];
#endif
	  
#if (GEOM > 1)
					REAL curl_x = (dyVz - dzVy)/(gdata.getCen_h1(i,j,k)*gdata.getCen_h2(i,j,k));
					REAL curl_y = (dzVx - dxVz)/(gdata.getCen_h0(i,j,k)*gdata.getCen_h2(i,j,k));
					REAL curl_z = (dxVy - dyVx)/(gdata.getCen_h0(i,j,k)*gdata.getCen_h1(i,j,k));
#else
					REAL curl_x = (dyVz - dzVy);
					REAL curl_y = (dzVx - dxVz);
					REAL curl_z = (dxVy - dyVx);
#endif
					REAL curl_val = (sqr(curl_x) + sqr(curl_y) + sqr(curl_z));
	    
					Vorticity += sqrt(curl_val)*dV;
					Enstrophy += (curl_val)*dV;
	  
					sqrdivv  += sqr(dxVx + dyVy + dzVz);
					sqrdissv += (sqr(dxVx) + sqr(dyVx) + sqr(dzVz) +
					             sqr(dxVy) + sqr(dyVy) + sqr(dzVy) +
					             sqr(dxVz) + sqr(dyVz) + sqr(dzVz));
	  
					v2ave += (sqr(gdata.om[q_sx](i,j,k)) +
					          sqr(gdata.om[q_sy](i,j,k)) +
					          sqr(gdata.om[q_sz](i,j,k)))*dV;
	  
					Vol += dV;
	  
				}
			}
		}
#ifdef parallel
		gdata.MpiSum(Vorticity);
		gdata.MpiSum(Enstrophy);
		gdata.MpiSum(sqrdivv);
		gdata.MpiSum(sqrdissv);
		gdata.MpiSum(Vol);
#endif
		v2ave /= Vol;
    
		Problem.computeFluct(gdata, Ekfluc, Ebfluc);

#ifdef parallel
		gdata.MpiSum(Ekfluc);
#if (FLUID_TYPE == CRONOS_MHD)
		gdata.MpiSum(Ebfluc);
#endif
		//     Ebfluc /= gdata.nproc[0]*gdata.nproc[1]*gdata.nproc[2];
		//     Ekfluc /= gdata.nproc[0]*gdata.nproc[1]*gdata.nproc[2];
#endif  

    
		// Compute temporal averages (in case no output at every timestep
		Vorticity_tave += gdata.dt*Vorticity;
		Enstrophy_tave += gdata.dt*Enstrophy;
		v2ave_tave += gdata.dt*v2ave;
		Ekfluc_tave += gdata.dt*Ekfluc;
#if (FLUID_TYPE == CRONOS_MHD)
		Ebfluc_tave += gdata.dt*Ebfluc;
#endif
		// Compute quantities in derived class (if necessary)
		Problem.computePhystest(gdata);
    
		// Writing to screen -- if desired
		if(Problem.get_Info() && Problem.checkout(3)) {
			double MachAve = sqrt(v2ave)/cs;
    
			if(gdata.rank == 0) {
#if(FLUID_TYPE == CRONOS_MULTIFLUID)
				cout << " Fluid: " << iFluid << endl;
#endif
				cout << setiosflags(ios::scientific) << setprecision(3) << setw(12);
				cout << "-------Fluctuations-----------------------------------" << endl;
				cout << " E_kin:  " << setw(9) << Ekfluc;
#if (FLUID_TYPE == CRONOS_MHD)
				if(Problem.mag) {
					cout << " E_mag:  " << setw(9) << Ebfluc;
				}
#endif
				cout << " E_dyn:  " << setw(9) << Ekfluc + Ebfluc << endl;
				cout << "-------Energertics:-----------------------------------" << endl;
				cout << " E_ges:  " << setw(9) << Eges;
				if(ENERGETICS == FULL){
					cout << " E_th:   " << setw(9) << Etherm;
				}
	
				if(FLUID_TYPE == CRONOS_MHD) {
					cout << " E_mag:  " << setw(9) << Emag;
				}
				cout << endl;
				cout << "-------Dynamics:--------------------------------------" << endl;
				cout << " Vorti.: " << setw(9) << Vorticity;
				cout << " Enstr.: " << setw(9) << Enstrophy;
				cout << " M_rms:  " << setw(9) << MachAve << endl;
				cout << "-------Conserved quantities:--------------------------" << endl;
				cout << " Mass: " << setw(9) << Mass << " err: ";
				cout << setw(9) << Mass-Mass0 << endl;
				cout << " p_x:  " << setw(9) << mv[0] << " err: ";
				cout << setw(9)<< (mv[0]-mv0[0])/(mvrms[0]+1.e-20) << endl;
				cout << " p_y:  " << setw(9) << mv[1] << " err: ";
				cout << setw(9)<< (mv[1]-mv0[1])/(mvrms[1]+1.e-20) << endl;
				cout << " p_z:  " << setw(9) << mv[2] << " err: ";
				cout << setw(9) << (mv[2]-mv0[2])/(mvrms[2]+1.e-20) << endl;
#if (FLUID_TYPE == CRONOS_MHD)
				cout << " Bx:   " << setw(9) << BAve[0] << " err: ";
				cout << setw(9) << (BAve[0]-BAve0[0])/(Brms[0]+1.e-20) << endl;
				cout << " By:   " << setw(9) << BAve[1] << " err: ";
				cout << setw(9) << (BAve[1]-BAve0[1])/(Brms[1]+1.e-20)<< endl;
				cout << " Bz:   " << setw(9) << BAve[2] << " err: ";
				cout << setw(9) << (BAve[2]-BAve0[2])/(Brms[2]+1.e-20)<< endl;
#endif
				if(ENERGETICS == FULL) {
					cout << " Eges: " << setw(9) << Eges << " err: ";
					cout << setw(9) << Eges-Eges0 << endl;
				}
				// 	if(gdata.time == 0.) {
				// 	  cout << "------------------------------------------------------";
				// 	  cout << endl;
				// 	}
				cout << resetiosflags(ios::scientific) << setprecision(6);
			}
			Problem.writePhystestInfo(gdata);
		}
    


		// Writing to output file -- if desired
		// Only here we use the temporal average:


		if(Problem.get_AsciiOut()) {
			Problem.set_AsciiOut(false);

			Ekfluc = Ekfluc_tave/time_ave;
#if (FLUID_TYPE == CRONOS_MHD)
			Ebfluc = Ebfluc_tave/time_ave;
#endif
			Vorticity = Vorticity_tave/time_ave;
			Etherm = Etherm_tave/time_ave;

			char fileenerg[255];
			sprintf(fileenerg,"%s/%s_efluct",getenv("poub"),getenv("pname"));
  
			ofstream outfile;
			if(gdata.rank == 0) {
				if(gdata.time == 0. && init) {
					outfile.open(fileenerg, ios::trunc | ios::out);
				} else {
					outfile.open(fileenerg, ios::app | ios::out);
				}
      
				outfile << gdata.time << " " << Ekfluc << " " << Ebfluc;
				outfile << " " << Ekfluc+Ebfluc << " " << Vorticity << " " << Etherm;
			}
			Problem.writePhystest(gdata, outfile);
			if(gdata.rank == 0) {
				outfile << endl;
				outfile.close();
			}

			Ekin_tave = 0.;
#if (FLUID_TYPE == CRONOS_MHD)
			Emag_tave = 0.;
#endif
			Etherm_tave = 0.;
			Eges_tave = 0.;
			cs_tave = 0.;
			Vorticity_tave = 0.;
			Enstrophy_tave = 0.;
			v2ave_tave = 0.;
			Ekfluc_tave = 0.;
#if (FLUID_TYPE == CRONOS_MHD)
			Ebfluc_tave = 0.;
#endif
			time_ave = 0.;

		}

	}

#if(ENERGETICS == FULL)
	if(!thermal) {
#if(FLUID_TYPE==CRONOS_MULTIFLUID)
		TrafoMulti[iFluid]->TransEth2T(gdata, gfunc, Problem);
#else
		Trafo->TransEth2T(gdata, gfunc, Problem);
#endif
	}
#endif



	// for(int q=0; q<N_OMINT; ++q) {
	// 	if(Problem.force_max(q)) {
	// 		gdata.om[q].set_min(Problem.min_Val(q));
	// 	}
	// 	if(Problem.force_min(q)) {
	// 		gdata.om[q].set_max(Problem.max_Val(q));
	// 	}
	// }

	// Errors will be evaluated every time
	TempNeg = 0;
	int erval(0), evval(0), epval(0);
	int e_nan[N_OMEGA] = {0};
	double Emin(1.e20);

	for (int k = 0; k <= gdata.mx[2]; k++) {
		for (int j = 0; j <= gdata.mx[1]; j++) {
			for (int i = 0; i <= gdata.mx[0]; i++) {

				REAL vabs = sqrt(sqr(gdata.om[q_sx](i,j,k)) +
				                 sqr(gdata.om[q_sy](i,j,k)) +
				                 sqr(gdata.om[q_sz](i,j,k)));
				
				if (gdata.om[q_rho](i,j,k) < rhomin) {
					cout << " Density below zero (" << gdata.om[q_rho](i,j,k)
						  << ") for rank " << gdata.rank
						  << ", ijk=(" << i << " " << j << " " << k << ") , pos=("
						  << gdata.get_x(i,0) << " " << gdata.get_y(j,0) << " "
						  << gdata.get_z(k,0) << ")" << endl;
					// exit(2); // commented out 19jun2019 JK
					erval += 1;
				}
				
				if (vabs > vmax) {
					cout << " Velocity too large ("
						  << gdata.om[q_sx](i,j,k) << " "
						  << gdata.om[q_sy](i,j,k) << " "
						  << gdata.om[q_sz](i,j,k) << ") for rank " << gdata.rank
						  << ", ijk=(" << i << " " << j << " " << k << ") , pos=("
						  << gdata.get_x(i,0) << " " << gdata.get_y(j,0) << " "
						  << gdata.get_z(k,0) << ")" << endl;
					evval += 1;
				}
				
				// If dual energy is used the minimum thermal energy
				// is set to zero.
#if(CRSWITCH_DUAL_ENERGY == CRONOS_ON)
				if(ENERGETICS == FULL && (gdata.om[q_Eges](i,j,k) < 0.)) {
#else
				if(ENERGETICS == FULL && (gdata.om[q_Eges](i,j,k) < Ethmin)) {
#endif
					TempNeg += 1;
					cout << " T < 0 for rank " << gdata.rank
						  << ", ijk=(" << i << " " << j << " " << k << "), "
						  << gdata.om[q_Eges](i,j,k) << endl;
					cout << gdata.om[q_rho](i,j,k) << " ";
					cout << gdata.om[q_sx](i,j,k) << " ";
					cout << gdata.om[q_sy](i,j,k) << " ";
					cout << gdata.om[q_sz](i,j,k) << " ";
					cout << gdata.om[q_Bx](i,j,k) << " ";
					cout << gdata.om[q_By](i,j,k) << " ";
					cout << gdata.om[q_Bz](i,j,k) << " ";
					cout << endl;
					epval += 1;
				}

				for (int q = 0; q < N_OMEGA; q++) {
					if (std::isnan(gdata.om[q](i,j,k))) {
						e_nan[q] += 1;	    
					}
				}
	
				if (ENERGETICS==FULL) {
					Emin = std::min(Emin, gdata.om[q_Eges](i,j,k));
				}
			}
		}
	}

	int errnum(0);
	errnum += erval;
	errnum += evval;
	errnum += epval;
	for (int q = 0; q < N_OMEGA; q++) {
		errnum += e_nan[q];
	}


	// Write error file only if errors did occur:
	if(errnum > 0) {

		char filecd[255];
		char filerr[255];

		sprintf(filecd,"%s",getenv("pname"));
#ifdef parallel
		sprintf(filerr,"%s.err.%2.2d",filecd,gdata.rank);
#else
		sprintf(filerr,"%s.err",filecd);
#endif

		FILE *fo;    
		if ((fo = fopen(filerr,"w")) == NULL) {
			throw CException("cannot open error file!");
		}

		fprintf(fo,"phystest called with the following parameters:\n");
		fprintf(fo," rhomin = %g\n",rhomin);
		fprintf(fo," Ethmin = %g\n",Ethmin);
		fprintf(fo," vmax   = %g\n",vmax);
		fprintf(fo," time   = %g\n",gdata.time);

		if(debug > 0) {
			for (int k = 0; k <= gdata.mx[2]; k++) {
				double zz = gdata.getCen_z(k);
				for (int j = 0; j <= gdata.mx[1]; j++) {
					double yy = gdata.getCen_y(j);
					for (int i = 0; i <= gdata.mx[0]; i++) {
						double xx = gdata.getCen_x(i);
	    
						REAL vabs = sqrt(sqr(gdata.om[q_sx](i,j,k)) +
						                 sqr(gdata.om[q_sy](i,j,k)) +
						                 sqr(gdata.om[q_sz](i,j,k)));
	    
						if (gdata.om[q_rho](i,j,k) < rhomin) {
							fprintf(fo," *** rho < rhomin ***\n");
							fprintf(fo," -> (i,j,k) = (%3d,%3d,%3d)\n",i,j,k);
							fprintf(fo," -> (x,y,z) = (%5.2f,%5.2f,%5.2f)\n",xx,yy,zz);
							fprintf(fo," -> value   = %g\n",gdata.om[q_rho](i,j,k));
						}

						if (vabs > vmax) {
							fprintf(fo," *** v   > vmax ***\n");
							fprintf(fo," -> (i,j,k) = (%3d,%3d,%3d)\n",i,j,k);
							fprintf(fo," -> (x,y,z) = (%5.2f,%5.2f,%5.2f)\n",xx,yy,zz);
							fprintf(fo," -> value   = %g\n",vabs);
							fprintf(fo," -> rho*vx  = %g\n",gdata.om[q_sx](i,j,k));
							fprintf(fo," -> rho*vy  = %g\n",gdata.om[q_sy](i,j,k));
							fprintf(fo," -> rho*vz  = %g\n",gdata.om[q_sz](i,j,k));
							fprintf(fo," -> rho     = %g\n",gdata.om[q_rho](i,j,k));
						}

						for (int q = 0; q < N_OMEGA; q++) {
							if (std::isnan(gdata.om[q](i,j,k))) {
								fprintf(fo,"err: NaN @(ijk) = (%3d %3d %3d)",i,j,k);
								fprintf(fo," # q=%2d\n",q);
								fprintf(fo," -> (x,y,z) = (%5.2f,%5.2f,%5.2f)\n",xx,yy,zz);
							}
						}
	
						if(ENERGETICS == FULL && (gdata.om[q_Eges](i,j,k) < Ethmin)) {
							fprintf(fo," *** Etherm < Ethmin ***\n");
							fprintf(fo," -> (i,j,k) = (%3d,%3d,%3d)\n",i,j,k);
							fprintf(fo," -> (x,y,z) = (%5.2f,%5.2f,%5.2f)\n",xx,yy,zz);
							fprintf(fo," -> value   = %g\n",gdata.om[q_Eges](i,j,k));
						}
	    
					}
				}
			}
		}
		fclose(fo);
	}
	
	bool abortProgram(true);
	int nan_sum(0);
#ifdef parallel
	gdata.MpiSum(TempNeg);
	gdata.MpiSum(erval);
	gdata.MpiSum(evval);
	gdata.MpiSum(epval);

	int e_nanglobal[N_OMEGA];
	MPI_Reduce(&e_nan, &e_nanglobal, N_OMEGA, MPI_INT,
	           MPI_SUM, 0, gdata.comm3d);
	MPI_Bcast(&e_nanglobal, N_OMEGA, MPI_INT, 0, gdata.comm3d);
	MPI_Barrier(gdata.comm3d);

	for (int q = 0; q < N_OMEGA; q++){
		e_nan[q] = e_nanglobal[q];
	}
#endif
	for (int q = 0; q < N_OMEGA; q++) {
		nan_sum += e_nan[q];  // Total NAN-count
	}
	// Total error count
	int esum = erval + epval + evval + nan_sum;


	if(erval > 0){
		if(gdata.rank == 0) {
			cout << " Negative densities: " << erval << endl;
		}
		abortProgram = true;
	}

	if(ENERGETICS == FULL && TempNeg != 0 && gdata.rank == 0){
		cout << " Negative Temps: " << TempNeg << endl;
	}

	if (esum > 0) {
		if(gdata.rank==0) {
			cout << "phystest failed!" << endl;
      
			printf("\nat time = %f:\n",gdata.time);
			if (erval > 0) printf("rho    < rmin at %5d grid points\n",erval);
			if (epval > 0) printf("pgas   < pmin at %5d grid points\n",epval);
			if (evval > 0) printf("v      > vmax at %5d grid points\n",evval);
			for (int q = 0; q < N_OMEGA; q++) {
				if (e_nan[q] > 0) {
					printf("om[%2d] = NaN at %5d grid points\n", q, e_nan[q]);
				}
			}
		}
		abortProgram = true;
	} else {
		abortProgram = false;
	}

	if(abortProgram && gdata.time == 0.) {
		throw CException(" Unphysical initial configuration ");
	} else if (abortProgram) {
		throw CException(" Unphysical configuration ");
	} else {
		return;
	}

	if (abortProgram) {
		exit(2);  // moved here from above  JK
	}
}



double HyperbolicSolver::power(double base, int expo)
{
	double ret = 1.;
	for (int mul = 1; mul <= abs(expo); mul++)
		ret *= base;
	return (expo > 0) ? ret : 1./ret;
}


void HyperbolicSolver::compute_B(Data &gdata, gridFunc &gfunc,
                                 ProblemType &Problem, bool initial_setup) {
  
	if(N_SUBS < 3) {
		cerr << " Need extra array for vector potential " << endl;
		exit(-27);
	}

#if(FLUID_TYPE == CRONOS_MULTIFLUID)
	int q_Ax = n_omIntAll+N_ADD;
	int q_Ay = n_omIntAll+N_ADD+1;
	int q_Az = n_omIntAll+N_ADD+2;

//	Pot & Ax = gdata.om[n_omIntAll+N_ADD];
//	Pot & Ay = gdata.om[n_omIntAll+N_ADD+1];
//	Pot & Az = gdata.om[n_omIntAll+N_ADD+2];
#else
	int q_Ax = n_omInt+N_ADD;
	int q_Ay = n_omInt+N_ADD+1;
	int q_Az = n_omInt+N_ADD+2;

//	Pot & Ax = gdata.om[n_omInt+N_ADD];
//	Pot & Ay = gdata.om[n_omInt+N_ADD+1];
//	Pot & Az = gdata.om[n_omInt+N_ADD+2];
#endif
	Pot & Ax = gdata.om[q_Ax];
	Pot & Ay = gdata.om[q_Ay];
	Pot & Az = gdata.om[q_Az];

#if (VEC_POT_BCS == TRUE)
	// If vector-potential boundaries are needed, we need to
	// (slightly) extend the corresponding fields:
	if(!bcVecPotResized) {

// 		Pot ABack;
// 		ABack.resize(gdata.mx);

// 		ABack = Ax;
// 		Ax.resize(Index::set(-B,-B,-B), 
// 		          Index::set(gdata.mx[0]+B,gdata.mx[1]+B+1,gdata.mx[2]+B+1));
		
// 		for (int k = -B; k <= gdata.mx[2]+B; ++k){
// 			for (int j = -B; j <= gdata.mx[1]+B; ++j){
// 				for (int i = -B; i <= gdata.mx[0]+B; ++i){

// 					ABack(i,j,k) = Ax(i,j,k);

// 				}
// 			}
// 		}

		Ax.resize_ll(Index::set(-B,-B,-B), 
		             Index::set(gdata.mx[0]+B,gdata.mx[1]+B+1,gdata.mx[2]+B+1));
		Ay.resize_ll(Index::set(-B,-B,-B), 
		             Index::set(gdata.mx[0]+B+1,gdata.mx[1]+B,gdata.mx[2]+B+1));
		Az.resize_ll(Index::set(-B,-B,-B), 
		             Index::set(gdata.mx[0]+B+1,gdata.mx[1]+B+1,gdata.mx[2]+B));
		bcVecPotResized = true;
	}



	// Need to compute boundary-conditions for vector-potential first:
	gfunc.boundary(gdata, Problem, gdata.om[q_Ax], B, q_Ax);
	gfunc.boundary(gdata, Problem, gdata.om[q_Ay], B, q_Ay);
	gfunc.boundary(gdata, Problem, gdata.om[q_Az], B, q_Az);

#endif

#ifdef parallel
	gfunc.boundary_MPI(gdata, Problem, gdata.om[q_Ax],B,q_Ax);
	gfunc.boundary_MPI(gdata, Problem, gdata.om[q_Ay],B,q_Ay);
	gfunc.boundary_MPI(gdata, Problem, gdata.om[q_Az],B,q_Az);
#endif


	Pot & Bx = gdata.om[q_Bx];
	Pot & By = gdata.om[q_By];
	Pot & Bz = gdata.om[q_Bz];


#if (VEC_POT_BCS == TRUE)
	int imin[3] = {-B,-B,-B};
	int imax[3] = {gdata.mx[0]+B, gdata.mx[1]+B, gdata.mx[2]+B};
#else
	int imin[3] = {-1, 0, 0};
	int imax[3] = {gdata.mx[0], gdata.mx[1], gdata.mx[2]};

	// Prepare to compute magnetic field within boundaries after initial conditions
	if(initial_setup) {
		Ax.resize_ll(Index::set(-B,-B,-B),
				Index::set(gdata.mx[0]+B,gdata.mx[1]+B+1,gdata.mx[2]+B+1));
		Ay.resize_ll(Index::set(-B,-B,-B),
				Index::set(gdata.mx[0]+B+1,gdata.mx[1]+B,gdata.mx[2]+B+1));
		Az.resize_ll(Index::set(-B,-B,-B),
				Index::set(gdata.mx[0]+B+1,gdata.mx[1]+B+1,gdata.mx[2]+B));
		for(int i_dir=0; i_dir<DIM; i_dir++) {
			imin[i_dir] = -B;
			imax[i_dir] = gdata.mx[i_dir]+B;
		}
	}
#endif

#if (NON_LINEAR_GRID == CRONOS_OFF)
	REAL idx[3]= {(1./gdata.dx[0]), (1./gdata.dx[1]), (1./gdata.dx[2])};
#endif

	for (int k = imin[2]; k <= imax[2]; ++k){
		for (int j = imin[1]; j <= imax[1]; ++j){
			for (int i = imin[0]; i <= imax[0]; ++i){

#if (NON_LINEAR_GRID == CRONOS_ON)
				REAL idx[3] = {(gdata.getCen_idx(0,i)), (gdata.getCen_idx(1,j)),
				               (gdata.getCen_idx(2,k))};
#endif

#if (GEOM > 1)

				REAL f_geom = 1./(gdata.h1(i,j,k,1,0,0)*gdata.h2(i,j,k,1,0,0));
				Bx(i,j,k) = (+(gdata.h2(i,j,k,1, 1,0)*Az(i+1,j+1,k  ) - 
				               gdata.h2(i,j,k,1,-1,0)*Az(i+1,j  ,k  ))*idx[1]
				             -(gdata.h1(i,j,k,1,0, 1)*Ay(i+1,j  ,k+1) -
				               gdata.h1(i,j,k,1,0,-1)*Ay(i+1,j  ,k  ))*idx[2])*f_geom;

#else

				Bx(i,j,k) = (+(Az(i+1,j+1,k  ) - Az(i+1,j  ,k  ))*idx[1]
				             -(Ay(i+1,j  ,k+1) - Ay(i+1,j  ,k  ))*idx[2]);

#endif
				
			}
		}
	}

#if (VEC_POT_BCS == FALSE)
	if(!initial_setup) {
		imin[0] = 0;
		imin[1] = -1;
	}
#endif

	for (int k = imin[2]; k <= imax[2]; ++k){
		for (int j = imin[1]; j <= imax[1]; ++j){
			for (int i = imin[0]; i <= imax[0]; ++i){

#if (NON_LINEAR_GRID == CRONOS_ON)
				REAL idx[3] = {(gdata.getCen_idx(0,i)), (gdata.getCen_idx(1,j)),
				               (gdata.getCen_idx(2,k))};
#endif
        
#if (GEOM > 1)

				REAL f_geom = 1./(gdata.h0(i,j,k,0,1,0)*gdata.h2(i,j,k,0,1,0));
				By(i,j,k) = (+(gdata.h0(i,j,k,0,1, 1)*Ax(i  ,j+1,k+1) -
				               gdata.h0(i,j,k,0,1,-1)*Ax(i  ,j+1,k  ))*idx[2]
				             -(gdata.h2(i,j,k, 1,1,0)*Az(i+1,j+1,k  ) -
				               gdata.h2(i,j,k,-1,1,0)*Az(i  ,j+1,k  ))*idx[0])*f_geom;
        
#else

				By(i,j,k) = (+(Ax(i  ,j+1,k+1) - Ax(i  ,j+1,k  ))*idx[2]
				             -(Az(i+1,j+1,k  ) - Az(i  ,j+1,k  ))*idx[0]);

#endif
			}
		}
	}

#if (VEC_POT_BCS == FALSE)
	if(!initial_setup) {
		imin[1] = 0;
		imin[2] = -1;
	}
#endif

	for (int k = imin[2]; k <= imax[2]; ++k){
		for (int j = imin[1]; j <= imax[1]; ++j){
			for (int i = imin[0]; i <= imax[0]; ++i){
        
#if (NON_LINEAR_GRID == CRONOS_ON)
				REAL idx[3] = {(gdata.getCen_idx(0,i)), (gdata.getCen_idx(1,j)),
				               (gdata.getCen_idx(2,k))};
#endif

#if (GEOM > 1)

				REAL f_geom = 1./(gdata.h0(i,j,k,0,0,1)*gdata.h1(i,j,k,0,0,1));
				Bz(i,j,k) = (+(gdata.h1(i,j,k, 1,0,1)*Ay(i+1,j  ,k+1) -
				               gdata.h1(i,j,k,-1,0,1)*Ay(i  ,j  ,k+1))*idx[0]
				             -(gdata.h0(i,j,k,0, 1,1)*Ax(i  ,j+1,k+1) -
				               gdata.h0(i,j,k,0,-1,1)*Ax(i  ,j  ,k+1))*idx[1])*f_geom;
				// REAL f_geom = 1./(gdata.h0(i,j,k+0.5)*gdata.h1(i,j,k+0.5));
				// Bz(i,j,k) = (+(gdata.h1(i+0.5,j,k+0.5)*Ay(i+1,j  ,k+1) -
				//                gdata.h1(i-0.5,j,k+0.5)*Ay(i  ,j  ,k+1))*gdata.idx[0]
				//              -(gdata.h0(i,j+0.5,k+0.5)*Ax(i  ,j+1,k+1) -
				//                gdata.h0(i,j-0.5,k+0.5)*Ax(i  ,j  ,k+1))*gdata.idx[1])*f_geom;

#else

				Bz(i,j,k) = (+(Ay(i+1,j  ,k+1) - Ay(i  ,j  ,k+1))*idx[0]
				             -(Ax(i  ,j+1,k+1) - Ax(i  ,j  ,k+1))*idx[1]);
//									if(i==0 && j>390 && j<410 && k==0) {
//										cout << " Mag " << gdata.om[q_Bz](i,j,k) << " " << q_Bz << endl;
//									}

#endif
        
			}
		}
	}

	// Store values at boundary points for constant BCs
	if(initial_setup) {
		gfunc.prep_boundaries(gdata, Problem);
	}
//	cout << " Bef BC values at (-1,0,0) ";
//	cout << gdata.om[4](-1,0,0) << " ";
//	cout << gdata.om[5](-1,0,0) << " ";
//	cout << gdata.om[6](-1,0,0) << " ";
//	cout << endl;
//	cout << " By at upper border " << " ";
//	cout << gdata.om[5](gdata.mx[0]-1,0,0) << " ";
//	cout << gdata.om[5](gdata.mx[0],0,0) << " ";
//	cout << gdata.mx[0] << " ";
//	cout << endl;
//	cout << " A at upper border " << " ";
//	cout << Ax(gdata.mx[0],1,1) << " ";
//	cout << Ax(gdata.mx[0],1,0) << " ";
//	cout << Az(gdata.mx[0]+1,1,0) << " ";
//	cout << Az(gdata.mx[0],1,0) << " ";
//	cout << Az(gdata.mx[0]-1,1,0) << " ";
//	cout << endl;

	gfunc.boundary(gdata, Problem, gdata.om[q_Bx], B, q_Bx);
	gfunc.boundary(gdata, Problem, gdata.om[q_By], B, q_By);
	gfunc.boundary(gdata, Problem, gdata.om[q_Bz], B, q_Bz);

//	cout << " Aft BC values at (-1,0,0) ";
//	cout << gdata.om[4](-1,0,0) << " ";
//	cout << gdata.om[5](-1,0,0) << " ";
//	cout << gdata.om[6](-1,0,0) << " ";
//	cout << endl;


}



void HyperbolicSolver::TransPot2Mag(Data &gdata, gridFunc &gfunc,
                                    ProblemType &Problem) {
  
	if(gdata.om[q_Bx].getName() == "B_x" || 
	   gdata.om[q_By].getName() == "B_y" ||
	   gdata.om[q_Bz].getName() == "B_z") {
		cerr << " Magnetic field already set " << endl;
		exit(2);
	}

	if(N_SUBS < 3) {
		cerr << " Need extra array for vector potential " << endl;
		exit(-27);
	}

#if(FLUID_TYPE == CRONOS_MULTIFLUID)
	Pot & A_x = gdata.om[n_omIntAll+N_ADD  ];
	Pot & A_y = gdata.om[n_omIntAll+N_ADD+1];
	Pot & A_z = gdata.om[n_omIntAll+N_ADD+2];
#else
	Pot & A_x = gdata.om[n_omInt+N_ADD  ];
	Pot & A_y = gdata.om[n_omInt+N_ADD+1];
	Pot & A_z = gdata.om[n_omInt+N_ADD+2];
#endif

	// Store type of BCs:
	int A_x_bc_TypeLow[3]  = {A_x.get_bcTypeLow(0), A_x.get_bcTypeLow(1), A_x.get_bcTypeLow(2)};
	int A_x_bc_TypeHigh[3] = {A_x.get_bcTypeHigh(0), A_x.get_bcTypeHigh(1), A_x.get_bcTypeHigh(2)};
	int A_y_bc_TypeLow[3]  = {A_y.get_bcTypeLow(0), A_y.get_bcTypeLow(1), A_y.get_bcTypeLow(2)};
	int A_y_bc_TypeHigh[3] = {A_y.get_bcTypeHigh(0), A_y.get_bcTypeHigh(1), A_y.get_bcTypeHigh(2)};
	int A_z_bc_TypeLow[3]  = {A_z.get_bcTypeLow(0), A_z.get_bcTypeLow(1), A_z.get_bcTypeLow(2)};
	int A_z_bc_TypeHigh[3] = {A_z.get_bcTypeHigh(0), A_z.get_bcTypeHigh(1), A_z.get_bcTypeHigh(2)};

	A_x = gdata.om[q_Bx];
	A_y = gdata.om[q_By];
	A_z = gdata.om[q_Bz];

	// Restort type of BCs
	A_x.set_bcType(A_x_bc_TypeLow, A_x_bc_TypeHigh);
	A_y.set_bcType(A_y_bc_TypeLow, A_y_bc_TypeHigh);
	A_z.set_bcType(A_z_bc_TypeLow, A_z_bc_TypeHigh);

//	cout << " BC Type for " << n_omInt+N_ADD << " " << A_x.get_bcTypeLow(0) << " ";
//	cout << A_x_bc_TypeLow[0] << " ";
//	cout << gdata.om[q_Bx].get_bcTypeLow(0) << endl;

//#if(FLUID_TYPE == CRONOS_MULTIFLUID)
//	gdata.om[n_omIntAll+N_ADD  ] = gdata.om[q_Bx];
//	gdata.om[n_omIntAll+N_ADD+1] = gdata.om[q_By];
//	gdata.om[n_omIntAll+N_ADD+2] = gdata.om[q_Bz];
//#else
//
//	gdata.om[n_omInt+N_ADD  ] = gdata.om[q_Bx];
//	gdata.om[n_omInt+N_ADD+1] = gdata.om[q_By];
//	gdata.om[n_omInt+N_ADD+2] = gdata.om[q_Bz];
//#endif

	// Renaming need to be done first because BCs in compute_B need to
	// know that they are dealing with B-field
	gdata.om[q_Bx].rename("B_x");
	gdata.om[q_By].rename("B_y");
	gdata.om[q_Bz].rename("B_z");
	
	// Indicate that this is first computation of magnetic field
	compute_B(gdata, gfunc, Problem, true);


}




void HyperbolicSolver::TransMag2Pot(Data &gdata, gridFunc &gfunc,
                                    ProblemType &Problem) {
	if(gdata.om[q_Bx].getName() == "A_x" || 
	   gdata.om[q_By].getName() == "A_y" ||
	   gdata.om[q_Bz].getName() == "A_z") {
		cerr << " Magnetic vector potential already set " << endl;
		exit(2);
	}


	REAL ChecksumAloc[3];
	//   REAL eps(1.e-14*(std::max(ChecksumA[0],std::max(ChecksumA[1],ChecksumA[2])) + 1.));
	REAL eps(1.e-10*(std::max(ChecksumA[0],std::max(ChecksumA[1],ChecksumA[2])) + 1.));
	//  REAL eps(1.e-14);
	ComputeChecksumA(gdata, ChecksumAloc);
	for(int dd=0; dd<DIM; ++dd) {
		if((ChecksumA[dd] > ChecksumAloc[dd]+eps) || 
		   (ChecksumA[dd] < ChecksumAloc[dd]-eps)) {
			cout << ChecksumA[dd] << " " << ChecksumAloc[dd] << " ";
			cout << ChecksumA[dd] - ChecksumAloc[dd] << " ";
			cout << dd << endl;
			cerr << " Vector Potential was corrupted -- exiting " << endl;
			exit(-56);
		}
	}
  
#if(FLUID_TYPE == CRONOS_MULTIFLUID)
	gdata.om[q_Bx] = gdata.om[n_omIntAll+N_ADD];
	gdata.om[q_By] = gdata.om[n_omIntAll+N_ADD+1];
	gdata.om[q_Bz] = gdata.om[n_omIntAll+N_ADD+2];
#else
	cerr << " specific::TransMag2Pot -> this looks wrong " << endl;
	exit(3);
	gdata.om[q_Bx] = gdata.om[n_omInt+1];
	gdata.om[q_By] = gdata.om[n_omInt+2];
	gdata.om[q_Bz] = gdata.om[n_omInt+3];
#endif

	gdata.om[q_Bx].rename("A_x");
	gdata.om[q_By].rename("A_y");
	gdata.om[q_Bz].rename("A_z");
  
}

void HyperbolicSolver::ComputeChecksumA(Data &gdata, REAL CheckSum[3]) {

	CheckSum[0] = 0.;
	CheckSum[1] = 0.;
	CheckSum[2] = 0.;
	for (int k = 0; k <= gdata.mx[2]; ++k){
		for (int j = 0; j <= gdata.mx[1]; ++j){
			for (int i = 0; i <= gdata.mx[0]; ++i){
#if(FLUID_TYPE == CRONOS_MULTIFLUID)
				CheckSum[0] += gdata.om[n_omInt+N_ADD](i,j,k);
				CheckSum[1] += gdata.om[n_omInt+N_ADD+1](i,j,k);
				CheckSum[2] += gdata.om[n_omInt+N_ADD+2](i,j,k);
#else
				CheckSum[0] += gdata.om[n_omInt+1](i,j,k);
				CheckSum[1] += gdata.om[n_omInt+2](i,j,k);
				CheckSum[2] += gdata.om[n_omInt+3](i,j,k);
#endif
			}
		}
	}  

}



void HyperbolicSolver::CheckNeg(NumMatrix<REAL,3> &field, int q, int rim)
{
	int lo[3], up[3];
	lo[0] = field.getLow(0);
	lo[1] = field.getLow(1);
	lo[2] = field.getLow(2);

	up[0] = field.getHigh(0);
	up[1] = field.getHigh(1);
	up[2] = field.getHigh(2);

	for (int k = lo[2]; k <= up[2]; k++){
		for (int j = lo[1]; j <= up[1]; j++){
			for (int i = lo[0]; i <= up[0]; i++){
				if(field(i,j,k) < 0) {
					throw CException("Field is < 0",q,i,j,k);
				}
			}
		}
	}
}



void HyperbolicSolver::CheckNan(NumMatrix<REAL,3> &field, 
                                int q, int rim, int pos,
                                string fieldname)
{
	string message = fieldname + " is NAN ";

	int lo[3], up[3];
	lo[0] = field.getLow(0);
	lo[1] = field.getLow(1);
	lo[2] = field.getLow(2);

	up[0] = field.getHigh(0);
	up[1] = field.getHigh(1);
	up[2] = field.getHigh(2);

	for (int k = lo[2]; k <= up[2]; k++){
		for (int j = lo[1]; j <= up[1]; j++){
			for (int i = lo[0]; i <= up[0]; i++){
				if(std::isnan(field(i,j,k))) {
					throw CException(message,pos,q,i,j,k);
				}
			}
		}
	}
}




