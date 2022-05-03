#include "specific.H"

#include <iomanip>

#include "gridgen.H"
#include "utils.H"


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

	init_constants(gdata);

	Trafo = std::make_unique<Transformations>(gdata.fluid, Problem, false);

	IntegrateA = true;
	bcVecPotResized = false;

	debug = static_cast<int>(value((char*)"debug"));
	
	string RiemannSolver = svalue((char*)"RiemannSolver");

	Riemann.resize(DirMax);

    //exclude more solver if unavailable in proxyapp
	if(RiemannSolver == "hll") {
		if(gdata.rank==0) {
			cout << " Using the HLL solver " << endl;
		}
		assert_fail_lazy() << "not implemented";

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
		assert_fail_lazy() << "not implemented";
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

	/*ReconstX.resize(TIME_SUBSTEPS);// = new Reconstruction * [TIME_SUBSTEPS];
	ReconstY.resize(TIME_SUBSTEPS);// = new Reconstruction * [TIME_SUBSTEPS];
	ReconstZ.resize(TIME_SUBSTEPS);// = new Reconstruction * [TIME_SUBSTEPS];

	for (int i=0; i < TIME_SUBSTEPS; i++) {
		ReconstX[i] = std::make_unique<Reconstruction>(gdata, 0, gdata.fluid, i);
		ReconstY[i] = std::make_unique<Reconstruction>(gdata, 1, gdata.fluid, i);
		ReconstZ[i] = std::make_unique<Reconstruction>(gdata, 2, gdata.fluid, i);
	}*/


#if (FLUID_TYPE == CRONOS_HYDRO)
	PhysFlux = std::make_unique<PhysFluxesHD>(gdata, gdata.fluid);
#endif

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
  
/*#ifdef SAVEMEM
	if(Save != NULL) {
		delete Save;
	}
#endif*/

#if (OMS_USER == TRUE)
	if(TimeIntegratorUser != NULL) {
		TimeIntegratorUser.reset();
	}  
#endif

	if (TimeIntegratorGeneric != NULL) {
		TimeIntegratorGeneric.reset();
	}
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


void HyperbolicSolver::addConstants_toH5(Data &gdata, Hdf5Stream &h5out)
{
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

void HyperbolicSolver::getConstants_fromH5(Data &gdata, Hdf5iStream &h5in)
{
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
#endif


		if(gdata.rank==0) {
			cout << "all ICs set." << endl;
		}

	}

	// Store values at boundaries for fixed boundary conditions
	gfunc.prep_boundaries(gdata, Problem);

	// if Vector-Potential is given -> Transform to magnetic field
	/*if (FLUID_TYPE == CRONOS_MHD || with_mag) {
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
	}*/


	init_general(gdata, gfunc, Problem);
}


void HyperbolicSolver::init_general(Data &gdata, gridFunc &gfunc,
                                    ProblemType &Problem)
{
	// General initialisation procedure
	// This routine is executed for normal start AND for a restart
	
	int i_magFluid = 0;

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
	Trafo->set_thermal(thermal);
#endif

	gfunc.boundary(gdata, Problem); 

#if (USE_COROTATION == CRONOS_ON)
		// Transform inertial frame for phystest
		Trafo->TransCorotToInert(gdata, gfunc, Problem);
#endif

	phystest(gdata, gfunc, Problem, -1);

#if (USE_COROTATION == CRONOS_ON)
	// Transform back to co-rotating frame
	Trafo->TransInertToCorot(gdata, gfunc, Problem);
#endif

}


void HyperbolicSolver::set_TimeIntegrator(const Data &gdata,
                                          gridFunc &gfunc) {

	TimeIntegratorGeneric.resize(n_omInt);

	for(int iom=0; iom<n_omInt; iom++) {

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
					qchange = n_omInt+N_ADD+q-q_Bx;
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


void HyperbolicSolver::restart(Data &gdata, gridFunc &gfunc,
                               ProblemType &Problem)
{
    // Store boundary values for fixed boundary conditions
	gfunc.prep_boundaries(gdata, Problem);
	init_general(gdata, gfunc, Problem);
}





double HyperbolicSolver::compute_divB(Data &gdata, gridFunc &gfunc,
                                    ProblemType &Problem)
{
	if(Problem.get_Info() && Problem.checkout(4)) {

		double divVal;
		if(N_ADD >= 2) {
			Pot & divB = gdata.om[n_omInt+1];
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



double HyperbolicSolver::compute_divB(Data &gdata, gridFunc &gfunc,
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
	double idx[3]= {(1./gdata.dx[0]), (1./gdata.dx[1]), (1./gdata.dx[2])};
#endif
	for (int k = min[2]; k <= max[2]; k++) {   // compute J = rot(B) and 
		for (int j = min[1]; j <= max[1]; j++) {   // div(B) as om fields
			for (int i = min[0]; i <= max[0]; i++) {

#if (NON_LINEAR_GRID == CRONOS_ON)
				double idx[3] = {(gdata.getCen_idx(0,i)), (gdata.getCen_idx(1,j)),
				               (gdata.getCen_idx(2,k))};
#endif


#if GEOM > 1
				double f_geom = 1./gdata.get_CellGeomTrafo(i,j,k);
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

	gfunc.boundary(gdata, Problem, divB,3);

	if(gdata.rank==0){
		divBSum /= BSum;
		cout << " Ave: " << divBSum << " " << divBSum*BSum/JSum << endl;
		//   cout << "Sum of B: " << BSum << endl;
		//     cout << "---------------------------------------------" << endl;
	}

	double divVal =  divBSum;
	return divVal;
}


void HyperbolicSolver::phystest(Data &gdata, gridFunc &gfunc,
                                ProblemType &Problem, int rkstep, int iFluid)
{
	bool init(false);
	if(rkstep < 0) {
		init = true;
		rkstep = TIME_SUBSTEPS-1;
	}
  
	double dV(gdata.get_CellVolume(0,0,0));
	double Vol(gdata.get_Volume());
	double Mass(0.), Eges(0.);
	double mv[3] = {0., 0., 0.};
	double Ekin(0.), Emag(0.), Etherm(0.);
	double cs(0.);
	double Ekfluc(0.), Ebfluc(0.);
	double Vorticity(0.), Enstrophy(0.);
	double sqrdivv(0.), sqrdissv(0.);
	double v2ave(0.);
	double mvrms[3] = {0., 0., 0.};
	int q_min = 0;

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
			Trafo->TransMomen2Vel(gdata, gfunc, Problem);
		}
		if(thermal) {
			Trafo->TransEth2E(gdata, gfunc, Problem);
		} else {
			Trafo->TransT2E(gdata, gfunc, Problem);
		}

	}
#endif
	if(gdata.om[q_sx].getName() == "v_x" && gdata.om[q_sy].getName() == "v_y" &&
	   gdata.om[q_sz].getName() == "v_z") {
		Trafo->TransVel2Momen(gdata, gfunc, Problem);
	}

	// Compute conservative quantities:
	Mass = gdata.computeInt(q_rho);

	mv[0] = gdata.computeInt(q_sx);
	mv[1] = gdata.computeInt(q_sy);
	mv[2] = gdata.computeInt(q_sz);  

	if(ENERGETICS == FULL) {
		Eges = gdata.computeInt(q_Eges);
	}
  
	mvrms[0] = gdata.computeRMS(q_sx);
	mvrms[1] = gdata.computeRMS(q_sy);
	mvrms[2] = gdata.computeRMS(q_sz);

	// Making certain to use primitve Variables for the rest of the
	// computations

	if(gdata.om[q_sx].getName() == "s_x" && gdata.om[q_sy].getName() == "s_y" &&
	   gdata.om[q_sz].getName() == "s_z") {
		Trafo->TransMomen2Vel(gdata, gfunc, Problem);
	}
// if(ENERGETICS == FULL) {
	// 	if(thermal) {
	// 		Trafo->TransE2Eth(gdata, gfunc, Problem);
	// 	} else {
	// 		Trafo->TransE2T(gdata, gfunc, Problem);
	// 	}
	// }
	if(ENERGETICS == FULL) {
		Trafo->TransE2Eth(gdata, gfunc, Problem);
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
    
		cs = cs/Vol;
        
		if(ENERGETICS != FULL) {
			Eges = Ekin;
		}
    
    
		// Set initial values at first step:
		if(init) {
			Mass0 = Mass;
			mv0[0] = mv[0];
			mv0[1] = mv[1];
			mv0[2] = mv[2];
			Eges0 = Eges;
		}
    
		// Saving time-averaged data:
		Ekin_tave += gdata.dt*Ekin;
		Etherm_tave += gdata.dt*Etherm;
		Eges_tave += gdata.dt*Eges;
		cs_tave += gdata.dt*cs;
		time_ave += gdata.dt;


		// Compute Velocity Stats:
		Vol = 0.;
#if (NON_LINEAR_GRID == CRONOS_OFF)
		double hx[3]= {(gdata.hx[0]), (gdata.hx[1]), (gdata.hx[2])};
#endif
		for (int k = 0; k < gdata.mx[2]; k++) {
			for (int j = 0; j < gdata.mx[1]; j++) {
				for (int i = 0; i < gdata.mx[0]; i++) {

#if ((GEOM > 1) || (NON_LINEAR_GRID == CRONOS_ON))
					dV = gdata.get_CellVolume(i,j,k);
#endif
#if (NON_LINEAR_GRID == CRONOS_ON)
					double hx[3]= {(gdata.getCen_hx(0,i)), (gdata.getCen_hx(1,j)),
					             (gdata.getCen_hx(2,k))};	
#endif

	  
#if (GEOM > 1)
					double f_geom = 1./gdata.get_CellGeomTrafo(i,j,k);

					double dxVx = (gdata.getCen_h1(i+1,j,k)*gdata.getCen_h2(i+1,j,k)*gdata.om[q_sx](i+1,j,k) -
					             gdata.getCen_h1(i-1,j,k)*gdata.getCen_h2(i-1,j,k)*gdata.om[q_sx](i-1,j,k))*hx[0]*f_geom;
					double dxVy = (gdata.getCen_h1(i+1,j,k)*gdata.om[q_sy](i+1,j,k) -
					             gdata.getCen_h1(i-1,j,k)*gdata.om[q_sy](i-1,j,k))*hx[0];
					double dxVz = (gdata.getCen_h2(i+1,j,k)*gdata.om[q_sz](i+1,j,k) -
					             gdata.getCen_h2(i-1,j,k)*gdata.om[q_sz](i-1,j,k))*hx[0];
					double dyVx = (gdata.getCen_h0(i,j+1,k)*gdata.om[q_sx](i,j+1,k) -
					             gdata.getCen_h0(i,j-1,k)*gdata.om[q_sx](i,j-1,k))*hx[1];
					double  dyVy = (gdata.getCen_h0(i,j+1,k)*gdata.getCen_h2(i,j+1,k)*gdata.om[q_sy](i,j+1,k) -
					              gdata.getCen_h0(i,j-1,k)*gdata.getCen_h2(i,j-1,k)*gdata.om[q_sy](i,j-1,k))*hx[1]*f_geom;
					double dyVz = (gdata.getCen_h2(i,j+1,k)*gdata.om[q_sz](i,j+1,k) -
					             gdata.getCen_h2(i,j-1,k)*gdata.om[q_sz](i,j-1,k))*hx[1];
					double dzVx = (gdata.getCen_h0(i,j,k+1)*gdata.om[q_sx](i,j,k+1) -
					             gdata.getCen_h0(i,j,k-1)*gdata.om[q_sx](i,j,k-1))*hx[2];
					double dzVy = (gdata.getCen_h1(i,j,k+1)*gdata.om[q_sy](i,j,k+1) -
					             gdata.getCen_h1(i,j,k-1)*gdata.om[q_sy](i,j,k-1))*hx[2];
					double dzVz = (gdata.getCen_h0(i,j,k+1)*gdata.getCen_h1(i,j,k+1)*gdata.om[q_sz](i,j,k+1) -
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
					double curl_x = (dyVz - dzVy)/(gdata.getCen_h1(i,j,k)*gdata.getCen_h2(i,j,k));
					double curl_y = (dzVx - dxVz)/(gdata.getCen_h0(i,j,k)*gdata.getCen_h2(i,j,k));
					double curl_z = (dxVy - dyVx)/(gdata.getCen_h0(i,j,k)*gdata.getCen_h1(i,j,k));
#else
					double curl_x = (dyVz - dzVy);
					double curl_y = (dzVx - dxVz);
					double curl_z = (dxVy - dyVx);
#endif
					double curl_val = (sqr(curl_x) + sqr(curl_y) + sqr(curl_z));
	    
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

		v2ave /= Vol;
    
		Problem.computeFluct(gdata, Ekfluc, Ebfluc);
    
		// Compute temporal averages (in case no output at every timestep
		Vorticity_tave += gdata.dt*Vorticity;
		Enstrophy_tave += gdata.dt*Enstrophy;
		v2ave_tave += gdata.dt*v2ave;
		Ekfluc_tave += gdata.dt*Ekfluc;

		// Compute quantities in derived class (if necessary)
		Problem.computePhystest(gdata);
    
		// Writing to screen -- if desired
		if(Problem.get_Info() && Problem.checkout(3)) {
			double MachAve = sqrt(v2ave)/cs;
    
			if(gdata.rank == 0) {
				cout << setiosflags(ios::scientific) << setprecision(3) << setw(12);
				cout << "-------Fluctuations-----------------------------------" << endl;
				cout << " E_kin:  " << setw(9) << Ekfluc;
				cout << " E_dyn:  " << setw(9) << Ekfluc + Ebfluc << endl;
				cout << "-------Energertics:-----------------------------------" << endl;
				cout << " E_ges:  " << setw(9) << Eges;
				if(ENERGETICS == FULL){
					cout << " E_th:   " << setw(9) << Etherm;
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
			Etherm_tave = 0.;
			Eges_tave = 0.;
			cs_tave = 0.;
			Vorticity_tave = 0.;
			Enstrophy_tave = 0.;
			v2ave_tave = 0.;
			Ekfluc_tave = 0.;
			time_ave = 0.;

		}

	}

#if(ENERGETICS == FULL)
	if(!thermal) {
		Trafo->TransEth2T(gdata, gfunc, Problem);
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

				double vabs = sqrt(sqr(gdata.om[q_sx](i,j,k)) +
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
		sprintf(filerr,"%s.err",filecd);

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
	    
						double vabs = sqrt(sqr(gdata.om[q_sx](i,j,k)) +
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


void HyperbolicSolver::CheckNan(NumMatrix<double,3> &field, 
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