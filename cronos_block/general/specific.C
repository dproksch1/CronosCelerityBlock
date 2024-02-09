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

	Trafo = std::make_unique<Transformations_Block>(gdata.fluid, Problem, false);

	IntegrateA = true;
	bcVecPotResized = false;

	debug = static_cast<int>(value((char*)"debug"));
	
	string RiemannSolver = svalue((char*)"RiemannSolver");

	Riemann.resize(DirMax);

	#if (GEOM != CARTESIAN)
		assert_fail_lazy() << "Non-Cartesian geometry is not implemented for this version";
	#endif

	#if (USE_COROTATION == CRONOS_ON)
		assert_fail_lazy() << "Corotating frame of reference is not implemented for this version";
	#endif

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
	// PhysFlux = std::make_unique<PhysFluxesHD>(gdata, gdata.fluid);
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
}


HyperbolicSolver::~HyperbolicSolver() {
  
/*#ifdef SAVEMEM
	if(Save != NULL) {
		delete Save;
	}
#endif*/

#if (OMS_USER == TRUE)
	if(TimeIntegratorUser != NULL) {
		delete [] TimeIntegratorUser;
	}  
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


void HyperbolicSolver::init(Queue &queue, Data &gdata, gridFunc &gfunc,
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


	init_general(queue, gdata, gfunc, Problem);
}


void HyperbolicSolver::init_general(Queue &queue, Data &gdata, gridFunc &gfunc,
                                    ProblemType &Problem)
{
	// General initialisation procedure
	// This routine is executed for normal start AND for a restart
	
	int i_magFluid = 0;

	// Introducing unique ids to all fields om, om_user - can be used
	// for identification (like, e.g., in bcs)
	gdata.set_fieldIds();

	set_TimeIntegrator(queue, gdata, gfunc);


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

}


void HyperbolicSolver::set_TimeIntegrator(Queue &queue, const Data &gdata,
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

	TimeIntegratorGeneric[0]->init_omBuffer(queue, gdata.mx);

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


void HyperbolicSolver::restart(Queue &queue, Data &gdata, gridFunc &gfunc,
                               ProblemType &Problem)
{
    // Store boundary values for fixed boundary conditions
	gfunc.prep_boundaries(gdata, Problem);
	init_general(queue, gdata, gfunc, Problem);
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