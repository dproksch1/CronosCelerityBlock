#include "gridgen.H"

void Environment::init_solvers(Data &gdata)
{
	// Specified user class
	RKSolver = new HyperbolicSolver(gdata, *Problem);
	try{
		RKSolver->init_constants(gdata);
	} catch (CException exep) {
		Abort(gdata, exep);
	}

	// EulerSolver = new ISMCooling(gdata);
	// try{
	// 	EulerSolver->init_constants(gdata, *gfunc, *Problem);
	// } catch (CException exep) {
	// 	Abort(gdata, exep);
	// }

}


void Environment::pdestep(Data &gdata)
{
	gdata.cfl = 0.;

	// REAL timeinit(gdata.time);

#ifdef RK_STEPS

	for (int n = 0; n < RK_STEPS; n++) {
		if(gdata.rank==0 && Problem->checkout(0)){
			cout << "  RKSTEP = " << n+1 << endl;
		}
		try {
			gdata.cfl = (std::max(RKSolver->singlestep(gdata, *gfunc, *Problem, n),
			                      gdata.cfl));
		} catch (CException exep) {

			Abort(gdata, exep);
		}

	}

#endif

  // gdata.time = timeinit;
  // gdata.cfl = (std::max(EulerSolver->singlestep(gdata, *gfunc, *Problem),
  // 		       gdata.cfl));
  // gdata.time += gdata.dt;

}
