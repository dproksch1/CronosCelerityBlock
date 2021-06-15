#include "gridgen.H"
#include "queue.H"

void Environment::init_solvers(Data &gdata)
{
	// Specified user class
	RKSolver = std::make_unique<HyperbolicSolver>(gdata, *Problem);
//	RKSolver = new HyperbolicSolver(gdata, *Problem);
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


void Environment::pdestep(Data &gdata, Queue& queue)
{
	gdata.cfl = 0.;

	// REAL timeinit(gdata.time);

#ifdef TIME_INTEGRATOR

	for (int n = 0; n < TIME_SUBSTEPS; n++) {
		if(gdata.rank==0 && Problem->checkout(0)){
			cout << "  RKSTEP = " << n+1 << endl;
		}
		try {
			gdata.cfl = (std::max(RKSolver->singlestep(gdata, *gfunc, *Problem, n, queue),
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
