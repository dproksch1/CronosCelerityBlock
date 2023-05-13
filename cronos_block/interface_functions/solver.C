#include "solver.H"

#include "gridgen.H"

#define CR_CELERITY CELERITY_OFF

void Environment::init_solvers(Data &gdata)
{
	// Specified user class
	RKSolver = std::make_unique<HyperbolicSolver>(gdata, *Problem);

	try{
		RKSolver->init_constants(gdata);
	} catch (CException exep) {
		Abort(gdata, exep);
	}

}


void Environment::pdestep(Data &gdata, Queue& queue)
{
	gdata.cfl = 0.;

	// double timeinit(gdata.time);

#ifdef TIME_INTEGRATOR

	for (int n = 0; n < TIME_SUBSTEPS; n++) {
		if(gdata.rank==0 && Problem->checkout(0)){
			cout << "  RKSTEP = " << n+1 << endl;
		}
		try {
			RKSolver->singlestep(gdata, *gfunc, *Problem, n, queue);
			// gdata.cfl = (std::max(RKSolver->singlestep(gdata, *gfunc, *Problem, n, queue), gdata.cfl));
		} catch (CException exep) {

			Abort(gdata, exep);
		}
	}

#endif

}
