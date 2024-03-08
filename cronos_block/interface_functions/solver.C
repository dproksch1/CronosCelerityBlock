#include "solver.H"

#include "gridgen.H"

#define CR_CELERITY CELERITY_OFF


//! @brief Initalizes the solver and its constants
void Environment::init_solvers(Data &gdata)
{
	// Specified user class
	rksolver = std::make_unique<HyperbolicSolver>(gdata, *Problem);

	try{
		rksolver->init_constants(gdata);
	} catch (CException exep) {
		Abort(gdata, exep);
	}

}


//! @brief Performs a single timestep using the partial-differential equation
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
			rksolver->singlestep(gdata, *gfunc, *Problem, n, queue);
			// gdata.cfl = (std::max(rksolver->singlestep(gdata, *gfunc, *Problem, n, queue), gdata.cfl));
		} catch (CException exep) {

			Abort(gdata, exep);
		}
	}

#endif

}
