#include "gridgen.H"

void Environment::init_solvers(Data &gdata)
{
	RKSolver = unique_ptr< HyperbolicSolver > ( new HyperbolicSolver(gdata, *Problem));
	try{
		RKSolver->init_constants(gdata);
	} catch (CException exep) {
		Abort(gdata, exep);
	}
}


void Environment::pdestep(Data &gdata)
{
	gdata.cfl = 1.;
	gdata.time += gdata.dt;
}
