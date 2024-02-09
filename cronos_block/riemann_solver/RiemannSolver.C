#include "RiemannSolver.H"
#include <iostream>

using namespace std;

double RiemannSolver::getAlphaCarbuncle() {
	return alpha_carbuncle;
}


void RiemannSolver::set_verbosity(int _verb) {
	verbosity = _verb;
}

RiemannSolverHD::RiemannSolverHD(const Data &gdata, int dir, int Fluid_Type) : RiemannSolver(gdata, dir, Fluid_Type) {
#if(FLUID_TYPE==CRONOS_HYDRO)
	this->q_rho = gdata.fluid.get_q_rho();
	this->q_sx  = gdata.fluid.get_q_sx();
	this->q_sy  = gdata.fluid.get_q_sy();
	this->q_sz  = gdata.fluid.get_q_sz();
	this->q_Eges = gdata.fluid.get_q_Eges();
	this->q_Eadd = gdata.fluid.get_q_Eadd();
#endif
}

#if (ENERGETICS != ISOTHERMAL)
HLLCSolver_Hydro::HLLCSolver_Hydro(const Data &gdata, int dir, int Fluid_Type) : RiemannSolverHD(gdata, dir, Fluid_Type) {
	veps = 1.e-120;

#if(FLUID_TYPE != CRONOS_MULTIFLUID)
	reset_Indices(gdata.fluid);
#else
	reset_Indices(gdata.fluids->fluids[0]);
#endif
	if(dir == 0) {
		this->qvPar = q_sx;
		this->qvP1  = q_sy;
		this->qvP2  = q_sz;
	} else if (dir == 1) {
		this->qvP2  = q_sx;
		this->qvPar = q_sy;
		this->qvP1  = q_sz;
	} else {
		this->qvP1  = q_sx;
		this->qvP2  = q_sy;
		this->qvPar = q_sz;
	}   
	gamma = value((char*)"Adiabatic_exponent");

}

#endif