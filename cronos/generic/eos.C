#include "eos.H"

EquationOfState::EquationOfState(const ProblemType &Problem) {

	this->gamma = Problem.gamma;
	this->eps   = 1.e-12;
	this->denominator = pow(Problem.rho0,1.-gamma);
	if(Problem.mag) {
		this->half_beta = value((char*)"Plasma_beta")/2.;
	} else {
		this->half_beta = 1;
	}
}

