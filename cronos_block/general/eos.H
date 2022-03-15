#ifndef CRONOS_EOS_H
#define CRONOS_EOS_H = 1

#include "problem.H"

/*
  Class supplying an equation of state with a speed of sound to be
  defined in the ProblemType-Class.
*/

class EquationOfState {
public:
	EquationOfState(const ProblemType &);
	inline double pressure(const double &, const double &, const double &,
	                     const double &);
	inline double pressure(const Data &, const ProblemType &, double, double, double, double) const;
protected:
	double eps, gamma, half_beta, denominator;
};




inline double EquationOfState::pressure(const Data &gdata, const ProblemType &Problem,
                                      double dens, double i, double j, double k) const
{
#ifdef SOUND
	if(gamma != 1.){
		return Problem.c2_iso(gdata,i,j,k)*pow(dens, gamma)*denominator;
	} else {
		return Problem.c2_iso(gdata,i,j,k)*dens;
	}
#else
	if(gamma != 1.){
		return Problem.c2_iso(gdata,i,j,k)*pow(dens, gamma)*denominator*half_beta;
	} else {
		return Problem.c2_iso(gdata,i,j,k)*dens*half_beta;
	}
#endif
}


inline double EquationOfState::pressure(const double &rhoinv, const double &psq,
                                      const double &Bsq, const double &Energy)
{
	assert(gamma != 1);
	double PresVal= (gamma-1)*(Energy - 0.5*psq*rhoinv - 0.5*Bsq);

	if (PresVal < 0.){
		PresVal = eps;
	}
	return PresVal;
}

#endif
