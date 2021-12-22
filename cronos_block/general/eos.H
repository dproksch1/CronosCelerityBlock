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
	inline REAL pressure(const REAL &, const REAL &, const REAL &,
	                     const REAL &);
	inline REAL pressure(const Data &, const ProblemType &, REAL, REAL, REAL, REAL) const;
protected:
	REAL eps, gamma, half_beta, denominator;
};




inline REAL EquationOfState::pressure(const Data &gdata, const ProblemType &Problem,
                                      REAL dens, REAL i, REAL j, REAL k) const
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


inline REAL EquationOfState::pressure(const REAL &rhoinv, const REAL &psq,
                                      const REAL &Bsq, const REAL &Energy)
{
	assert(gamma != 1);
	REAL PresVal= (gamma-1)*(Energy - 0.5*psq*rhoinv - 0.5*Bsq);

	if (PresVal < 0.){
		PresVal = eps;
	}
	return PresVal;
}

#endif
