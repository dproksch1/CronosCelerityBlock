#ifndef CRONOS_SUPERBEE_H
#define CRONOS_SUPERBEE_H = 1
#include <iostream>

#include "constants.H"

using namespace std;

class Superbee /*: public limiter_func*/ {
public:
	Superbee() {
		sigma1 = sigma2 = 0.;
	}
	REAL compute(REAL&, REAL&, REAL&);
private:
	REAL minmod(double aVal, double bVal);
	REAL maxmod(double aVal, double bVal);
	REAL my_abs(double value);
	REAL sigma1, sigma2;
};


inline REAL Superbee::my_abs(double value) {
	return sqrt(value*value);
}

inline REAL Superbee::minmod(double aVal, double bVal) {
	if(aVal*bVal > 0.) {
		if(my_abs(aVal) < my_abs(bVal)) {
			return aVal;
		} else {
			return bVal;
		}
	} else {
		return 0.;
	}
}

inline REAL Superbee::maxmod(double aVal, double bVal) {
	if(aVal*bVal > 0.) {
		if(my_abs(aVal) > my_abs(bVal)) {
			return aVal;
		} else {
			return bVal;
		}
	} else {
		return 0.;
	}
}

inline REAL Superbee::compute(REAL& up, REAL& uc, REAL& um)
{
	sigma1 = minmod(up,2*um);
	sigma2 = minmod(2*up,um);

	return maxmod(sigma1, sigma2);

}


#endif
