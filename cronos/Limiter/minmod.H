#ifndef CRONOS_MINMOD_H
#define CRONOS_MINMOD_H = 1
#include <iostream>

#include "constants.H"

using namespace std;

class minmod /*: public limiter_func*/ {
public:
	minmod() {
		Theta = value((char*)"Theta");
	}

	minmod(const minmod&) = default;

	REAL compute(REAL&, REAL&, REAL&) const;
private:
	REAL Theta;
};

inline REAL minmod::compute(REAL& up, REAL& uc, REAL& um) const
{
	up *= Theta;
	um *= Theta;
	if(up*um > 0.) {
		if(up > 0){
			return fmin(up,(fmin(uc,um)));
		} else if(up < 0){
			return fmax(up,(fmax(uc,um)));
		} else {
			return 0.0;
			//cerr << " Error " << endl;
			//exit(1);
		}
	} else {
		return 0.;
	}
}


#endif
