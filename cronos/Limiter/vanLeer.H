#ifndef CRONOS_VANLEER_H
#define CRONOS_VANLEER_H = 1

#include "constants.H"
#include <algorithm>

class vanLeer /*: public limiter_func*/ {
public:
	vanLeer() {
	}
	REAL compute(REAL&, REAL&, REAL&);
};

inline REAL vanLeer::compute(REAL& up, REAL& uc, REAL& um)
{
	REAL prod = up*um;
	if (prod > 0.) {
		return prod / uc;
	} else {
		return 0.;
	}
}


#endif
