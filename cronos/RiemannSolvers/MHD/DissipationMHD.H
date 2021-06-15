#ifndef CRONOS_DISSIPATIONMHD_H
#define CRONOS_DISSIPATIONMHD_H = 1

#include "Dissipation.H"

class DissipationMHD: public Dissipation {
public:
	DissipationMHD(){}
	void modifyEmf(Data &, ProblemType &, fields_2D &, int, int);
};

#endif
