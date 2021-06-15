#ifndef CRONOS_DISSIPATION_H
#define CRONOS_DISSIPATION_H = 1

#include "data.H"
#include "problem.H"
#include "fields_local.H"

class Dissipation {
public:
	Dissipation(){};
	void modifyFlux(Data &, ProblemType &, fields_1D &, const int &dir);
};

#endif
