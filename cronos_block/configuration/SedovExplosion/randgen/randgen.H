#ifndef CRONOS_RANDGEN_H
#define CRONOS_RANDGEN_H = 1
#include <stdlib.h>
#include <math.h>
#include "util.H"
#include "MersenneTwister.h"

class RandomNumberGenerator{
public:
	RandomNumberGenerator();
	RandomNumberGenerator(const long int &);
	double getRand();
	double getRand(const double &);
	double getRandGauss(const double &);
	long int getCount();
	void getRandState(MTRand::uint32 [MTRand::SAVE]);
	void setRandState(MTRand::uint32 [MTRand::SAVE]);
	void flushGauss();
private:
	long int seed;
	bool iset;
	double gset;
	MTRand *RandGen;
	long int count;
};


#endif
