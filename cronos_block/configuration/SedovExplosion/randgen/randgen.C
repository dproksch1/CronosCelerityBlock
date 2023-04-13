#include "randgen.H"
#include<ctime>

RandomNumberGenerator::RandomNumberGenerator(): RandomNumberGenerator(std::time(0)) {}


RandomNumberGenerator::RandomNumberGenerator(const long int & seed) {
	this->seed = seed;
	this->count = 0;
	this->gset = 0.;
	this->iset = false;
	RandGen = new MTRand(seed);
	for(int i=0; i<1500; ++i) {
		getRand();
	}
}

double RandomNumberGenerator::getRand() {
	count++;
	return RandGen->rand();
}

double RandomNumberGenerator::getRand(const double &max) {
	return RandGen->rand(max);
}



double RandomNumberGenerator::getRandGauss(const double &sigma) {
	double fac;
  
	if(!iset){
		double radsq(0.);
		double v1(0.), v2(0.);
		while(radsq >= 1. || radsq == 0.){
			v1 = getRand(2.)-1.;
			v2 = getRand(2.)-1.;
			radsq = sqr(v1)+sqr(v2);}
		fac = sqrt(-2.*sqr(sigma)*log(radsq)/radsq);
		gset = v1*fac;
		iset = 1;
		return v2*fac;
	} else{
		iset = false;
		return gset;
	}
}


void RandomNumberGenerator::flushGauss() {
	iset = false;
}


void RandomNumberGenerator::getRandState(MTRand::uint32 RandState[MTRand::SAVE]) 
{
	RandGen->save(RandState);
}


void RandomNumberGenerator::setRandState(MTRand::uint32 RandState[MTRand::SAVE]) 
{
	RandGen->load(RandState);
}

long int RandomNumberGenerator::getCount() {
	return count;
}
