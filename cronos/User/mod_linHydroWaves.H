#ifndef CRONOS_MOD_LINHYDROWAVES_H
#define CRONOS_MOD_LINHYDROWAVES_H = 1

#include "problem.H"
#include "data.H"
#include "gridgen.H"


class linHydroWaves: public ProblemType {
public:
	linHydroWaves(const Data &);
	virtual void init_fields(Data &, int [], int []);
	virtual void writePhystest(Data &, ofstream &);
private:
	REAL amplitude, yama;
	REAL cs2, v0;
};


linHydroWaves::linHydroWaves(const Data &gdata) : ProblemType(gdata) {
	this->name = "Linear Hydro Wave Test ";
	this->amplitude = value((char*)"Amplitude");
	this->yama = value((char*)"Adiabatic_exponent");
	this->cs2 = sqr(value((char*)"Isothermal_Soundspeed"));
	this->v0 = value((char*)"Background_Velocity");
}

void linHydroWaves::init_fields(Data &gdata, int ibeg[3], int iend[3])
{
	//REAL P0
	for (int k = ibeg[2]; k <= iend[2]; k++){
		for (int j = ibeg[1]; j <= iend[1]; j++){
			for (int i = ibeg[0]; i <= iend[0]; i++){
				REAL xx = gdata.getCen_x(i);
	
				gdata.om[q_rho](i,j,k) = 1. + amplitude*sin(2.*pi*xx);
				gdata.om[q_sx](i,j,k) = v0 - amplitude*sin(2.*pi*xx);

	
			}
		}
	}
}


void linHydroWaves::writePhystest(Data &gdata, ofstream &outfile)
{
	REAL rhoVal = gdata.om[q_rho](gdata.mx[0]/4,0,0);
	if(gdata.rank == 0) {
		outfile << " " << rhoVal;
	}
}



#endif
