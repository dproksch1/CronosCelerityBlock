#ifndef CRONOS_MOD_TEST_NONLIN_GRID_H
#define CRONOS_MOD_TEST_NONLIN_GRID_H

#include "problem.H"
#include "data.H"
#include "gridgen.H"


class GridFunction_Sin: public Interface_GridFunction {
public:
	GridFunction_Sin() {}
	REAL get_gridFunc(REAL);
};

// function gives 0 at ratio=0 and 1 at ratio=1
REAL GridFunction_Sin::get_gridFunc(REAL ratio) {
	//	return sin(0.5*pi*ratio);
	// return 0.5+0.5*sqr(ratio);
	// return 0.5+0.5*sin(0.375*pi*ratio + 0.0625*pi);
	return sqrt(2)*sin(0.25*pi*ratio);
}

class GridFunction_modEXP: public Interface_GridFunction {
    public:
	GridFunction_modEXP(REAL, REAL, REAL);
  	REAL get_gridFunc(REAL);
  	REAL domain_begin, domain_end, domain_len;
};

GridFunction_modEXP::GridFunction_modEXP(REAL domain_begin,
		REAL domain_end, REAL domain_len) {
	this->domain_begin = domain_begin;
	this->domain_end = domain_end;
	this->domain_len = domain_len;
}


REAL GridFunction_modEXP::get_gridFunc(REAL ratio) {
	REAL a=0.5;
	return (pow(domain_end/domain_begin,a*ratio)-1.0)/(pow(domain_end/domain_begin,a)-1.0);
}

class TestNonlinGrid: public ProblemType {
public:
	TestNonlinGrid(const Data &);
	virtual void init_fields(Data &, int [], int []);
	virtual void writePhystest(Data &, ofstream &);
private:
	REAL amplitude, yama;
	REAL cs2, v0;
};


TestNonlinGrid::TestNonlinGrid(const Data &gdata) : ProblemType(gdata) {
	this->name = " Test nonlinear grid  ";
	this->amplitude = value((char*)"Amplitude");
	this->yama = value((char*)"Adiabatic_exponent");
	this->cs2 = sqr(value((char*)"Isothermal_Soundspeed"));
	this->v0 = value((char*)"Background_Velocity");
}

void TestNonlinGrid::init_fields(Data &gdata, int ibeg[3], int iend[3])
{
	//REAL P0
	for (int k = ibeg[2]; k <= iend[2]; k++){
		for (int j = ibeg[1]; j <= iend[1]; j++){
			for (int i = ibeg[0]; i <= iend[0]; i++){
				REAL xx = gdata.getCen_x(i);
	
				gdata.om[q_rho](i,j,k) = 1. + amplitude*sin(2.*pi*xx);
				gdata.om[q_sx](i,j,k) = v0 - amplitude*sin(2.*pi*xx);

				// if(k==0 && j==0 && i>=0) {
				// 	cout << i << " " << gdata.getCen_x(i) << " ";
				// 	// cout << gdata.getEdgL_x(i) << " ";
				// 	// cout << gdata.getEdgR_x(i);
				// 	// cout << endl;
				// 	// cout << gdata.getCen_dx(0,i)<< endl;
				// 	cout << gdata.om[q_rho](i,j,k) << endl;
				// }
			}
		}
	}
}


void TestNonlinGrid::writePhystest(Data &gdata, ofstream &outfile)
{
	REAL rhoVal = gdata.om[q_rho](gdata.mx[0]/4,0,0);
	if(gdata.rank == 0) {
		outfile << " " << rhoVal;
	}
}



#endif
