#ifndef CRONOS_MOD_SOD_H
#define CRONOS_MOD_SOD_H = 1

#include "problem.H"
#include "data.H"
#include "gridgen.H"


class SOD: public ProblemType {
public:
	SOD(const Data &);
	virtual void init_fields(Data &, int[], int[]);
private:
	REAL gamma;
	REAL rhoL, rhoR;
	REAL uL, uR;
	REAL pL, pR;
	REAL position;
	int shockDir;
};


SOD::SOD(const Data &gdata) : ProblemType(gdata) {


	gamma = value((char*)"Adiabatic_exponent");
	rhoL = value((char*)"rho_left"); 
	rhoR = value((char*)"rho_right");
	uL = value((char*)"u_left");
	uR = value((char*)"u_right");
	pL = value((char*)"p_left");
	pR =  value((char*)"p_right");
	position = value((char*)"Shock_Position");
	shockDir = static_cast<int>(value((char*)"Shock_Direction"));

}


void SOD::init_fields(Data &gdata, int ibeg[0], int iend[0]) {

	for(int k=ibeg[2]; k<=iend[2]; ++k) {
		for(int j=ibeg[1]; j<=iend[1]; ++j) {
			for(int i=ibeg[0]; i<=iend[0]; ++i) {
				REAL xx = gdata.getCen_x(i);
				if(shockDir == 1) {
					xx = gdata.getCen_y(j);
				} else if (shockDir == 2) {
					xx = gdata.getCen_z(k);
				}


				// Velocity
				gdata.om[q_sx](i,j,k) = 0.;
				gdata.om[q_sy](i,j,k) = 0.;
				gdata.om[q_sz](i,j,k) = 0.;


				if (xx <= position) {
					// Density
					gdata.om[q_rho](i,j,k) = rhoL;
					gdata.om[shockDir+1](i,j,k) = uL;
				} else {
					gdata.om[q_rho](i,j,k) = rhoR;
					gdata.om[shockDir+1](i,j,k) = uL;
				}
				


		
				// Thermal energy
				if (xx <= position) {
					gdata.om[q_Eges](i,j,k) = pL/(gamma-1);
				} else {
					gdata.om[q_Eges](i,j,k) = pR/(gamma-1);
				}

				if(i==0 && j==0) {
					cout << k << " " << gdata.om[0](i,j,k) << " ";
					cout << gdata.om[4](i,j,k) << endl;
				}
								
			}
		}
	}
}



#endif











