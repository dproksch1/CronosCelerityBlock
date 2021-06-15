#ifndef CRONOS_MHD_BLAST3D_POLAR_H
#define CRONOS_MHD_BLAST3D_POLAR_H = 1

#include "problem.H"

// using namespace std;



class MHDBlast3DPolar: public ProblemType {
public:
	MHDBlast3DPolar(Data &);
	void init_fields(Data &, int [], int []);
private:
	REAL InnerPressure, OuterPressure;
	REAL Blast3D_Radius, Blast3D_Volume;
	REAL Mag0;
	RandomNumberGenerator* random;
};


MHDBlast3DPolar::MHDBlast3DPolar(Data &gdata) : ProblemType(gdata)
{
	name = " 3D MHD blast wave on cyl grid";
#if (GEOM != 2)
	cerr << " Error: Supposed to be a cylindrical problem ";
	cerr << endl;
	exit(2);
#endif
	this->OuterPressure = value((char*)"OuterPressure");
	this->InnerPressure = value((char*)"InnerPressure");
	this->Blast3D_Radius = 3.5*gdata.dx[0];
	this->Blast3D_Radius = value((char*)"Blast_Radius");
	this->Mag0 = value((char*)"Magnetic_Field")/(sqrt(4.*pi));
	this->Mag0 = value((char*)"Magnetic_Field");

}



void MHDBlast3DPolar::init_fields(Data &gdata, int ibeg[3], int iend[3]) 
{

	for(int k=ibeg[2]; k<=iend[2]; ++k) {
		REAL zval = gdata.getCen_z(k);
		for(int j=ibeg[1]; j<=iend[1]; ++j) {
			REAL phi = gdata.getCen_y(j);
			for(int i=ibeg[0]; i<=iend[0]; ++i) {
				REAL r_cyl = gdata.getCen_x(i);
				REAL r_sph = sqrt(sqr(r_cyl) + sqr(zval));
				
				gdata.om[q_rho](i,j,k) = rho0;
				gdata.om[q_sx](i,j,k) = 0.;
				gdata.om[q_sy](i,j,k) = 0.;
				gdata.om[q_sz](i,j,k) = 0.;
				gdata.om[q_Eges](i,j,k) = OuterPressure;

				if(r_sph <= Blast3D_Radius) {
					gdata.om[q_Eges](i,j,k) = InnerPressure;
				}
				gdata.om[q_Eges](i,j,k) *= 1./(gamma-1.);
				
				// Vector potential:
				gdata.om[q_Bx](i,j,k) = 0.;
				gdata.om[q_By](i,j,k) = 0.;
			}
		}
	}

	// z-component of Vector potential -- shift taken into account!
	for(int k=ibeg[2]; k<=iend[2]; ++k) {
		for(int j=ibeg[1]; j<=iend[1]; ++j) {
			REAL phi = gdata.getEdgL_y(j);
			for(int i=ibeg[0]; i<=iend[0]; ++i) {
				REAL r_cyl = gdata.getEdgL_x(i);
				
				gdata.om[q_Bz](i,j,k) = Mag0*r_cyl*sin(phi);

			}
		}
	}


	// cout << " Mag0 " << Mag0 << endl;
	// exit(2);

	gdata.om[q_Bx].rename("A_x");
	gdata.om[q_By].rename("A_y");
	gdata.om[q_Bz].rename("A_z");

}


#endif
