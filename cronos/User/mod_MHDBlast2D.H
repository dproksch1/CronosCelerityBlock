#ifndef CRONOS_MHD_BLAST_2D_H
#define CRONOS_MHD_BLAST_2D_H = 1

#include "problem.H"

// using namespace std;


class MHDBlast2D: public ProblemType {
public:
	MHDBlast2D(Data &);
	void init_fields(Data &, int [], int []);
private:
	REAL InnerPressure, OuterPressure;
	REAL Blast_Radius, Blast_Volume;
	REAL Mag0;
	RandomNumberGenerator* random;
};


MHDBlast2D::MHDBlast2D(Data &gdata) : ProblemType(gdata)
{
	name = " 2D MHD blast wave ";
#if (GEOM == 1)
	name += "(Cartesian) ";
#elif (GEOM == 2)
	name += "(Cylindrical) ";
#endif
#if (GEOM > 2)
	cerr << " Error: Supposed to be a cylindrical or Cartesian problem ";
	cerr << endl;
	exit(2);
#endif
	this->OuterPressure = value((char*)"OuterPressure");
	this->InnerPressure = value((char*)"InnerPressure");
	this->Blast_Radius = 3.5*gdata.dx[0];
	this->Blast_Radius = value((char*)"Blast_Radius");
	this->Mag0 = value((char*)"Magnetic_Field")/(sqrt(4.*pi));
	this->Mag0 = value((char*)"Magnetic_Field");

}



void MHDBlast2D::init_fields(Data &gdata, int ibeg[3], int iend[3]) 
{

#if (GEOM < 3)

	for(int k=ibeg[2]; k<=iend[2]; ++k) {
		for(int j=ibeg[1]; j<=iend[1]; ++j) {
#if (GEOM == 1)
			REAL yval = gdata.getCen_y(j);
#elif (GEOM == 2)
			REAL phi = gdata.getCen_y(j);
#endif
			for(int i=ibeg[0]; i<=iend[0]; ++i) {
#if (GEOM == 1)
				REAL xval = gdata.getCen_x(i);
				REAL r_cyl = sqrt(sqr(xval) + sqr(yval));
				REAL phi = atan2(yval, xval);
#elif (GEOM == 2)
				REAL r_cyl = gdata.getCen_x(i);
#endif
				
				gdata.om[q_rho](i,j,k) = rho0;
				gdata.om[q_sx](i,j,k) = 0.;
				gdata.om[q_sy](i,j,k) = 0.;
				gdata.om[q_sz](i,j,k) = 0.;
				gdata.om[q_Eges](i,j,k) = OuterPressure;

				if(r_cyl <= Blast_Radius) {
					gdata.om[q_Eges](i,j,k) = InnerPressure;
				}
#if (ENERGETICS == FULL)
				gdata.om[q_Eges](i,j,k) *= 1./(gamma-1.);
#endif
				
				// Vector potential:
				gdata.om[q_Bx](i,j,k) = 0.;
				gdata.om[q_By](i,j,k) = 0.;
			}
		}
	}

	// z-component of Vector potential -- shift taken into account!
	for(int k=ibeg[2]; k<=iend[2]; ++k) {
		for(int j=ibeg[1]; j<=iend[1]; ++j) {
			for(int i=ibeg[0]; i<=iend[0]; ++i) {
#if (GEOM == 1)
				REAL yval = gdata.getEdgL_y(j);

#elif (GEOM == 2)
				REAL phi = gdata.getEdgL_y(j);
				REAL r_cyl = gdata.getEdgL_x(i);
				REAL yval = r_cyl*sin(phi);
#endif
				gdata.om[q_Bz](i,j,k) = Mag0*yval;

			}
		}
	}


	cout << " Mag0 " << Mag0 << endl;
	// exit(2);

	gdata.om[q_Bx].rename("A_x");
	gdata.om[q_By].rename("A_y");
	gdata.om[q_Bz].rename("A_z");

#endif

}


#endif
