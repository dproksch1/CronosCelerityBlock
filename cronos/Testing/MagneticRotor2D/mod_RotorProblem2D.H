#ifndef CRONOS_ROTOR_PROB2D_H
#define CRONOS_ROTOR_PROB2D_H = 1

#include "problem.H"

// using namespace std;

//#if (GEOM < 3)
//#define GEOM 1
class RotorProblem2D: public ProblemType {
public:
	RotorProblem2D(Data &);
	void init_fields(Data &, int [], int []);
private:
	REAL r_inner, r_outer;
	REAL Mag0;
};


RotorProblem2D::RotorProblem2D(Data &gdata) : ProblemType(gdata)
{
#if (FLUID_TYPE == CRONOS_MHD)
	name = " 2D magnetic rotator problem ";
#elif (FLUID_TYPE == CRONOS_HYDRO)
	name = " 2D hydro rotator problem ";
#endif
#if (GEOM == 1)
	name += "(Cartesian) ";
#elif (GEOM == 2)
	name += "(Cylindrical) ";
#endif
	this->r_inner = value((char*)"InnerRadius");
	this->r_outer = value((char*)"OuterRadius");
	this->Mag0 = 5./(sqrt(4.*pi));
	// this->Mag0 = 1.;

}



void RotorProblem2D::init_fields(Data &gdata, int ibeg[3], int iend[3]) 
{
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
				
				REAL func_rad(0.);
				if(r_cyl < r_inner) {
					func_rad = 1.;
				} else if (r_cyl < r_outer) {
					func_rad = (r_outer - r_cyl)/r_outer;
				} else {
					func_rad = 0.;
				}

				// rho:
				gdata.om[q_rho](i,j,k) = rho0*(1. + 9.*func_rad);
				// v_R
				gdata.om[q_sx](i,j,k) = 0.;
				// v_phi
				REAL Omega = 2.*func_rad/r_inner;
				// gdata.om[q_sy](i,j,k) = 2.*func_rad;
				REAL v_phi = r_cyl*Omega;
#if (GEOM == 1)
				gdata.om[q_sx](i,j,k) = -v_phi*sin(phi);
				gdata.om[q_sy](i,j,k) =  v_phi*cos(phi);
#elif (GEOM == 2)
				gdata.om[q_sy](i,j,k) = v_phi;
#endif
				gdata.om[q_sz](i,j,k) = 0.;
				gdata.om[q_Eges](i,j,k) = 1./(gamma - 1.);

#if (FLUID_TYPE == CRONOS_MHD)
				// Vector potential:
				gdata.om[q_Bx](i,j,k) = 0.;
				gdata.om[q_By](i,j,k) = 0.;
#endif
			}
		}
	}
	// exit(2);


#if (FLUID_TYPE == CRONOS_MHD)

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


	gdata.om[q_Bx].rename("A_x");
	gdata.om[q_By].rename("A_y");
	gdata.om[q_Bz].rename("A_z");

#endif

}

#endif

//#endif
