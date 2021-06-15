#ifndef CRONOS_MOD_MAGLOOPCART_H
#define CRONOS_MOD_MAGLOOPCARTH = 1

#include "problem.H"
#include "data.H"
#include "gridgen.H"


class MagLoopCart: public ProblemType {
public:
	MagLoopCart(const Data &);
	virtual void init_fields(Data &, int [], int []);
private:

	REAL rho_init, P_init, v0;
	REAL A0, radA;
	REAL cos, sin;
	int dir;
	bool isothermal;
};


MagLoopCart::MagLoopCart(const Data &gdata) : ProblemType(gdata) {

#if (FLUID_TYPE != CRONOS_MHD)
	cerr << endl;
	cerr << "!---------------------------------------------!"<< endl;
	cerr << "! ERROR: Need to use FLUID_TYPE CRONOS_MHD    !"<< endl;
	cerr << "!    --> correct in constants.H! Exiting...   !"<< endl;
	cerr << "!---------------------------------------------!"<< endl;
	cerr << endl;
	exit(2);
#endif

	if(ENERGETICS == FULL) {
		name = " magnetic loop advection (2D) ";
		cerr << " Has to be isothermal " << endl;
		exit(3);
		isothermal = false;
	} else {
		name = " magnetic loop adection - isothermal (2D) ";
		if(gamma > 1.) {
			cerr << " Has to be isothermal ";
			exit(2);
		}
		isothermal = true;
	}


	rho_init  = value((char*)"InitialDensity");
	P_init= value((char*)"InitialPressure");
	A0 = value((char*)"VectorPotential");
	radA = value((char*)"radVecPot");
	v0 = value((char*)"AdvectionSpeed");
	dir = value((char*)"direction");
	v0 *= sqrt(5.);
	cos = 2./sqrt(5.);
	sin = 1./sqrt(5.);

}

void MagLoopCart::init_fields(Data &gdata, int ibeg[3], int iend[3])
{

	for (int k = ibeg[2]; k <= iend[2]; k++){
		REAL zz = gdata.getCen_z(k);
		for (int j = ibeg[1]; j <= iend[1]; j++){
			REAL yy = gdata.getCen_y(j);
			for (int i = ibeg[0]; i <= iend[0]; i++){
				REAL xx = gdata.getCen_x(i);

				REAL r_sph;
				if(dir==0) {
					r_sph = sqrt(xx*xx + yy*yy);
				} else if (dir==1) {
					r_sph = sqrt(xx*xx + zz*zz);
				} else {
					r_sph = sqrt(yy*yy + zz*zz);
				}

				gdata.om[q_rho](i,j,k) = rho_init;
				if(dir==0) {
					gdata.om[q_sx](i,j,k) = v0*cos;
					gdata.om[q_sy](i,j,k) = v0*sin;
					gdata.om[q_sz](i,j,k) = 0.;
				} else if (dir==1) {
					gdata.om[q_sx](i,j,k) = v0*cos;
					gdata.om[q_sy](i,j,k) = 0.;
					gdata.om[q_sz](i,j,k) = v0*sin;
				} else {
					gdata.om[q_sx](i,j,k) = 0.;
					gdata.om[q_sy](i,j,k) = v0*cos;
					gdata.om[q_sz](i,j,k) = v0*sin;
				}

				gdata.om[q_sx](i,j,k) = 0;
				gdata.om[q_sy](i,j,k) = v0;
				gdata.om[q_sz](i,j,k) = v0;


				gdata.om[q_Bx](i,j,k) = 0.;
				gdata.om[q_By](i,j,k) = 0.;
				gdata.om[q_Bz](i,j,k) = 0.;

				if(r_sph <= radA) {
					if(dir==0) {
						gdata.om[q_Bz](i,j,k) = A0*(radA - r_sph);
					} else if (dir==1) {
						gdata.om[q_By](i,j,k) = A0*(radA - r_sph);
					} else {
						gdata.om[q_Bx](i,j,k) = A0*(radA - r_sph);
					}

				}

			}
		}
	}
	gdata.om[q_Bx].rename("A_x");
	gdata.om[q_By].rename("A_y");
	gdata.om[q_Bz].rename("A_z");
}





#endif
