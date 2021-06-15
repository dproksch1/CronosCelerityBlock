#ifndef CRONOS_MOD_BRIOWU_H
#define CRONOS_MOD_BRIOWU_H = 1

#include "problem.H"
#include "data.H"
#include "gridgen.H"


class BrioWuTest: public ProblemType {
public:
	BrioWuTest(const Data &);
	virtual void init_fields(Data &, int [], int []);
private:
	REAL rhoLInit, rhoRInit;
	REAL VLXInit, VLYInit, VLZInit;
	REAL VRXInit, VRYInit, VRZInit;
	REAL BLParInit, BLPerpInit;
	REAL BRParInit, BRPerpInit;
	REAL BLXInit, BLYInit, BLZInit;
	REAL BRXInit, BRYInit, BRZInit;
	REAL PLInit, PRInit;
	int shockdir;
	bool isothermal;
};


BrioWuTest::BrioWuTest(const Data &gdata) : ProblemType(gdata) {

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
		name = " Shock Tube problem (MHD) 1D ";
		isothermal = false;
	} else {
		name = " Shock Tube problem (MHD) 1D -- isothermal ";
		if(gamma > 1.) {
			cerr << " Has to be isothermal ";
			exit(2);
		}
		isothermal = true;
	}

	rhoLInit = value((char*)"rho0L");
	rhoRInit = value((char*)"rho0R");

	VLXInit = value((char*)"v0xL");
	VLYInit = value((char*)"v0yL");
	VLZInit = value((char*)"v0zL");
	VRXInit = value((char*)"v0xR");
	VRYInit = value((char*)"v0yR");
	VRZInit = value((char*)"v0zR");


	if(value_exists("B0xL")) {
		BLXInit = value((char*)"B0xL");
		BLYInit = value((char*)"B0yL");
		BLZInit = value((char*)"B0zL");
		BRXInit = value((char*)"B0xR");
		BRYInit = value((char*)"B0yR");
		BRZInit = value((char*)"B0zR");
	}

	PLInit  = value((char*)"PL0");
	PRInit  = value((char*)"PR0");

	shockdir = static_cast<int>(value((char*)"shockdir"));

	if(value_exists("B0Lpar")) {
		BLParInit  = value((char*)"B0Lpar");
		BLPerpInit = value((char*)"B0Lperp");
		BRParInit  = value((char*)"B0Rpar");
		BRPerpInit = value((char*)"B0Rperp");
		if(shockdir == 1) {
			BLXInit = BLParInit;
			BRXInit = BRParInit;
			BLYInit = BLPerpInit;
			BRYInit = BRPerpInit;
			BLZInit = 0.;
			BRZInit = 0.;
		} else if (shockdir == 2) {
			BLXInit = 0.;
			BRXInit = 0.;
			BLYInit = BLParInit;
			BRYInit = BRParInit;
			BLZInit = BLPerpInit;
			BRZInit = BRPerpInit;
		} else {
			BLZInit = BLParInit;
			BRZInit = BRParInit;
			BLXInit = BLPerpInit;
			BRXInit = BRPerpInit;
			BLYInit = 0.;
			BRYInit = 0.;
		}
	}

	if(value_exists("mod_Mag_CGS")) {
		BLXInit *= 1./sqrt(4.*M_PI);
		BLYInit *= 1./sqrt(4.*M_PI);
		BLZInit *= 1./sqrt(4.*M_PI);
		BRXInit *= 1./sqrt(4.*M_PI);
		BRYInit *= 1./sqrt(4.*M_PI);
		BRZInit *= 1./sqrt(4.*M_PI);
	}
}

void BrioWuTest::init_fields(Data &gdata, int ibeg[3], int iend[3])
{
	REAL xPar;
	// gdata.om[q_Eges].rename("Temp");

	for (int k = ibeg[2]; k <= iend[2]; k++){
		cout << k << endl;
		REAL zz = gdata.getCen_z(k);
		for (int j = ibeg[1]; j <= iend[1]; j++){
			REAL yy = gdata.getCen_y(j);
			for (int i = ibeg[0]; i <= iend[0]; i++){
				REAL xx = gdata.getCen_x(i);

				if(shockdir == 1){
					xPar = xx;
				} else if (shockdir == 2){
					xPar= yy;
				} else {
					xPar = zz;
				}

				if(xPar < 0.5){
					gdata.om[q_rho](i,j,k) = rhoLInit;
					gdata.om[q_sx](i,j,k) = VLXInit;
					gdata.om[q_sy](i,j,k) = VLYInit;
					gdata.om[q_sz](i,j,k) = VLZInit;
					// 					gdata.om[q_Bx](i,j,k) = BLXInit;
					// 					gdata.om[q_By](i,j,k) = BLYInit;
					// 					gdata.om[q_Bz](i,j,k) = BLZInit;
					gdata.om[q_Bx](i,j,k) = BLYInit*zz;
					gdata.om[q_By](i,j,k) = BLZInit*xx;
					gdata.om[q_Bz](i,j,k) = BLXInit*yy;

					if(!isothermal) {
						gdata.om[q_Eges](i,j,k) = PLInit;
						if(gdata.om[q_Eges].getName() == "Temp") {
							gdata.om[q_Eges](i,j,k) /= gdata.om[q_rho](i,j,k);
						} else {
							gdata.om[q_Eges](i,j,k) *= 1./(gamma-1.);
						}
					}
				} else {
					gdata.om[q_rho](i,j,k) = rhoRInit;
					gdata.om[q_sx](i,j,k) = VRXInit;
					gdata.om[q_sy](i,j,k) = VRYInit;
					gdata.om[q_sz](i,j,k) = VRZInit;
					// 					gdata.om[q_Bx](i,j,k) = BRXInit;
					// 					gdata.om[q_By](i,j,k) = BRYInit;
					// 					gdata.om[q_Bz](i,j,k) = BRZInit;
					gdata.om[q_Bx](i,j,k) = BRYInit*zz;
					gdata.om[q_By](i,j,k) = BRZInit*xx;
					gdata.om[q_Bz](i,j,k) = BRXInit*yy;

					if(!isothermal) {
						gdata.om[q_Eges](i,j,k) = PRInit;
						if(gdata.om[q_Eges].getName() == "Temp") {
							gdata.om[q_Eges](i,j,k) /= gdata.om[q_rho](i,j,k);
						} else {
							 gdata.om[q_Eges](i,j,k) *= 1./(gamma-1.);
						}
					}
				}
#if (CRSWITCH_DUAL_ENERGY == CRONOS_ON)
				gdata.om[q_Eadd](i,j,k) = gdata.om[q_Eges](i,j,k)*pow(gdata.om[q_rho](i,j,k),1-gamma);
#endif
	
			}
		}
	}
	gdata.om[q_Bx].rename("A_x");
	gdata.om[q_By].rename("A_y");
	gdata.om[q_Bz].rename("A_z");
}





#endif
