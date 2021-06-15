#ifndef CRONOS_MOD_BRIOWU2D_H
#define CRONOS_MOD_BRIOWU2D_H = 1

#include "problem.H"
#include "data.H"
#include "gridgen.H"


class BrioWuTest2D: public ProblemType {
public:
	BrioWuTest2D(const Data &);
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
	REAL alpha;
	int shockdir;
	bool isothermal;
};


BrioWuTest2D::BrioWuTest2D(const Data &gdata) : ProblemType(gdata) {

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

	BLParInit  = value((char*)"B0Lpar");
	BLPerpInit = value((char*)"B0Lperp");
	BRParInit  = value((char*)"B0Rpar");
	BRPerpInit = value((char*)"B0Rperp");

	PLInit  = value((char*)"PL0");
	PRInit  = value((char*)"PR0");

//	shock_angle = value((char*)"shock_angle");
	alpha = M_PI*0.25;

	shockdir = static_cast<int>(value((char*)"shockdir"));


	BLXInit = cos(alpha)*BLParInit - sin(alpha)*BLPerpInit;
	BLYInit = sin(alpha)*BLParInit - cos(alpha)*BLPerpInit;

	BLXInit = 0.;
	BLYInit = sqrt(2);


	BRXInit = cos(alpha)*BRParInit - sin(alpha)*BRPerpInit;
	BRYInit = sin(alpha)*BRParInit - cos(alpha)*BRPerpInit;
	BRXInit = sqrt(2)*0.5;
	BRYInit = sqrt(2)*0.5;
	BLZInit = 0.;
	BRZInit = 0.;



}

void BrioWuTest2D::init_fields(Data &gdata, int ibeg[3], int iend[3])
{
	REAL xPar;
	// gdata.om[q_Eges].rename("Temp");

	// do a rotation around the z-axis

	for (int k = ibeg[2]; k <= iend[2]; k++){
		cout << k << endl;
		REAL zz = gdata.getCen_z(k);
		for (int j = ibeg[1]; j <= iend[1]; j++){
			REAL yy = gdata.getCen_y(j);
			REAL yR = gdata.getEdgR_y(j);
			for (int i = ibeg[0]; i <= iend[0]; i++){
				REAL xx = gdata.getCen_x(i);
				REAL xR = gdata.getEdgR_x(i);

				xPar = xx + yy;

//				if(xPar < 0.5 || xPar >= 1.5){
				if((xPar > 0 && xPar < 0.5) || (xPar >= 1. && xPar <1.5)){
//				if(xPar < 1){
					gdata.om[q_rho](i,j,k) = rhoLInit;
					gdata.om[q_sx](i,j,k) = cos(alpha)*VLXInit - sin(alpha)*VLYInit;
					gdata.om[q_sy](i,j,k) = sin(alpha)*VLXInit + cos(alpha)*VLYInit;
					gdata.om[q_sz](i,j,k) = VLZInit;
					// 					gdata.om[q_Bx](i,j,k) = BLXInit;
					// 					gdata.om[q_By](i,j,k) = BLYInit;
					// 					gdata.om[q_Bz](i,j,k) = BLZInit;
					gdata.om[q_Bx](i,j,k) = BLYInit*zz;
					gdata.om[q_By](i,j,k) = BLZInit*xx;
					gdata.om[q_Bz](i,j,k) = BLXInit*yy;
					gdata.om[q_Bx](i,j,k) = BLXInit;
					gdata.om[q_By](i,j,k) = BLYInit;
					gdata.om[q_Bz](i,j,k) = 0.;

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
					gdata.om[q_Bx](i,j,k) = BRXInit;
					gdata.om[q_By](i,j,k) = BRYInit;
					gdata.om[q_Bz](i,j,k) = 0.;

					if(!isothermal) {
						gdata.om[q_Eges](i,j,k) = PRInit;
						if(gdata.om[q_Eges].getName() == "Temp") {
							gdata.om[q_Eges](i,j,k) /= gdata.om[q_rho](i,j,k);
						} else {
							 gdata.om[q_Eges](i,j,k) *= 1./(gamma-1.);
						}
					}
				}

				// Special treatment for magnetic field
				// Bx
				xPar = xR + yy;
//				if(xPar < 1){
//				if(xPar < 0.5 || xPar >= 1.5){
				if((xPar > 0 && xPar < 0.5) || (xPar >= 1. && xPar <1.5)){
					gdata.om[q_Bx](i,j,k) = BLXInit;
				} else {
					gdata.om[q_Bx](i,j,k) = BRXInit;
				}

				// By
				xPar = xx + yR;
//				if(xPar < 1.){
//				if(xPar < 0.5 || xPar >= 1.5){
				if((xPar > 0 && xPar < 0.5) || (xPar >= 1. && xPar <1.5)){
					gdata.om[q_By](i,j,k) = BLYInit;
				} else {
					gdata.om[q_By](i,j,k) = BRYInit;
				}

				gdata.om[q_Bz](i,j,k) = 0.;

#if (CRSWITCH_DUAL_ENERGY == CRONOS_ON)
				gdata.om[q_Eadd](i,j,k) = gdata.om[q_Eges](i,j,k)*pow(gdata.om[q_rho](i,j,k),1-gamma);
#endif
	
			}
		}
	}

//	// Hack: repair y-magnetic field:
//	for (int k = ibeg[2]; k <= iend[2]; k++){
//		for (int j = ibeg[1]+1; j <= iend[1]; j++){
//			for (int i = ibeg[0]+1; i <= iend[0]; i++){
//				gdata.om[q_By](i,j,k) = gdata.om[q_By](i,j-1,k) -
//						gdata.dx[1]/gdata.dx[0]*(gdata.om[q_Bx](i,j,k) - gdata.om[q_Bx](i-1,j,k));
//			}
//		}
//	}

//	gdata.om[q_Bx].rename("A_x");
//	gdata.om[q_By].rename("A_y");
//	gdata.om[q_Bz].rename("A_z");
}





#endif
