#ifndef CRONOS_MOD_MULT_TEST_H
#define CRONOS_MOD_MULT_TEST_H = 1

#include "problem.H"
#include "data.H"
#include "gridgen.H"


class MultTest: public ProblemType {
public:
	MultTest(const Data &);
	virtual void init_fields(Data &, int [], int []);
private:
	REAL rhoLInit_MHD, rhoRInit_MHD;
	REAL VLXInit_MHD, VLYInit_MHD, VLZInit_MHD;
	REAL VRXInit_MHD, VRYInit_MHD, VRZInit_MHD;
	REAL BLParInit, BLPerpInit;
	REAL BRParInit, BRPerpInit;
	REAL BLXInit, BLYInit, BLZInit;
	REAL BRXInit, BRYInit, BRZInit;
	REAL PLInit_MHD, PRInit_MHD;
	REAL rhoLInit_HD, rhoRInit_HD;
	REAL VLXInit_HD, VLYInit_HD, VLZInit_HD;
	REAL VRXInit_HD, VRYInit_HD, VRZInit_HD;
	REAL PLInit_HD, PRInit_HD;
	int shockdir;
	bool isothermal;
};


MultTest::MultTest(const Data &gdata) : ProblemType(gdata) {

#if (FLUID_TYPE != CRONOS_MULTIFLUID)
	cerr << endl;
	cerr << "!---------------------------------------------!"<< endl;
	cerr << "! ERROR: Need to use FLUID_TYPE CRONOS_MULTIFLUID    !"<< endl;
	cerr << "!    --> correct in constants.H! Exiting...   !"<< endl;
	cerr << "!---------------------------------------------!"<< endl;
	cerr << endl;
	exit(2);
#endif

	if(ENERGETICS == FULL) {
		name = " Shock Tube problem (MHD & HD) 1D ";
		isothermal = false;
	} else {
		name = " Shock Tube problem (MHD & HD) 1D -- isothermal ";
		if(gamma > 1.) {
			cerr << " Has to be isothermal ";
			exit(2);
		}
		isothermal = true;
	}

	rhoLInit_MHD = value((char*)"rho0L_MHD");
	rhoRInit_MHD = value((char*)"rho0R_MHD");

	VLXInit_MHD = value((char*)"v0xL_MHD");
	VLYInit_MHD = value((char*)"v0yL_MHD");
	VLZInit_MHD = value((char*)"v0zL_MHD");
	VRXInit_MHD = value((char*)"v0xR_MHD");
	VRYInit_MHD = value((char*)"v0yR_MHD");
	VRZInit_MHD = value((char*)"v0zR_MHD");

	BLParInit  = value((char*)"B0Lpar");
	BLPerpInit = value((char*)"B0Lperp");
	BRParInit  = value((char*)"B0Rpar");
	BRPerpInit = value((char*)"B0Rperp");

#if(ENERGETICS == FULL)
	PLInit_MHD  = value((char*)"PL0_MHD");
	PRInit_MHD  = value((char*)"PR0_MHD");
#endif

	rhoLInit_HD = value((char*)"rho0L_HD");
	rhoRInit_HD = value((char*)"rho0R_HD");


	int u0L_HD = value((char*)"u0L_HD");
	int u0R_HD = value((char*)"u0R_HD");


#if(ENERGETICS == FULL)
	PLInit_HD  = value((char*)"p0L_HD");
	PRInit_HD  = value((char*)"p0R_HD");
#endif

	shockdir = static_cast<int>(value((char*)"shockdir"));

	if(shockdir == 1) {
		BLXInit = BLParInit;
		BRXInit = BRParInit;
		BLYInit = BLPerpInit;
		BRYInit = BRPerpInit;
		BLZInit = 0.;
		BRZInit = 0.;

		VLXInit_HD = u0L_HD;
		VLYInit_HD = u0R_HD;
		VLZInit_HD = 0;
		VRXInit_HD = 0;
		VRYInit_HD = 0;
		VRZInit_HD = 0;

	} else if (shockdir == 2) {
		BLXInit = 0.;
		BRXInit = 0.;
		BLYInit = BLParInit;
		BRYInit = BRParInit;
		BLZInit = BLPerpInit;
		BRZInit = BRPerpInit;

		VLXInit_HD = 0;
		VLYInit_HD = 0;
		VLZInit_HD = u0L_HD;
		VRXInit_HD = u0R_HD;
		VRYInit_HD = 0;
		VRZInit_HD = 0;

	} else {
		BLZInit = BLParInit;
		BRZInit = BRParInit;
		BLXInit = BLPerpInit;
		BRXInit = BRPerpInit;
		BLYInit = 0.;
		BRYInit = 0.;

		VLXInit_HD = 0;
		VLYInit_HD = 0;
		VLZInit_HD = 0;
		VRXInit_HD = 0;
		VRYInit_HD = u0L_HD;
		VRZInit_HD = u0R_HD;
	}
}

void MultTest::init_fields(Data &gdata, int ibeg[3], int iend[3])
{
	REAL xPar;
	// gdata.om[q_Eges].rename("Temp");

	int numFluids = gdata.fluids->get_numFluids();
	q_Bx = gdata.fluids->get_q_Bx();
	q_By = gdata.fluids->get_q_By();
	q_Bz = gdata.fluids->get_q_Bz();

	for(int iFluid=0; iFluid<numFluids; ++iFluid) {

		int fluidType = gdata.fluids->get_fluidType(iFluid);

		cout << " Init for fluid " << iFluid << " of type " << fluidType << endl;

		q_rho = gdata.fluids->fluids[iFluid].get_q_rho_global();
		q_sx = gdata.fluids->fluids[iFluid].get_q_sx_global();
		q_sy = gdata.fluids->fluids[iFluid].get_q_sy_global();
		q_sz = gdata.fluids->fluids[iFluid].get_q_sz_global();
		if(!isothermal) {
			q_Eges = gdata.fluids->fluids[iFluid].get_q_Eges_global();
		}
#if (CRSWITCH_DUAL_ENERGY == CRONOS_ON)
		q_Add = gdata.fluids->fluids[iFluid].get_q_EAdd_global();
#endif

		for (int k = ibeg[2]; k <= iend[2]; k++){
//			cout << k << " " << fluidType << " " << CRONOS_MHD << endl;
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

						if(fluidType == CRONOS_MHD) {
							gdata.om[q_rho](i,j,k) = rhoLInit_MHD;
							gdata.om[q_sx](i,j,k) = VLXInit_MHD;
							gdata.om[q_sy](i,j,k) = VLYInit_MHD;
							gdata.om[q_sz](i,j,k) = VLZInit_MHD;
							// 					gdata.om[q_Bx](i,j,k) = BLXInit;
							// 					gdata.om[q_By](i,j,k) = BLYInit;
							// 					gdata.om[q_Bz](i,j,k) = BLZInit;
							gdata.om[q_Bx](i,j,k) = BLYInit*zz;
							gdata.om[q_By](i,j,k) = BLZInit*xx;
							gdata.om[q_Bz](i,j,k) = BLXInit*yy;

							if(!isothermal) {
								gdata.om[q_Eges](i,j,k) = PLInit_MHD;
							}
						} else {
							gdata.om[q_rho](i,j,k) = rhoLInit_HD;
							gdata.om[q_sx](i,j,k) = VLXInit_HD;
							gdata.om[q_sy](i,j,k) = VLYInit_HD;
							gdata.om[q_sz](i,j,k) = VLZInit_HD;

							if(!isothermal) {
								gdata.om[q_Eges](i,j,k) = PLInit_HD;
							}
						}


					} else {
						if(fluidType == CRONOS_MHD) {
							gdata.om[q_rho](i,j,k) = rhoRInit_MHD;
							gdata.om[q_sx](i,j,k) = VRXInit_MHD;
							gdata.om[q_sy](i,j,k) = VRYInit_MHD;
							gdata.om[q_sz](i,j,k) = VRZInit_MHD;
							// 					gdata.om[q_Bx](i,j,k) = BRXInit;
							// 					gdata.om[q_By](i,j,k) = BRYInit;
							// 					gdata.om[q_Bz](i,j,k) = BRZInit;
							gdata.om[q_Bx](i,j,k) = BRYInit*zz;
							gdata.om[q_By](i,j,k) = BRZInit*xx;
							gdata.om[q_Bz](i,j,k) = BRXInit*yy;

							if(!isothermal) {
								gdata.om[q_Eges](i,j,k) = PRInit_MHD;
							}
						} else {
							gdata.om[q_rho](i,j,k) = rhoRInit_HD;
							gdata.om[q_sx](i,j,k) = VRXInit_HD;
							gdata.om[q_sy](i,j,k) = VRYInit_HD;
							gdata.om[q_sz](i,j,k) = VRZInit_HD;

							if(!isothermal) {
								gdata.om[q_Eges](i,j,k) = PRInit_HD;
							}
						}
					}

					if(!isothermal) {
						if(gdata.om[q_Eges].getName() == "Temp") {
							gdata.om[q_Eges](i,j,k) /= gdata.om[q_rho](i,j,k);
						} else {
							gdata.om[q_Eges](i,j,k) *= 1./(gamma-1.);
						}
					}

#if (CRSWITCH_DUAL_ENERGY == CRONOS_ON)
					gdata.om[q_Eadd](i,j,k) = gdata.om[q_Eges](i,j,k)*pow(gdata.om[q_rho](i,j,k),1-gamma);
#endif

//					if(i==0 && j>390 && j<410 && k==0) {
//						cout << " Mag " << gdata.om[q_Bz](i,j,k) << " " << q_Bz << endl;
//					}

				}
			}
		}
	}
	gdata.om[q_Bx].rename("A_x");
	gdata.om[q_By].rename("A_y");
	gdata.om[q_Bz].rename("A_z");
}





#endif
