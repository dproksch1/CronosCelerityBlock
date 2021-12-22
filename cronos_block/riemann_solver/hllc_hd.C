#include "RiemannSolverHD.H"

//#define HLLCSOLVER_HYDRO_VEPS 1.e-120

#if (ENERGETICS != ISOTHERMAL)
HLLCSolver_Hydro::HLLCSolver_Hydro(const Data &gdata, int dir, int Fluid_Type) : RiemannSolverHD(gdata, dir, Fluid_Type) {
	veps = 1.e-120;

#if(FLUID_TYPE != CRONOS_MULTIFLUID)
	reset_Indices(gdata.fluid);
#else
	reset_Indices(gdata.fluids->fluids[0]);
#endif
	if(dir == 0) {
		this->qvPar = q_sx;
		this->qvP1  = q_sy;
		this->qvP2  = q_sz;
	} else if (dir == 1) {
		this->qvP2  = q_sx;
		this->qvPar = q_sy;
		this->qvP1  = q_sz;
	} else {
		this->qvP1  = q_sx;
		this->qvP2  = q_sy;
		this->qvPar = q_sz;
	}   
	gamma = value((char*)"Adiabatic_exponent");

}


void HLLCSolver_Hydro::get_NumFlux(const Data &gdata, const phys_fields_0D &pfM,
	const phys_fields_0D &pfP, num_fields_0D &f_num, int dir, int iFluid) const {
	//! Compute numerical flux from phys fluxes using hllc solver
	/*
	 * @PARAM pfM variables and fluxes on left-hand side of interface
	 * @PARAM pfP variables and fluxes on right-hand side of interface
	 */


	// In case of multifluid run, need to re-evaluate indices:
#if(FLUID_TYPE == CRONOS_MULTIFLUID)
	if(dir == 0) {
		this->qvPar = q_sx;
		this->qvP1  = q_sy;
		this->qvP2  = q_sz;
	} else if (dir == 1) {
		this->qvP2  = q_sx;
		this->qvPar = q_sy;
		this->qvP1  = q_sz;
	} else {
		this->qvP1  = q_sx;
		this->qvP2  = q_sy;
		this->qvPar = q_sz;
	}
#endif

	int n_omInt(pfM.get_num());

	/*-------------------------------------------------------
		  Compute fluxes for regions outside the Riemann fan:
		  -----------------------------------------------------*/

	if(f_num.v_ch_m <= 0.) { // Riemann fan going right

		for(int q=0; q<n_omInt; ++q) {
			f_num.flux_num[q] = f_num.flux_num[q] = pfM.flux_phys[q];
		}
		//#if(FLUID_TYPE==CRONOS_MULTIFLUID)
		//			fl.ptotal[iFluid](i) = pfR.ptotal(i-1);
		//#else
		f_num.ptotal_num = pfM.ptotal;
		//#endif

	} else if (f_num.v_ch_p <= 0.) { // Riemann fan going left

		for(int q=0; q<n_omInt; ++q) {
			f_num.flux_num[q] = pfP.flux_phys[q];
		}
		//#if(FLUID_TYPE==CRONOS_MULTIFLUID)
		//			fl.ptotal[iFluid](i) = pfL.ptotal(i);
		//#else
		f_num.ptotal_num = pfP.ptotal;
		//#endif

	} else {

		/*-------------------------------------------------------
			  The rest is taking place inside the Riemann fan
			  -------------------------------------------------------*/

		// If carbuncle test is used -> check whether cells i and i-1 are
		// flagged as lying within the vicinity of a strong shock
		if(gdata.use_carbuncleFlag) {
			// if both cells are flagged -> use hll for all variables instead
			if(pfM.carbuncle_flag*pfP.carbuncle_flag > 0) {
				for(int q=0; q<pfM.get_num(); ++q) {
					REAL fac = 1./(f_num.v_ch_p + f_num.v_ch_m+veps);

					f_num.flux_num[q] = (f_num.v_ch_m*pfP.flux_phys[q] +
							f_num.v_ch_p*pfM.flux_phys[q] -
							f_num.v_ch_m*f_num.v_ch_p*(pfP.uCon[q] - pfM.uCon[q]))*fac;
				}
				return;
			}
		}



		//REAL uConL[n_omInt], uConR[n_omInt];
		std::vector<REAL> uConL(n_omInt);
		std::vector<REAL> uConR(n_omInt);
		//REAL uPriL[n_omInt], uPriR[n_omInt];
		std::vector<REAL> uPriL(n_omInt);
		std::vector<REAL> uPriR(n_omInt);
		//REAL uConSL[n_omInt];
		std::vector<REAL> uConSL(n_omInt);
		//REAL uConSR[n_omInt];
		std::vector<REAL> uConSR(n_omInt);

		// Saving array values:

		for(int q=0; q<n_omInt; ++q) {
			/*
				In this context "L" and "M" mean (L)eft from the cell-face
				or on the (M)inus side
				The same holds for "R" and "P"
			 */
			uConR[q] = pfP.uCon[q];
			uConL[q] = pfM.uCon[q];
			uPriR[q] = pfP.uPri[q];
			uPriL[q] = pfM.uPri[q];
		}

#if (USE_COROTATION == CRONOS_ON)
		// Overwrite settings for co-rotation case:
		// uConR[qvPar] = uConL[q_rho]*uPriL[qvPar];
		// uConR[qvP1]  = uConL[q_rho]*uPriL[qvP1];
		// uConR[qvP2]  = uConL[q_rho]*uPriL[qvP2];
		uConR[qvPar] = uConR[q_rho]*uPriR[qvPar];
		uConR[qvP1]  = uConR[q_rho]*uPriR[qvP1];
		uConR[qvP2]  = uConR[q_rho]*uPriR[qvP2];
		uConL[qvPar] = uConL[q_rho]*uPriL[qvPar];
		uConL[qvP1]  = uConL[q_rho]*uPriL[qvP1];
		uConL[qvP2]  = uConL[q_rho]*uPriL[qvP2];
#endif


		REAL rhoL = uConL[q_rho];
		REAL rhoR = uConR[q_rho];

		REAL pThermL(pfM.ptherm);
		REAL pThermR(pfP.ptherm);

		REAL vParL = uPriL[qvPar];
		REAL vParR = uPriR[qvPar];

		REAL sParL = uConL[qvPar];
		REAL sParR = uConR[qvPar];

		// Characteristic velocities like in Toro book
		REAL vCharL = -f_num.v_ch_m;
		REAL vCharR =  f_num.v_ch_p;

		// Just the signal velocities relative to background flow
		REAL vSigL = vCharL - vParL;
		REAL vSigR = vCharR - vParR;

		// Velocity of contact discontinuity S*
		REAL idenom = 1./(rhoL*vSigL - rhoR*vSigR + veps);
		// Eq. (10.70)
		REAL vCharS = (pThermR - pThermL +
				sParL*vSigL - sParR*vSigR)*idenom;

		// Starred quantities according to Eq. (10.73)
		uConSL[q_rho] = uConL[q_rho]*vSigL/(vCharL - vCharS);
		uConSR[q_rho] = uConR[q_rho]*vSigR/(vCharR - vCharS);

		// if(i==54) {
		// 	cout << " uCon: " << uConSL[q_rho] << " ";
		// 	cout << uConSR[q_rho] << " " << vCharS << " " << vCharL << " " << vCharR << endl;
		// }

		// Velocities in middle region
		uConSL[qvPar] = uConSL[q_rho]*vCharS;
		// uConSL[qvP1 ] = uConSL[q_rho]*uPriSL[qvP1];
		// uConSL[qvP2 ] = uConSL[q_rho]*uPriSL[qvP2];
		uConSL[qvP1 ] = uConSL[q_rho]*uPriL[qvP1];
		uConSL[qvP2 ] = uConSL[q_rho]*uPriL[qvP2];

		uConSR[qvPar] = uConSR[q_rho]*vCharS;
		// uConSR[qvP1 ] = uConSR[q_rho]*uPriSR[qvP1];
		// uConSR[qvP2 ] = uConSR[q_rho]*uPriSR[qvP2];
		uConSR[qvP1 ] = uConSR[q_rho]*uPriR[qvP1];
		uConSR[qvP2 ] = uConSR[q_rho]*uPriR[qvP2];

		// Overall energy in middle region
		REAL egesL = uConL[q_Eges];
		uConSL[q_Eges] = uConSL[q_rho]*(egesL/rhoL + (vCharS - vParL)*
				(vCharS + pThermL/(rhoL*vSigL)));
		REAL egesR = uConR[q_Eges];
		uConSR[q_Eges] = uConSR[q_rho]*(egesR/rhoR + (vCharS - vParR)*
				(vCharS + pThermR/(rhoR*vSigR)));

#if(CRSWITCH_DUAL_ENERGY == CRONOS_ON)
	// Compute the entropy in the starred region for dual
		// energy description

		// Pressure in starred region (according to Eq. (10.36) in
		// Toro book) - should be identical for L and R values
		REAL pThermS = pThermL + rhoL*(vSigL - vParL)*(vCharS - vParL);

		// With this compute entropy:
		uConSL[q_Eadd] = pThermS/pow(uConSL[q_rho], gamma-1.);
		uConSR[q_Eadd] = pThermS/pow(uConSR[q_rho], gamma-1.);
#endif


		// Fluxes from Eq (10.71)
		if(vCharS >= 0.) {

			for(int q=0; q<n_omInt; ++q) {
				f_num.flux_num[q] = pfM.flux_phys[q] + vCharL*(uConSL[q] -
						uConL[q]);
			}

			// fl.ptotal(i) = ptotalL(i);
			f_num.ptotal_num = pfM.ptherm;

		} else {

			for(int q=0; q<n_omInt; ++q) {
				f_num.flux_num[q] = pfP.flux_phys[q] + vCharR*(uConSR[q] -
						uConR[q]);
			}
			// fl.ptotal(i) = ptotalR(i);
			f_num.ptotal_num = pfP.ptherm;

		}

	}
}

void get_NumFlux2(const Data &gdata, const phys_fields_0D &pfM,
	const phys_fields_0D &pfP, num_fields_0D &f_num, int dir, int iFluid) {
	//! Compute numerical flux from phys fluxes using hllc solver
	/*
	 * @PARAM pfM variables and fluxes on left-hand side of interface
	 * @PARAM pfP variables and fluxes on right-hand side of interface
	 */

	int qvPar, qvP1, qvP2;
	int q_sx = gdata.fluid.get_q_sx();
	int q_sy = gdata.fluid.get_q_sy();
	int q_sz = gdata.fluid.get_q_sz();

	// In case of multifluid run, need to re-evaluate indices:
//#if(FLUID_TYPE == CRONOS_MULTIFLUID)
	if(dir == 0) {
		qvPar = q_sx;
		qvP1  = q_sy;
		qvP2  = q_sz;
	} else if (dir == 1) {
		qvP2  = q_sx;
		qvPar = q_sy;
		qvP1  = q_sz;
	} else {
		qvP1  = q_sx;
		qvP2  = q_sy;
		qvPar = q_sz;
	}
//#endif

	int n_omInt(pfM.get_num());

	/*-------------------------------------------------------
		  Compute fluxes for regions outside the Riemann fan:
		  -----------------------------------------------------*/

	if(f_num.v_ch_m <= 0.) { // Riemann fan going right

		for(int q=0; q<n_omInt; ++q) {
			f_num.flux_num[q] = f_num.flux_num[q] = pfM.flux_phys[q];
		}
		//#if(FLUID_TYPE==CRONOS_MULTIFLUID)
		//			fl.ptotal[iFluid](i) = pfR.ptotal(i-1);
		//#else
		f_num.ptotal_num = pfM.ptotal;
		//#endif

	} else if (f_num.v_ch_p <= 0.) { // Riemann fan going left

		for(int q=0; q<n_omInt; ++q) {
			f_num.flux_num[q] = pfP.flux_phys[q];
		}
		//#if(FLUID_TYPE==CRONOS_MULTIFLUID)
		//			fl.ptotal[iFluid](i) = pfL.ptotal(i);
		//#else
		f_num.ptotal_num = pfP.ptotal;
		//#endif

	} else {

		/*-------------------------------------------------------
			  The rest is taking place inside the Riemann fan
			  -------------------------------------------------------*/

		// If carbuncle test is used -> check whether cells i and i-1 are
		// flagged as lying within the vicinity of a strong shock
		if(gdata.use_carbuncleFlag) {
			// if both cells are flagged -> use hll for all variables instead
			if(pfM.carbuncle_flag*pfP.carbuncle_flag > 0) {
				for(int q=0; q<pfM.get_num(); ++q) {
					REAL fac = 1./(f_num.v_ch_p + f_num.v_ch_m+ HLLCSOLVER_HYDRO_VEPS);

					f_num.flux_num[q] = (f_num.v_ch_m*pfP.flux_phys[q] +
							f_num.v_ch_p*pfM.flux_phys[q] -
							f_num.v_ch_m*f_num.v_ch_p*(pfP.uCon[q] - pfM.uCon[q]))*fac;
				}
				return;
			}
		}



		//REAL uConL[n_omInt], uConR[n_omInt];
		std::vector<REAL> uConL(n_omInt);
		std::vector<REAL> uConR(n_omInt);
		//REAL uPriL[n_omInt], uPriR[n_omInt];
		std::vector<REAL> uPriL(n_omInt);
		std::vector<REAL> uPriR(n_omInt);
		//REAL uConSL[n_omInt];
		std::vector<REAL> uConSL(n_omInt);
		//REAL uConSR[n_omInt];
		std::vector<REAL> uConSR(n_omInt);

		// Saving array values:

		for(int q=0; q<n_omInt; ++q) {
			/*
				In this context "L" and "M" mean (L)eft from the cell-face
				or on the (M)inus side
				The same holds for "R" and "P"
			 */
			uConR[q] = pfP.uCon[q];
			uConL[q] = pfM.uCon[q];
			uPriR[q] = pfP.uPri[q];
			uPriL[q] = pfM.uPri[q];
		}

		int q_rho = gdata.fluid.get_q_rho();

#if (USE_COROTATION == CRONOS_ON)
		// Overwrite settings for co-rotation case:
		// uConR[qvPar] = uConL[q_rho]*uPriL[qvPar];
		// uConR[qvP1]  = uConL[q_rho]*uPriL[qvP1];
		// uConR[qvP2]  = uConL[q_rho]*uPriL[qvP2];
		uConR[qvPar] = uConR[q_rho]*uPriR[qvPar];
		uConR[qvP1]  = uConR[q_rho]*uPriR[qvP1];
		uConR[qvP2]  = uConR[q_rho]*uPriR[qvP2];
		uConL[qvPar] = uConL[q_rho]*uPriL[qvPar];
		uConL[qvP1]  = uConL[q_rho]*uPriL[qvP1];
		uConL[qvP2]  = uConL[q_rho]*uPriL[qvP2];
#endif


		REAL rhoL = uConL[q_rho];
		REAL rhoR = uConR[q_rho];

		REAL pThermL(pfM.ptherm);
		REAL pThermR(pfP.ptherm);

		REAL vParL = uPriL[qvPar];
		REAL vParR = uPriR[qvPar];

		REAL sParL = uConL[qvPar];
		REAL sParR = uConR[qvPar];

		// Characteristic velocities like in Toro book
		REAL vCharL = -f_num.v_ch_m;
		REAL vCharR =  f_num.v_ch_p;

		// Just the signal velocities relative to background flow
		REAL vSigL = vCharL - vParL;
		REAL vSigR = vCharR - vParR;

		// Velocity of contact discontinuity S*
		REAL idenom = 1./(rhoL*vSigL - rhoR*vSigR + HLLCSOLVER_HYDRO_VEPS);
		// Eq. (10.70)
		REAL vCharS = (pThermR - pThermL +
				sParL*vSigL - sParR*vSigR)*idenom;

		// Starred quantities according to Eq. (10.73)
		uConSL[q_rho] = uConL[q_rho]*vSigL/(vCharL - vCharS);
		uConSR[q_rho] = uConR[q_rho]*vSigR/(vCharR - vCharS);

		// if(i==54) {
		// 	cout << " uCon: " << uConSL[q_rho] << " ";
		// 	cout << uConSR[q_rho] << " " << vCharS << " " << vCharL << " " << vCharR << endl;
		// }

		// Velocities in middle region
		uConSL[qvPar] = uConSL[q_rho]*vCharS;
		// uConSL[qvP1 ] = uConSL[q_rho]*uPriSL[qvP1];
		// uConSL[qvP2 ] = uConSL[q_rho]*uPriSL[qvP2];
		uConSL[qvP1 ] = uConSL[q_rho]*uPriL[qvP1];
		uConSL[qvP2 ] = uConSL[q_rho]*uPriL[qvP2];

		uConSR[qvPar] = uConSR[q_rho]*vCharS;
		// uConSR[qvP1 ] = uConSR[q_rho]*uPriSR[qvP1];
		// uConSR[qvP2 ] = uConSR[q_rho]*uPriSR[qvP2];
		uConSR[qvP1 ] = uConSR[q_rho]*uPriR[qvP1];
		uConSR[qvP2 ] = uConSR[q_rho]*uPriR[qvP2];

		// Overall energy in middle region
		int q_Eges = gdata.fluid.get_q_Eges();
		REAL egesL = uConL[q_Eges];
		uConSL[q_Eges] = uConSL[q_rho]*(egesL/rhoL + (vCharS - vParL)*
				(vCharS + pThermL/(rhoL*vSigL)));
		REAL egesR = uConR[q_Eges];
		uConSR[q_Eges] = uConSR[q_rho]*(egesR/rhoR + (vCharS - vParR)*
				(vCharS + pThermR/(rhoR*vSigR)));

#if(CRSWITCH_DUAL_ENERGY == CRONOS_ON)
	// Compute the entropy in the starred region for dual
		// energy description

		// Pressure in starred region (according to Eq. (10.36) in
		// Toro book) - should be identical for L and R values
		REAL pThermS = pThermL + rhoL*(vSigL - vParL)*(vCharS - vParL);

		// With this compute entropy:
		REAL gamma = value((char*)"Adiabatic_exponent");
		uConSL[q_Eadd] = pThermS/pow(uConSL[q_rho], gamma-1.);
		uConSR[q_Eadd] = pThermS/pow(uConSR[q_rho], gamma-1.);
#endif


		// Fluxes from Eq (10.71)
		if(vCharS >= 0.) {

			for(int q=0; q<n_omInt; ++q) {
				f_num.flux_num[q] = pfM.flux_phys[q] + vCharL*(uConSL[q] -
						uConL[q]);
			}

			// fl.ptotal(i) = ptotalL(i);
			f_num.ptotal_num = pfM.ptherm;

		} else {

			for(int q=0; q<n_omInt; ++q) {
				f_num.flux_num[q] = pfP.flux_phys[q] + vCharR*(uConSR[q] -
						uConR[q]);
			}
			// fl.ptotal(i) = ptotalR(i);
			f_num.ptotal_num = pfP.ptherm;

		}

	}
}

#endif
