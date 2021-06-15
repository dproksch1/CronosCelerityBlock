#include "RiemannSolverMHD.H"
#include <iomanip>
#include <stdlib.h>

int sign(REAL number) {
	return ((number > 0) - (number < 0));
}

//HLLDSolver::HLLDSolver(const Data &gdata, int dir, int Fluid_Type) : RiemannSolverMHD(gdata, dir, Fluid_Type) {
//
//#if(FLUID_TYPE != CRONOS_MULTIFLUID)
//	reset_Indices(gdata.fluid);
//#else
//	reset_Indices(gdata.fluids->fluids[0]);
//#endif
////	this->q_Eges = q_Eges;
//	this->qBMin = q_Bx;
//	this->qBMax = q_Bz;
//	this->veps = 1.e-120;
//	this->loceps = 1.e-4;
//	gamma = value((char*)"Adiabatic_exponent");
//
//	if(dir == 0) {
//		this->qvPar = q_sx;
//		this->qvP1  = q_sy;
//		this->qvP2  = q_sz;
//		this->qBPar = q_Bx;
//		this->qBP1  = q_By;
//		this->qBP2  = q_Bz;
//	} else if (dir == 1) {
//		this->qvP2  = q_sx;
//		this->qvPar = q_sy;
//		this->qvP1  = q_sz;
//		this->qBP2  = q_Bx;
//		this->qBPar = q_By;
//		this->qBP1  = q_Bz;
//	} else {
//		this->qvP1  = q_sx;
//		this->qvP2  = q_sy;
//		this->qvPar = q_sz;
//		this->qBP1  = q_Bx;
//		this->qBP2  = q_By;
//		this->qBPar = q_Bz;
//	}		
//
//}
//
//void HLLDSolver::get_NumFlux(Queue /*queue*/, Data &gdata,
//                             phys_fields_1D &pfL,
//                             phys_fields_1D &pfR,
//                             fields_1D &fl, const int &dir, int iFluid) const {
//
////	if(dir == 0) {
////		this->qvPar = q_sx;
////		this->qvP1  = q_sy;
////		this->qvP2  = q_sz;
////		this->qBPar = q_Bx;
////		this->qBP1  = q_By;
////		this->qBP2  = q_Bz;
////	} else if (dir == 1) {
////		this->qvP2  = q_sx;
////		this->qvPar = q_sy;
////		this->qvP1  = q_sz;
////		this->qBP2  = q_Bx;
////		this->qBPar = q_By;
////		this->qBP1  = q_Bz;
////	} else {
////		this->qvP1  = q_sx;
////		this->qvP2  = q_sy;
////		this->qvPar = q_sz;
////		this->qBP1  = q_Bx;
////		this->qBP2  = q_By;
////		this->qBPar = q_Bz;
////	}
//
//	// In case of multifluid run, need to re-evaluate indices:
//#if(FLUID_TYPE == CRONOS_MULTIFLUID)
//	if(dir == 0) {
//		this->qvPar = q_sx;
//		this->qvP1  = q_sy;
//		this->qvP2  = q_sz;
//		this->qBPar = q_Bx;
//		this->qBP1  = q_By;
//		this->qBP2  = q_Bz;
//	} else if (dir == 1) {
//		this->qvP2  = q_sx;
//		this->qvPar = q_sy;
//		this->qvP1  = q_sz;
//		this->qBP2  = q_Bx;
//		this->qBPar = q_By;
//		this->qBP1  = q_Bz;
//	} else {
//		this->qvP1  = q_sx;
//		this->qvP2  = q_sy;
//		this->qvPar = q_sz;
//		this->qBP1  = q_Bx;
//		this->qBP2  = q_By;
//		this->qBPar = q_Bz;
//	}
//#endif
//	int n_omInt(fl.get_num());
//
//	if(ENERGETICS == FULL) {
//
//		double eps = 1.e-8;
//
//		for (int i = -1; i <= gdata.mx[dir]+1; ++i){
//
//
//			//REAL uConL[n_omInt], uConR[n_omInt];
//			std::vector<REAL> uConL(n_omInt);
//			std::vector<REAL> uConR(n_omInt);
//			//REAL uPriL[n_omInt], uPriR[n_omInt];
//			std::vector<REAL> uPriL(n_omInt);
//			std::vector<REAL> uPriR(n_omInt);
//			NumArray<REAL> uConSL(n_omInt), uPriSL(n_omInt);
//			NumArray<REAL> uConSR(n_omInt), uPriSR(n_omInt);
////			REAL uConSL[n_omInt], uPriSL[n_omInt];
////			REAL uConSR[n_omInt], uPriSR[n_omInt];
//
//			// Saving array values:
//			for(int q=0; q<n_omInt; ++q) {
//				/*
//				   Be aware of the two differnt meanings of "L"
//				   For hlld "L" means "on the Left hand side of the cell interface"
//				   For the reconstruction it means on the left side of the cell
//				*/
//				uConR[q] = pfL.uConORIG[q](i  );
//				uConL[q] = pfR.uConORIG[q](i-1);
//				uPriR[q] = pfL.uPriORIG[q](i  );
//				uPriL[q] = pfR.uPriORIG[q](i-1);
//			}
//
//			REAL rhoL = uConL[q_rho];
//			REAL rhoR = uConR[q_rho];
//
//			REAL vParL = uPriL[qvPar];
//			REAL vParR = uPriR[qvPar];
//
//			REAL sParL = uConL[qvPar];
//			REAL sParR = uConR[qvPar];
//
//			REAL BPar = uConL[qBPar];
//
//			REAL sqrBL = sqr(uConL[qBPar]) + sqr(uConL[qBP1 ]) + sqr(uConL[qBP2 ]);
//			REAL sqrBR = sqr(uConR[qBPar]) + sqr(uConR[qBP1 ]) + sqr(uConR[qBP2 ]);
//
//			REAL presL = pfR.pthermORIG(i-1);
//			REAL presR = pfL.pthermORIG(i  );
//
//			ptotalL(i) = pfR.pthermORIG(i-1) + 0.5*sqrBL;
//			ptotalR(i) = pfL.pthermORIG(i  ) + 0.5*sqrBR;
//
//
//			//if((verbosity==52 && i==50)) {
//			//	cout << endl << " ptotal left ";
//			//	cout << ptotalL(i) << " ";
//			//	cout << sqrBL << " ";
//			//	cout << presL << " ";
//			//	cout << uPriL[q_Eges] << " ";
//			//	cout << endl;
//			//}
//			//if((verbosity==62 && i==50)) {
//			//	cout << endl << " ptotal righ  ";
//			//	cout << ptotalR(i) << " ";
//			//	cout << sqrBR << " ";
//			//	cout << presR << " ";
//			//	cout << uPriR[q_Eges] << " ";
//			//	cout << endl;
//			//}
//
//
//
//			/*-------------------------------------------------------
//				Compute maximum speeds - see Eq. (67) in Miyoshi & Kusano (2005)
//			 -----------------------------------------------------*/
//
//			double sqrBPL = sqrBL + gamma*presL;
//			double sqrBPR = sqrBR + gamma*presR;
//
//			// fast mode speed estimates
//			double cfL = sqrt( (sqrBPL + sqrt( sqr(sqrBPL) - 4.*gamma*presL*sqr(uConL[qBPar])) )/(2.*rhoL) );
//			double cfR = sqrt( (sqrBPR + sqrt( sqr(sqrBPR) - 4.*gamma*presR*sqr(uConR[qBPar])) )/(2.*rhoR) );
//
//			REAL vCharL = std::min(vParL, vParR) - std::max(cfL, cfR);
//			REAL vCharR = std::max(vParL, vParR) + std::max(cfL, cfR);
//
//			//if((verbosity==34 && i==66)) {
//			//	cout << " vChars kaputt ";
//			//	cout << cfL << " " << cfR << " " << vParL << " ";
//			//	cout << sqr(BPar) << " "<< sqrBL << " ";
//			//	cout << gamma << " " << presL << " " << gamma*presL << " ";
//			//	cout << endl;
//			//}
//
//
//			/*-------------------------------------------------------
//			  Compute fluxes for regions outside the Riemann fan:
//			  -----------------------------------------------------*/
//		
//			if(fl.v_ch_mORIG(i) <= 0.) { // Riemann fan going right
//
//				for(int q=0; q<n_omInt; ++q) {
//					fl.fluxORIG[q](i) = pfR.fluxORIG[q](i-1);
//				}
//				fl.ptotalORIG(i) = pfR.ptotalORIG(i-1);
//				//if(verbosity > 0 && false) {
//				//	cout << " cell " << i << " in region l " << endl;
//				//}
////				if(gdata.time>0.0999 && gdata.time<0.1001 && dir==0) {
////					cout << i << " " << gdata.getCen_x(i) << " left " << endl;
////				}
//
//			} else if (fl.v_ch_pORIG(i) <= 0.) { // Riemann fan going left
//
//				for(int q=0; q<n_omInt; ++q) {
//					fl.fluxORIG[q](i) = pfL.fluxORIG[q](i);
//				}
//				fl.ptotalORIG(i) = pfL.ptotalORIG(i);
//				//if(verbosity > 0 && false) {
//				//	cout << " cell " << i << " in region r " << endl;
//				//}
////				if(gdata.time>0.0999 && gdata.time<0.1001 && dir==0) {
////					cout << i << " " << gdata.getCen_x(i) << " right " << endl;
////				}
//
//			} else {
//
//				/*-------------------------------------------------------
//				  The rest is taking place inside the Riemann fan
//				  -------------------------------------------------------*/
//
//				// If carbuncle test is used -> check whether cells i and i-1 are
//				// flagged as lying within the vicinity of a strong shock
//				if(gdata.use_carbuncleFlag) {
//					// if both cells are flagged -> use hll for all variables instead
//					if(fl.carbuncle_flag(i)*fl.carbuncle_flag(i-1) > 0) {
//						for(int q=0; q<fl.get_num(); ++q) {
//							REAL fac = 1./(fl.v_ch_pORIG(i)+fl.v_ch_mORIG(i)+veps);
//
//							fl.fluxORIG[q](i) = (fl.v_ch_mORIG(i)*pfL.fluxORIG[q](i) +
//									fl.v_ch_pORIG(i)*pfR.fluxORIG[q](i-1) -
//									fl.v_ch_mORIG(i)*fl.v_ch_pORIG(i)*(pfL.uConORIG[q](i  )-
//											pfR.uConORIG[q](i-1)))*fac;
//						}
//						continue;
//					}
//				}
//
//			
////				REAL uConL[n_omInt], uConR[n_omInt];
////				REAL uPriL[n_omInt], uPriR[n_omInt];
////				REAL uConSL[n_omInt], uPriSL[n_omInt];
////				REAL uConSR[n_omInt], uPriSR[n_omInt];
////
////				// Saving array values:
////
////				for(int q=0; q<n_omInt; ++q) {
////					/*
////					   Be aware of the two differnt meanings of "L"
////					   For hlld "L" means "on the Left hand side of the cell interface"
////					   For the reconstruction it means on the left side of the cell
////					*/
////					uConR[q] = pfL.uCon[q](i  );
////					uConL[q] = pfR.uCon[q](i-1);
////					uPriR[q] = pfL.uPri[q](i  );
////					uPriL[q] = pfR.uPri[q](i-1);
////				}
//      
//#if (USE_COROTATION == CRONOS_ON)
//				// Overwrite settings for co-rotation case:
//				// uConR[qvPar] = uConL[q_rho]*uPriL[qvPar];
//				// uConR[qvP1]  = uConL[q_rho]*uPriL[qvP1];
//				// uConR[qvP2]  = uConL[q_rho]*uPriL[qvP2];
//				uConR[qvPar] = uConR[q_rho]*uPriR[qvPar];
//				uConR[qvP1]  = uConR[q_rho]*uPriR[qvP1];
//				uConR[qvP2]  = uConR[q_rho]*uPriR[qvP2];
//				uConL[qvPar] = uConL[q_rho]*uPriL[qvPar];
//				uConL[qvP1]  = uConL[q_rho]*uPriL[qvP1];
//				uConL[qvP2]  = uConL[q_rho]*uPriL[qvP2];
//#endif
//
////				REAL rhoL = uConL[q_rho];
////				REAL rhoR = uConR[q_rho];
////
////				REAL vParL = uPriL[qvPar];
////				REAL vParR = uPriR[qvPar];
////
////				REAL sParL = uConL[qvPar];
////				REAL sParR = uConR[qvPar];
////
////				REAL BPar = uConL[qBPar];
////
////				REAL pBL = sqr(uConL[qBPar]) + sqr(uConL[qBP1 ]) + sqr(uConL[qBP2 ]);
////				REAL pBR = sqr(uConR[qBPar]) + sqr(uConR[qBP1 ]) + sqr(uConR[qBP2 ]);
//
//				ptotalL(i) = pfR.pthermORIG(i-1) + 0.5*sqrBL;
//				ptotalR(i) = pfL.pthermORIG(i  ) + 0.5*sqrBR;
//
////				REAL vCharL = -fl.v_ch_m(i);
////				REAL vCharR =  fl.v_ch_p(i);
//
//				// Need to compute the
//
//
//				REAL vSigL = vCharL - vParL;
//				REAL vSigR = vCharR - vParR;
//
//
//				REAL idenom = 1./(vSigR*rhoR - vSigL*rhoL + veps);
//
//
//				// Eq. (38)
//				// S_M according to eq 38 in hlld paper
//				REAL vCharM = (vSigR*sParR - vSigL*sParL -
//				               ptotalR(i) + ptotalL(i))*idenom;
//
//				//if((verbosity==42 && i==65)) {
//				//	cout << " vCharM ";
//				//	cout << vSigR << " " << vSigL << " ";
//				//	cout << sParR << " " << sParL << " ";
//				//	cout << ptotalR(i) << " " << ptotalL(i) << " ";
//				//	cout << idenom << " " << vCharM << " ";
//				//	cout << vSigR*sParR - vSigL*sParL -
//				//               ptotalR(i) + ptotalL(i) << " ";
//				//	cout << endl;
//				//}
//				//if((verbosity==24 && i==gdata.mx[1]-65+1)) {
//				//	cout << " vCharM ";
//				//	cout << vSigL << " " << vSigR << " ";
//				//	cout << sParL << " " << sParR << " ";
//				//	cout << ptotalL(i) << " " << ptotalR(i) << " ";
//				//	cout << idenom << " " << vCharM << " ";
//				//	cout << vSigR*sParR - vSigL*sParL -
//				//               ptotalR(i) + ptotalL(i) << " ";
//				//	cout << endl;
//				//}
//
//				//if((verbosity==52 && i==50)) {
//				//	cout << " vCharM left ";
//				//	cout << vSigR << " " << vSigL << " ";
//				//	cout << sParR << " " << sParL << " ";
//				//	cout << ptotalR(i) << " " << ptotalL(i) << " ";
//				//	cout << idenom << " " << vCharM << " ";
//				//	cout << vSigR*sParR - vSigL*sParL -
//				//               ptotalR(i) + ptotalL(i) << " ";
//				//	cout << endl;
//				//}
//				//if((verbosity==62 && i==50)) {
//				//	cout << " vCharM righ ";
//				//	cout << vSigL << " " << vSigR << " ";
//				//	cout << sParL << " " << sParR << " ";
//				//	cout << ptotalL(i) << " " << ptotalR(i) << " ";
//				//	cout << idenom << " " << vCharM << " ";
//				//	cout << vSigR*sParR - vSigL*sParL -
//				//               ptotalR(i) + ptotalL(i) << " ";
//				//	cout << endl;
//				//}
//
//
//
//				/*----------------------------------------------------------
//				  Start computation of plasma state in Regions L* and R*
//				  --------------------------------------------------------*/
//
//				// Total pressure in Regions L* and R* - see Eq. (41)
//				REAL ptotalS = (vSigR*rhoR*ptotalL(i) - vSigL*rhoL*ptotalR(i) +
//				                rhoL*rhoR*vSigL*vSigR*(vParR - vParL))*idenom;
//
//
//				// Mass denisty in Regions L* and R* - Eq. (43)
//				uConSL[q_rho] = uConL[q_rho]*vSigL/(vCharL - vCharM);
//				uConSR[q_rho] = uConR[q_rho]*vSigR/(vCharR - vCharM);
//
//
//
//				//if((verbosity==52 && i==50)) {
//				//	cout << " rho star left ";
//				//	cout << uConSR[q_rho] << " ";
//				//	cout << uConR[q_rho] << " ";
//				//	cout << vSigR << " ";
//				//	cout << vCharM << " ";
//				//	cout << endl;
//				//}
//				//if((verbosity==62 && i==50)) {
//				//	cout << " rho star right ";
//				//	cout << uConSL[q_rho] << " ";
//				//	cout << uConL[q_rho] << " ";
//				//	cout << vSigL << " ";
//				//	cout << vCharM << " ";
//				//	cout << endl;
//				//}
//
//
//				// Needed lateron
//				REAL sqrtRhoSL = sqrt(uConSL[q_rho]);
//				REAL sqrtRhoSR = sqrt(uConSR[q_rho]);
//
//				// Characteristic velocities of rim of regions L* and R* - see Eq. (51)
//				REAL vCharSL = vCharM - std::abs(BPar)/sqrtRhoSL;
//				REAL vCharSR = vCharM + std::abs(BPar)/sqrtRhoSR;
//
//
//				//if(verbosity==52 && i==50) {
//				//	cout << " speeds bot " << i << " ";
//				//	cout << vCharL << " " << vCharSL << " " << vCharM << " " << vCharSR << " " << vCharR << " " << std::abs(BPar) << endl;
//				//}
//				//if(verbosity==22 && i==50) {
//				//	cout << " speeds " << i << " ";
//				//	cout << vCharL << " " << vCharSL << " " << vCharM << " " << vCharSR << " " << vCharR << " " << std::abs(BPar) << endl;
//				//}
//				//if(verbosity==62 && i==50) {
//				//	cout << " speeds top " << i << " ";
//				//	cout << vCharR << " " << vCharSR << " " << vCharM << " " << vCharSL << " " << vCharL << " " << std::abs(BPar) << endl;
//				//}
//
//
////				if(verbosity>0) {
//				//if((verbosity==42 && i==65+1)) {
//				//	cout << " speeds " << i << " ";
//				//	cout << vCharL << " " << vCharSL << " " << vCharM << " " << vCharSR << " " << vCharR << " " << std::abs(BPar) << endl;
//				//}
//				//if((verbosity==24 && i==gdata.mx[1]-65)) {
//				//	cout << " speeds " << i << " ";
//				//	cout << -vCharR << " " << -vCharSR << " " << -vCharM << " " << -vCharSL << " " << -vCharL << " " << std::abs(BPar) << endl;
//				//}
////				if((verbosity==34 && i==51)) {
////					cout << " speeds kaputt " << i << " ";
////					cout << -vCharR << " " << -vCharSR << " " << -vCharM << " " << -vCharSL << " " << -vCharL << " " << std::abs(BPar) << endl;
////					cout << "               " << sqrtRhoSL << " ";
////					cout << uConSL[q_rho] << " " << uConL[q_rho] << " ";
////					cout << uConSR[q_rho] << " " << uConR[q_rho] << " ";
////					cout << endl;
////				}
//
////				if((verbosity==34 && i==51)) {
////					cout << " Left and right stuff ";
////					cout << uConL[q_rho] << " " << uConR[q_rho] << " ";
////					cout << uConL[q_sx] << " " << uConR[q_sx] << " ";
////					cout << uConL[q_sy] << " " << uConR[q_sy] << " ";
////					cout << uConL[q_sz] << " " << uConR[q_sz] << " ";
////					cout << endl;
////					cout << "                     ";
////					cout << uConL[q_Bx] << " " << uConR[q_Bx] << " ";
////					cout << uConL[q_By] << " " << uConR[q_By] << " ";
////					cout << uConL[q_Bz] << " " << uConR[q_Bz] << " ";
////					cout << uConL[q_Eges] << " " << uConR[q_Eges] << " ";
////					cout << ptotalL[i] << "  "<< ptotalR[i] << " ";
////					cout << endl;
////				}
//
//				/*
//
//				  According to Mignone a destinction has to be made in the
//				  following cases:
//
//				  Here we just compute intermediate values for the *-values of
//				  the individual components of the magnetic induction
//
//				*/
//
//				bool use_hllc = false;
//				if((vCharSL - vCharL) <  loceps*(vCharM - vCharL)) {
//					use_hllc = true;
//				}
//
//				if((vCharSR - vCharR) > -loceps*(vCharR - vCharM)) {
//					use_hllc = true;
//				}
//				
////				if(verbosity>0 && i==27) {
////					cout << " use " << use_hllc << " " << vCharL << " " << vCharM << " " << vCharR << endl;
////				}
//
//				use_hllc = false;
//
//				if(use_hllc) {
//
//					idenom = 1./(vCharR - vCharL);
//
//					for(int q=qBMin; q<=qBMax; ++q) {
//
//						REAL value = (vCharR*uConR[q] - vCharL*uConL[q] +
//						              pfR.fluxORIG[q](i-1) - pfL.fluxORIG[q](i  ))*idenom;
//						uConSL[q] = value;
//						uConSR[q] = value;
//
//					}
//
//					// By this a computation of Regions ** is avoided
//					vCharSL = vCharSR = vCharM;
//
//				} else {
//
//					// Using hlld
//
//					REAL facL = (rhoL*sqr(vSigL) - sqr(BPar)) / 
//						(rhoL*vSigL*(vCharL - vCharM) - sqr(BPar));
//
//
//					// Compute components of mag induction 
//					// see also Eqs. (45) and (47)
//					uConSL[qBPar] = BPar;
//
//					if(std::abs(vCharM) < 1.e-120 && std::abs(vCharL - vSigL) < 1.e-120) {
//						uConSL[qBP1 ] = uConL[qBP1 ];
//						uConSL[qBP2 ] = uConL[qBP2 ];
//					} else {
//						uConSL[qBP1 ] = uConL[qBP1 ]*facL;
//						uConSL[qBP2 ] = uConL[qBP2 ]*facL;
//					}
//
//
//					REAL facR = (rhoR*sqr(vSigR) - sqr(BPar)) / 
//						(rhoR*vSigR*(vCharR - vCharM) - sqr(BPar));
//
//					uConSR[qBPar] = BPar;
//
//					if(std::abs(vCharM) < 1.e-120 && std::abs(vCharR - vSigR) < 1.e-120) {
//						uConSR[qBP1 ] = uConR[qBP1 ];
//						uConSR[qBP2 ] = uConR[qBP2 ];
//					} else {
//						uConSR[qBP1 ] = uConR[qBP1 ]*facR;
//						uConSR[qBP2 ] = uConR[qBP2 ]*facR;
//					}
//					
//					//if((verbosity==34 && i==51 )) {
//					//	cout << " hlld ";
//					//	cout << uConSR[qBP1] << " ";
//					//	cout << uConSL[qBP1] << " ";
//					//	cout << facL << " ";
//					//	cout << endl;
//					//}
//
//				}
//
//				if(verbosity>0 && i==27) {
//					cout << " Comps mag " << uConSR[qBP1 ] << " " << uConSR[qBP2 ] << " ";
//					cout << endl;
//				}
//
//				// Compute perp components of velocity in regions L* and R*
//				// (this is simplified by using the above) see Eqs. (44) &
//				// (46)
//
//				REAL facL = BPar/(rhoL*vSigL);
//				REAL facR = BPar/(rhoR*vSigR);
//
//				uPriSL[qvP1] = uPriL[qvP1] - (uConSL[qBP1] - uConL[qBP1])*facL;
//				uPriSR[qvP1] = uPriR[qvP1] - (uConSR[qBP1] - uConR[qBP1])*facR;
//
//
//				uPriSL[qvP2] = uPriL[qvP2] - (uConSL[qBP2] - uConL[qBP2])*facL;
//				uPriSR[qvP2] = uPriR[qvP2] - (uConSR[qBP2] - uConR[qBP2])*facR;
//
//
////				if((verbosity==34 && i==51)) {
////					cout << " innermost ";
////					cout << uPriSL[qvP1] << " ";
////					cout << uPriSR[qvP1] << " ";
////					cout << uConSR[qBP1] << " ";
////					cout << uConSL[qBP1] << " ";
////					cout << endl;
////				}
//
//				// DEBUG
//#if 1
//				// check for degenerate case (nominator of Eq (44) etc is zero)
//
//				if(std::abs(rhoL*vSigL*(vCharL - vCharM) - sqr(BPar)) < eps*ptotalS) { // degenerate case
//					uPriSL[qvP1] = uPriL[qvP1];
//					uPriSL[qvP2] = uPriL[qvP2];
//
//					uConSL[qBP1] = uConL[qBP1];
//					uConSL[qBP2] = uConL[qBP2];
//				} else {
//					// Direct implementation as given in paper
//					// start with perpendicular velocity
//					facL = BPar*(vCharM - vParL) /
//							(rhoL*vSigL*(vCharL - vCharM) - sqr(BPar));
//					uPriSL[qvP1] = uPriL[qvP1] - (uConL[qBP1])*facL;
//					uPriSL[qvP2] = uPriL[qvP2] - (uConL[qBP2])*facL;
////					uPriSL[qvP1] = uPriL[qvP1] - (uConSL[qBP1])*facL;
////					uPriSL[qvP2] = uPriL[qvP2] - (uConSL[qBP2])*facL;
//
//					// now, magnetic field:
//					facL = (rhoL*sqr(vSigL) - sqr(BPar)) /
//							(rhoL*vSigL*(vCharL - vCharM) - sqr(BPar));
//					uConSL[qBP1 ] = uConL[qBP1 ]*facL;
//					uConSL[qBP2 ] = uConL[qBP2 ]*facL;
//				}
//
//
//				if(std::abs(rhoR*vSigR*(vCharR - vCharM) - sqr(BPar)) < eps*ptotalS) { // degenerate case
//					uPriSR[qvP1] = uPriR[qvP1];
//					uPriSR[qvP2] = uPriR[qvP2];
//
//					uConSR[qBP1] = uConR[qBP1];
//					uConSR[qBP2] = uConR[qBP2];
//				} else {
//					facR = BPar*(vCharM - vParR) /
//							(rhoR*vSigR*(vCharR - vCharM) - sqr(BPar));
//					uPriSR[qvP1] = uPriR[qvP1] - (uConR[qBP1])*facR;
//					uPriSR[qvP2] = uPriR[qvP2] - (uConR[qBP2])*facR;
////					uPriSR[qvP1] = uPriR[qvP1] - (uConSR[qBP1])*facR;
////					uPriSR[qvP2] = uPriR[qvP2] - (uConSR[qBP2])*facR;
//
//					//if(verbosity==22 && (i>34 && i<38)) {
//					//	cout << " Rechte Velo " << i << " ";
//					//	cout << facR << " ";
//					//	cout << (rhoR*vSigR*(vCharR - vCharM) - sqr(BPar)) << " ";
//					//	cout << eps*ptotalS << " ";
//					//	cout << uPriSR[qvP2] << " ";
//					//	cout << endl;
//
//					//}
//
//					// now magnetic field
//					facR = (rhoR*sqr(vSigR) - sqr(BPar)) /
//							(rhoR*vSigR*(vCharR - vCharM) - sqr(BPar));
//					uConSR[qBP1 ] = uConR[qBP1 ]*facR;
//					uConSR[qBP2 ] = uConR[qBP2 ]*facR;
//
//
//				}
//#endif
//
//				// Get corresponding values for the momentum:
//				uConSL[qvPar] = uConSL[q_rho]*vCharM;
//				uConSL[qvP1 ] = uConSL[q_rho]*uPriSL[qvP1];
//				uConSL[qvP2 ] = uConSL[q_rho]*uPriSL[qvP2];
//
//				uConSR[qvPar] = uConSR[q_rho]*vCharM;
//				uConSR[qvP1 ] = uConSR[q_rho]*uPriSR[qvP1];
//				uConSR[qvP2 ] = uConSR[q_rho]*uPriSR[qvP2];
//				
//
//				// Compute Energy in stared state (see Eq. (48))
//				REAL egesL = uConL[q_Eges];
//				REAL vScalBL  = (uPriL[qvPar]*BPar +
//				                 uPriL[qvP1 ]*uConL[qBP1 ] +
//				                 uPriL[qvP2 ]*uConL[qBP2 ]);
//				REAL vScalBLS = (vCharM*BPar +
//				                 uPriSL[qvP1 ]*uConSL[qBP1 ] +
//				                 uPriSL[qvP2 ]*uConSL[qBP2 ]);
//				               
//
//				uConSL[q_Eges] = (vSigL*egesL - ptotalL(i)*vParL + ptotalS*vCharM +
//				                 BPar*(vScalBL - vScalBLS))/(vCharL - vCharM);
//
//
//
//				REAL egesR = uConR[q_Eges];
//				REAL vScalBR  = (uPriR[qvPar]*BPar +
//				                 uPriR[qvP1 ]*uConR[qBP1 ] +
//				                 uPriR[qvP2 ]*uConR[qBP2 ]);
//				REAL vScalBRS = (vCharM*BPar +
//				                 uPriSR[qvP1 ]*uConSR[qBP1 ] +
//				                 uPriSR[qvP2 ]*uConSR[qBP2 ]);
//				               
//
//				uConSR[q_Eges] = (vSigR*egesR - ptotalR(i)*vParR + ptotalS*vCharM +
//				                 BPar*(vScalBR - vScalBRS))/(vCharR - vCharM);
//
//
//				if(verbosity>0 && i==27 && false) {
//					cout << " Ethingy  " << uConSL[q_Eges] << " " << uConSR[q_Eges] << endl;
//				}
//
//
//				// In case of dual energy approach we need to compute
//				// the entropy at the different positions - for this
//				// we just use the fact that entropy s is given as s =
//				// p/rho^{gamma-1}
//#if(CRSWITCH_DUAL_ENERGY == CRONOS_ON)
//
//				// First compute thermal pressure from total pressure
//				// (constant in starred region -- see Eq (40)
//
//				REAL BsqSL = (sqr(uConSL[qBPar]) + sqr(uConSL[qBP1]) +
//				              sqr(uConSL[qBP2]));
//				REAL BsqSR = (sqr(uConSR[qBPar]) + sqr(uConSR[qBP1]) +
//				              sqr(uConSR[qBP2]));
//
//				REAL pthSL = ptotalS - 0.5*BsqSL;
//				REAL pthSR = ptotalS - 0.5*BsqSR;
//
//				REAL fac = 1./(gamma-1);
//				REAL EthSL = pthSL*fac;
//				REAL EthSR = pthSR*fac;
//
//				// Compute entropy:
//				uConSL[q_Eadd] = EthSL/pow(uConSL[q_rho], gamma-1.);
//				uConSR[q_Eadd] = EthSR/pow(uConSR[q_rho], gamma-1.);
////				uConSL[q_Eadd] = pthSL/pow(uConSL[q_rho], gamma-1.);
////				uConSR[q_Eadd] = pthSR/pow(uConSR[q_rho], gamma-1.);
//
//				if(verbosity==22 && (i>48 && i<52)) {
//					cout << "    Local entropy out " << i << " ";
//					cout << uConSL[q_Eadd] << " ";
//					cout << uConSR[q_Eadd] << " ";
//					cout << endl;
//				}
//
//#endif
//
//				/*-------------------------------------------------------
//				  With these the computation for the stared regions is
//				  completed. Therefore, we can now obtain the corresponding
//				  fluxes:
//				  -----------------------------------------------------*/
//
//				//if(verbosity>0 && i==27 && false) {
//				//	cout << "    relevant region ";
//				//	cout << vCharSL << " " << vCharSR << endl;
//				//}
//
//				if(vCharSL >= 0.) {
//
//					for(int q=0; q<n_omInt; ++q) {
//						fl.fluxORIG[q](i) = pfR.fluxORIG[q](i-1) + vCharL*(uConSL[q] -
//						                                           uConL[q]);
//					}
//
//					//if(verbosity==62 && i==50) {
//					//	cout << " left " << fl.fluxORIG[0](i) << " " << vCharL << " ";
//					//	cout << pfR.fluxORIG[0](i-1) << " ";
//					//	cout << uConSL[0] << " " << uConL[0] << " ";
//					//	cout << endl;
//					//}
//
//					fl.ptotalORIG(i) = ptotalL(i);
//					//if(gdata.time>0.0999 && gdata.time<0.1001 && dir==0) {
//					//	cout << i << " " << gdata.getCen_x(i) << " left star " << endl;
//					//}
//					//if(verbosity > 0 && false) {
//					//	cout << " cell " << i << " in region sl " << endl;
//					//}
//
//
//				} else if (vCharSR <= 0.) {
//
//					for(int q=0; q<n_omInt; ++q) {
//						fl.fluxORIG[q](i) = pfL.fluxORIG[q](i  ) + vCharR*(uConSR[q] -
//						                                           uConR[q]);
//					}
//					fl.ptotalORIG(i) = ptotalR(i);
//					//if(verbosity > 0 && false) {
//					//	cout << " cell " << i << " in region sr " << endl;
//					//}
//					//if(verbosity==62 && i==50) {
//					//	cout << " more left " << fl.fluxORIG[0](i) << endl;
//					//}
//
//					//if(verbosity==22 && (i>48 && i<52)) {
//					//	cout << endl << " right flux ";
//					//	cout << fl.fluxORIG[q_Eadd](i) << " ";
//					//	cout << pfL.fluxORIG[q_Eadd](i) << " ";
//					//	cout << uConSR[q_Eadd] << " ";
//					//	cout << uConR[q_Eadd] << " ";
//					//	cout << uConL[q_Eadd] << " ";
//					//	cout << endl;
//					//}
//
//
////					if(gdata.time>0.0999 && gdata.time<0.1001 && dir==0) {
////						cout << i << " " << gdata.getCen_x(i) << " right star" << endl;
////					}
//
//				} else {
//
//					/*------------------------------------------------------- I
//					  not inside the above regions we now have to compute state
//					  u**
//
//					  This is only computed if necessary. For Bx = 0 this state
//					  does not exist
//
//					  -----------------------------------------------------*/
//
//					// Define NumArray to hold conservative variables in double star region:
//					NumArray<double> uConSSL(n_omInt), uConSSR(n_omInt);
//
//					// Exclude case of zero parallel magnetic field:
//					if( (0.5*sqr(BPar) < eps*ptotalS) ) {
//
//						uConSSL = uConSL;
//						uConSSR = uConSR;
//
//					} else {
//
//						//					REAL uConSSL[n_omInt], uConSSR[n_omInt];
//						//REAL uPriSS[n_omInt];
//						std::vector<REAL> uPriSS(n_omInt);
//
//
//						// Density does not change - see Eq. (49):
//						uConSSL[q_rho] = uConSL[q_rho];
//						uConSSR[q_rho] = uConSR[q_rho];
//
//
//						REAL idenom = 1./(sqrtRhoSL + sqrtRhoSR);
//
//						// Computation for velocity (only one ** state) - see
//						// Eqs. (59) and (60)
//						uPriSS[qvP1] = (sqrtRhoSL*uPriSL[qvP1] + sqrtRhoSR*uPriSR[qvP1] +
//								(uConSR[qBP1] - uConSL[qBP1])*sign(BPar))*idenom;
//
//						uPriSS[qvP2] = (sqrtRhoSL*uPriSL[qvP2] + sqrtRhoSR*uPriSR[qvP2] +
//								(uConSR[qBP2] - uConSL[qBP2])*sign(BPar))*idenom;
//
//						// Corresponding Computation for Momentum:
//						uConSSL[qvPar] = uConSSL[q_rho]*vCharM;
//						uConSSR[qvPar] = uConSSR[q_rho]*vCharM;
//
//						uConSSL[qvP1 ] = uConSSL[q_rho]*uPriSS[qvP1];
//						uConSSR[qvP1 ] = uConSSR[q_rho]*uPriSS[qvP1];
//
//						uConSSL[qvP2 ] = uConSSL[q_rho]*uPriSS[qvP2];
//						uConSSR[qvP2 ] = uConSSR[q_rho]*uPriSS[qvP2];
//
//
//						// Computation for magnetic field (only one state **) - see
//						// Eqs. (61) and (62)
//						// *** HERE WE DELETED THE RELEVANT BUG (sign(qBPar) -> sign(BPar)) ***
//						uPriSS[qBP1] = (sqrtRhoSL*uConSR[qBP1] + sqrtRhoSR*uConSL[qBP1] +
//								sqrtRhoSL*sqrtRhoSR*(uPriSR[qvP1] -
//										uPriSL[qvP1])*sign(BPar))*idenom;
//
//						uPriSS[qBP2] = (sqrtRhoSL*uConSR[qBP2] + sqrtRhoSR*uConSL[qBP2] +
//								sqrtRhoSL*sqrtRhoSR*(uPriSR[qvP2] -
//										uPriSL[qvP2])*sign(BPar))*idenom;
//
//						//if(verbosity==22 && (i>34 && i<38)) {
//						//	cout << " qBP2 " << i << " ";
//						//	cout << idenom << " ";
//						//	cout << uPriSR[qvP2] << " ";
//						//	cout << uPriSL[qvP2] << " ";
//						//	cout << endl;
//
//						//}
//
//
//						// Set relevant values
//						uConSSL[qBPar] = uConSSR[qBPar] = BPar;
//
//						uConSSL[qBP1 ] = uConSSR[qBP1 ] = uPriSS[qBP1 ];
//
//						uConSSL[qBP2 ] = uConSSR[qBP2 ] = uPriSS[qBP2 ];
//
//						//if((verbosity==42 && i==65+1)) {
//						//	cout << " By bot ";
//						//	cout << uPriSS[qBP2] << " ";
//						//	cout << uPriSL[qBP2] << " ";
//						//	cout << uPriSR[qBP2] << " ";
//						//	cout << uPriL[qBP2] << " ";
//						//	cout << uPriR[qBP2] << " ";
//						//	cout << uPriSR[qvP2] -
//						//			uPriSL[qvP2] << " ";
//						//	cout << endl;
//						//}
//						//if((verbosity==24 && i==gdata.mx[1]-65)) {
//						//	cout << " By top ";
//						//	cout << uPriSS[qBP2] << " ";
//						//	cout << uPriSL[qBP2] << " ";
//						//	cout << uPriSR[qBP2] << " ";
//						//	cout << uPriL[qBP2] << " ";
//						//	cout << uPriR[qBP2] << " ";
//						//	cout << uPriSR[qvP2] -
//						//			uPriSL[qvP2] << " ";
//						//	cout << endl;
//						//}
//
//
//						//if(verbosity>0 && i==27 && false) {
//						//	cout << "   Deep ";
//						//	cout << uConSSL[qBPar] << " ";
//						//	cout << uConSSL[qBP1 ] << " ";
//						//	cout << uConSSL[qBP2 ] << " ";
//						//	cout << uPriSR[qvP1] << " ";
//						//	cout << uPriSL[qvP1] << " ";
//						//	cout << uPriSR[qvP2] << " ";
//						//	cout << uPriSL[qvP2] << " ";
//						//	cout << endl;
//						//}
//
//
//						// Computation for overall energy -- see Eq. (63)
//
//						REAL vScalBSS = (vCharM*BPar + uPriSS[qvP1]*uPriSS[qBP1] +
//								uPriSS[qvP2]*uPriSS[qBP2]);
//
//						uConSSL[q_Eges] = uConSL[q_Eges] - sqrtRhoSL*(vScalBLS -
//								vScalBSS)*sign(BPar);
//
//						uConSSR[q_Eges] = uConSR[q_Eges] + sqrtRhoSR*(vScalBRS -
//								vScalBSS)*sign(BPar);
//
//
//
//
//#if(CRSWITCH_DUAL_ENERGY == CRONOS_ON)
//						// Compute local entropy for dual energy approach:
//						// Get local mag field (square):
//						REAL BsqSS = (sqr(uConSSL[qBPar]) + sqr(uConSSL[qBP1 ]) +
//								sqr(uConSSL[qBP2 ]));
//						// Get local thermal pressure:
//						REAL pthSS = ptotalS - 0.5*BsqSS;
//						REAL EthSS = pthSS/(gamma-1);
//
//						// Compute resulting entropy
////						REAL entSS = pthSS/(pow(uConSSL[q_rho], gamma-1.));
//						REAL entSS = EthSS/(pow(uConSSL[q_rho], gamma-1.));
//
//						// Set values:
//						uConSSL[q_Eadd] = uConSSR[q_Eadd] = entSS;
//
//						if(verbosity==22 && (i>48 && i<52)) {
//							cout << "    Local entropy " << i << " ";
//							cout << entSS << " " << pthSS << " ";
//							cout << ptotalS << " ";
//							cout << BsqSS << " ";
//							cout << uConSSL[qBPar] << " " << uConSSL[qBP1 ] << "  "<< uConSSL[qBP2 ] << "  ";
//							cout << qBP2 << " ";
//							cout << endl;
//						}
//
//#endif
//
//
//					}
//					//if((verbosity==42 && i==65+1)) {
//					//		cout << endl << " innermost bottom " << i << " ";
//					//		cout << uConSSL[q_Eges] << " ";
//					//		cout << uConSSR[q_Eges] << " ";
//					//		cout << uConSL[q_Eges] << " ";
//					//		cout << uConSR[q_Eges] << " ";
//					//		cout << uConL[q_Eges] << " ";
//					//		cout << uConR[q_Eges] << " ";
//					//		cout << endl << "    ";
//					//		cout << ptotalL[i] << " ";
//					//		cout << ptotalR[i] << " ";
//					//		cout << ptotalS << " ";
//					//		cout << endl << "    ";
//					//		cout << vCharSL*(uConSSL[7] - uConL[7]) << " ";
//					//		cout << pfR.fluxORIG[7](i-1) + vCharSL*(uConSSL[7] - uConL[7]) << " ";
//					//		cout << pfL.fluxORIG[7](i  ) + vCharSR*(uConSSR[7] - uConR[7]) << " ";
//					//		cout << endl;
//					//	}
//					//	if((verbosity==24 && i==gdata.mx[1]-65)) {
//					//		cout << endl << " innermost top " << i << " ";
//					//		cout << uConSSR[q_Eges] << " ";
//					//		cout << uConSSL[q_Eges] << " ";
//					//		cout << uConSR[q_Eges] << " ";
//					//		cout << uConSL[q_Eges] << " ";
//					//		cout << uConR[q_Eges] << " ";
//					//		cout << uConL[q_Eges] << " ";
//					//		cout << endl << "    ";
//					//		cout << pfR.fluxORIG[7](i-1) + vCharSL*(uConSSL[7] - uConL[7]) << " ";
//					//		cout << pfL.fluxORIG[7](i  ) + vCharSR*(uConSSR[7] - uConR[7]) << " ";
//					//		cout << endl;
//
//					//	}
//
//					// Now compute the corresponding fluxes:
//					if(vCharM >= 0.) {
//						for(int q=0; q<n_omInt; ++q) {
//
//							fl.fluxORIG[q](i) = pfR.fluxORIG[q](i-1) + (vCharSL*(uConSSL[q] -
//							                                             uConSL[q]) +
//							                                    vCharL *(uConSL[q] -
//							                                             uConL[q]));
//
//						}
//						fl.ptotalORIG(i) = ptotalL(i);
//						if(verbosity > 0 && false) {
//							cout << " cell " << i << " in region ssl " << endl;
//						}
//
//						if(verbosity==22 && (i>48 && i<52)) {
//							cout <<  " hier fluss " << i << " ";
//							cout << fl.fluxORIG[q_Eadd](i) << " ";
//							cout << endl;
//						}
//						//						if(gdata.time>0.0999 && gdata.time<0.1001 && dir==0) {
////							cout << i << " " << gdata.getCen_x(i) << " left star star " << endl;
////						}
//
//					} else {
//						for(int q=0; q<n_omInt; ++q) {
//
//							fl.fluxORIG[q](i) = pfL.fluxORIG[q](i  ) + (vCharSR*(uConSSR[q] -
//							                                             uConSR[q]) +
//							                                    vCharR *(uConSR[q] -
//							                                             uConR[q]));
//
//						}
//						fl.ptotalORIG(i) = ptotalR(i);
//						if(verbosity > 0 && false) {
//							cout << " cell " << i << " in region ssr " << endl;
//						}
//						if(verbosity==22 && (i>48 && i<52)) {
//							cout <<  " hier fluss " << i << " ";
//							cout << fl.fluxORIG[q_Eadd](i) << " ";
////							cout << pfL.flux[q_Eadd](i  ) << " ";
////							cout << vCharSR << " " << vCharR << " ";
//							cout << uConSSR[q_Eadd] << " " << uConR[q_Eadd] << " ";
//							cout << uConL[q_Eadd] << " " << uConSR[q_Eadd] << " ";
//							cout << uConSL[q_Eadd] << " " << uConSSL[q_Eadd] << " ";
//							cout << endl;
//						}
////						if(gdata.time>0.0999 && gdata.time<0.1001 && dir==0) {
////							cout << i << " " << gdata.getCen_x(i) << " right star star " << endl;
////						}
//
//					}
//
//				}
//
//			}
//		}
//
//
//		if(verbosity==42) {
//			cout << " resulting flux 66 " << fl.fluxORIG[7](65+1) << endl;
//		}
//		if(verbosity==24) {
//			cout << " resulting flux oben " << fl.fluxORIG[7](gdata.mx[1]-65) << endl;
//		}
//
//
//
//	} else {
//
//		for (int i = -1; i <= gdata.mx[dir]+1; ++i){
//
//			// Derivation for isothermal case.
//
//			/*-------------------------------------------------------
//			  Compute fluxes for regions outside the Riemann fan:
//			  -----------------------------------------------------*/
//			
//			if(fl.v_ch_mORIG(i) <= 0.) { // Riemann fan going right
//				
//				for(int q=0; q<n_omInt; ++q) {
//					fl.fluxORIG[q](i) = pfR.fluxORIG[q](i-1);
//				}
//				fl.ptotalORIG(i) = pfR.ptotalORIG(i-1);
//
//			} else if (fl.v_ch_pORIG(i) <= 0.) { // Riemann fan going left
//				
//				for(int q=0; q<n_omInt; ++q) {
//					fl.fluxORIG[q](i) = pfL.fluxORIG[q](i);
//				}
//				fl.ptotalORIG(i) = pfL.ptotalORIG(i);
//				
//			} else {
//				
//			
//				/*-------------------------------------------------------
//				  The rest is taking place inside the Riemann fan
//				  -------------------------------------------------------*/
//				
//				//REAL uConL[n_omInt], uConR[n_omInt];
//				std::vector<REAL> uConL(n_omInt);
//				std::vector<REAL> uConR(n_omInt);
//				//REAL uPriL[n_omInt], uPriR[n_omInt];
//				std::vector<REAL> uPriL(n_omInt);
//				std::vector<REAL> uPriR(n_omInt);
//				//REAL uConSL[n_omInt], uConSR[n_omInt];
//				std::vector<REAL> uConSL(n_omInt);
//				std::vector<REAL> uConSR(n_omInt);
//				
//				// Saving array values:
//				
//				for(int q=0; q<n_omInt; ++q) {
//					/* 
//					   Be aware of the two differnt meanings of "L" For
//					   hlld "L" means "on the Left hand side of the cell
//					   interface" For the reconstruction it means on the
//					   left side of the cell
//					*/
//					uConR[q] = pfL.uConORIG[q](i  );
//					uConL[q] = pfR.uConORIG[q](i-1);
//					uPriR[q] = pfL.uPriORIG[q](i  );
//					uPriL[q] = pfR.uPriORIG[q](i-1);
//				}
//				
//				REAL rhoL = uConL[q_rho];
//				REAL rhoR = uConR[q_rho];
//
//				REAL vParL = uPriL[qvPar];
//				REAL vParR = uPriR[qvPar];
//				
//				REAL BPar = uConL[qBPar];
//			
//				REAL vCharL = -fl.v_ch_mORIG(i);
//				REAL vCharR =  fl.v_ch_pORIG(i);
//
//				REAL vSigL = vCharL - vParL;
//				REAL vSigR = vCharR - vParR;
//
//				REAL idenom = 1./(vCharR - vCharL);
//				REAL rhoHll = (rhoR*vSigR - rhoL*vSigL)*idenom;
//				REAL rhoS = rhoHll;
//				REAL fluxRhoHll = (vCharL*rhoR*vSigR - vCharR*rhoL*vSigL)*idenom;
//
//			
//				// See text below Eq. (23) in Mignone (2007)
//				REAL uS = fluxRhoHll/rhoHll;
//				REAL sqrtRho = sqrt(rhoHll);
//				
//				REAL vCharSL = uS - std::abs(BPar)/sqrtRho;
//				REAL vCharSR = uS + std::abs(BPar)/sqrtRho;
//
//
//				/*------------------------------------------------------
//				  Whenever a degeneracy occurrs (if vCharSL -> vCharL or
//				  vCharSR -> vCharR) we revert to hll
//				  -----------------------------------------------------*/
//
//				bool use_hll = false;
//				
//				if((vCharSL - vCharL) <  loceps*(vCharR - vCharL)) {
//					use_hll = true;
//				}
//
//				if((vCharSR - vCharR) > -loceps*(vCharR - vCharL)) {
//					use_hll = true;
//				}
//
//
//
//				idenom = 1./(vCharR - vCharL + veps);
//
//#if EXTRACT_PRESSURE == TRUE
//				fl.ptotal(i) = (vCharR*pfL.ptotal(i  ) -
//				                vCharL*pfR.ptotal(i-1))*idenom;
//#endif
//				
//
//				if(use_hll) {
//
//					for(int q=0; q<n_omInt; ++q) {
//
//						fl.fluxORIG[q](i) = (fl.v_ch_mORIG(i)*pfL.fluxORIG[q](i) +
//						                 fl.v_ch_pORIG(i)*pfR.fluxORIG[q](i-1) -
//						                 fl.v_ch_mORIG(i)*fl.v_ch_pORIG(i)*(pfL.uConORIG[q](i  )-
//						                                            pfR.uConORIG[q](i-1)))*idenom;
//					
//					}
//					
//				} else {
//
//					REAL idenom = 1./(vCharR - vCharL + veps);
//
//					/*----------------------------------------------------
//					  Set fluxes, which are constant throughout the fan:
//					  these are: density, parallel velocity and parallel
//					  magnetic field
//					  --------------------------------------------------*/
//
//					fl.fluxORIG[q_rho](i) = fluxRhoHll;
//				
//					fl.fluxORIG[qvPar](i) = (vCharR*pfR.fluxORIG[qvPar](i-1) -
//					                     vCharL*pfL.fluxORIG[qvPar](i  ) +
//					                     vCharR*vCharL*(uConR[qvPar] -
//					                                    uConL[qvPar]))*idenom;
//
//					fl.fluxORIG[qBPar](i) = 0.;
//
//					/*----------------------------------------------------
//					  Compute state u*
//					  ---------------------------------------------------*/
//
//					REAL idenomL = 1./((vCharL - vCharSL)*(vCharL - vCharSR));
//
//					// See Eq. (30) in Mignone (2007)
//					uConSL[qvP1] = rhoS*uPriL[qvP1] - 
//						BPar*uConL[qBP1]*(uS - uPriL[qvPar])*idenomL;
//					// See Eq. (31) in Mignone (2007)
//					uConSL[qvP2] = rhoS*uPriL[qvP2] - 
//						BPar*uConL[qBP2]*(uS - uPriL[qvPar])*idenomL;
//
//
//					REAL idenomR = 1./((vCharR - vCharSL)*(vCharR - vCharSR));
//
//					uConSR[qvP1] = rhoS*uPriR[qvP1] - 
//						BPar*uConR[qBP1]*(uS - uPriR[qvPar])*idenomR;
//					uConSR[qvP2] = rhoS*uPriR[qvP2] - 
//						BPar*uConR[qBP2]*(uS - uPriR[qvPar])*idenomR;
//
//
//					REAL rhoSInv = 1./rhoS;
//					// See Eq. (32) in Mignone (2007)
//					uConSL[qBP1] = uConL[qBP1]*rhoSInv*(uConL[q_rho]*sqr(vSigL) -
//					                                    sqr(BPar))*idenomL;
//					// See Eq. (33) in Mignone (2007)
//					uConSL[qBP2] = uConL[qBP2]*rhoSInv*(uConL[q_rho]*sqr(vSigL) -
//					                                    sqr(BPar))*idenomL;
//
//					uConSR[qBP1] = uConR[qBP1]*rhoSInv*(uConR[q_rho]*sqr(vSigR) -
//					                                    sqr(BPar))*idenomR;
//					uConSR[qBP2] = uConR[qBP2]*rhoSInv*(uConR[q_rho]*sqr(vSigR) -
//					                                    sqr(BPar))*idenomR;
//
//
//					// Set velocity & magnetic fluxes (see Eq. (38) in Mignone (2007))
//					if(vCharSL >= 0.) {
//				
//						fl.fluxORIG[qvP1](i) = pfR.fluxORIG[qvP1](i-1) + vCharL*(uConSL[qvP1] -
//						                                                 uConL[qvP1]);
//						fl.fluxORIG[qvP2](i) = pfR.fluxORIG[qvP2](i-1) + vCharL*(uConSL[qvP2] -
//						                                                 uConL[qvP2]);
//
//						fl.fluxORIG[qBP1](i) = pfR.fluxORIG[qBP1](i-1) + vCharL*(uConSL[qBP1] -
//						                                                 uConL[qBP1]);
//						fl.fluxORIG[qBP2](i) = pfR.fluxORIG[qBP2](i-1) + vCharL*(uConSL[qBP2] -
//						                                                 uConL[qBP2]);
//					} else if (vCharSR <= 0.) {
//
//						fl.fluxORIG[qvP1](i) = pfL.fluxORIG[qvP1](i  ) + vCharR*(uConSR[qvP1] -
//						                                                 uConR[qvP1]);
//						fl.fluxORIG[qvP2](i) = pfL.fluxORIG[qvP2](i  ) + vCharR*(uConSR[qvP2] -
//						                                                 uConR[qvP2]);
//
//						fl.fluxORIG[qBP1](i) = pfL.fluxORIG[qBP1](i  ) + vCharR*(uConSR[qBP1] -
//						                                                 uConR[qBP1]);
//						fl.fluxORIG[qBP2](i) = pfL.fluxORIG[qBP2](i  ) + vCharR*(uConSR[qBP2] -
//						                                                 uConR[qBP2]);
//					} else {
//					
//						/*----------------------------------------------------------
//						  Doing computation for central fan (only if necessary):
//						  --------------------------------------------------------*/
//
//						int sBPar = sign(BPar);
//						REAL fac = sBPar*sqrtRho;
//
//						//REAL uConSC[n_omInt];
//						std::vector<REAL> uConSC(n_omInt);
//
//						// See Eq. (34) in Mignone (2007)
//						uConSC[qvP1] = 0.5*((uConSL[qvP1] + uConSR[qvP1]) + 
//						                    (uConSR[qBP1] - uConSL[qBP1])*fac);
//
//						// See Eq. (35) in Mignone (2007)
//						uConSC[qvP2] = 0.5*((uConSL[qvP2] + uConSR[qvP2]) + 
//						                    (uConSR[qBP2] - uConSL[qBP2])*fac);
//
//						// See Eq. (36) in Mignone (2007)
//						uConSC[qBP1] = 0.5*((uConSL[qBP1] + uConSR[qBP1]) +
//						                    (uConSR[qvP1] - uConSL[qvP1])/fac);
//
//						// See Eq. (37) in Mignone (2007)
//						uConSC[qBP2] = 0.5*((uConSL[qBP2] + uConSR[qBP2]) +
//						                    (uConSR[qvP2] - uConSL[qvP2])/fac);
//
//						// Compute numerical fluxes from physical flux prescription
//						// (see Eq. (24) in Mignone (2007))
//						fl.fluxORIG[qvP1](i) = uConSC[qvP1]*uS - BPar*uConSC[qBP1];
//						fl.fluxORIG[qvP2](i) = uConSC[qvP2]*uS - BPar*uConSC[qBP2];
//
//						fl.fluxORIG[qBP1](i) = uConSC[qBP1]*uS - BPar*uConSC[qvP1]/rhoS;
//						fl.fluxORIG[qBP2](i) = uConSC[qBP2]*uS - BPar*uConSC[qvP2]/rhoS;
//
//					}
//				}
//			}
//		}
//	}
//}
//
//
//void HLLDSolver::get_NumFlux(Data &gdata, phys_fields_0D &pfM,
//		phys_fields_0D &pfP, num_fields_0D &f_num, int dir, int iFluid) const {
//	// In case of multifluid run, need to re-evaluate indices:
//#if(FLUID_TYPE == CRONOS_MULTIFLUID)
//	if(dir == 0) {
//		this->qvPar = q_sx;
//		this->qvP1  = q_sy;
//		this->qvP2  = q_sz;
//		this->qBPar = q_Bx;
//		this->qBP1  = q_By;
//		this->qBP2  = q_Bz;
//	} else if (dir == 1) {
//		this->qvP2  = q_sx;
//		this->qvPar = q_sy;
//		this->qvP1  = q_sz;
//		this->qBP2  = q_Bx;
//		this->qBPar = q_By;
//		this->qBP1  = q_Bz;
//	} else {
//		this->qvP1  = q_sx;
//		this->qvP2  = q_sy;
//		this->qvPar = q_sz;
//		this->qBP1  = q_Bx;
//		this->qBP2  = q_By;
//		this->qBPar = q_Bz;
//	}
//#endif
//	int n_omInt(pfM.get_num());
//
//	double ptotalL, ptotalR;
//
//	if(ENERGETICS == FULL) {
//
//
//		/*-------------------------------------------------------
//			  Compute fluxes for regions outside the Riemann fan:
//			  -----------------------------------------------------*/
//
//		if(f_num.v_ch_m <= 0.) { // Riemann fan going right
//
//			for(int q=0; q<n_omInt; ++q) {
//				f_num.flux_num[q] = pfM.flux_phys[q];
//			}
//			f_num.ptotal_num = pfM.ptotal;
//
//		} else if (f_num.v_ch_p <= 0.) { // Riemann fan going left
//
//			for(int q=0; q<n_omInt; ++q) {
//				f_num.flux_num[q] = pfP.flux_phys[q];
//			}
//			f_num.ptotal_num = pfP.ptotal;
//
//		} else {
//
//			/*-------------------------------------------------------
//				  The rest is taking place inside the Riemann fan
//				  -------------------------------------------------------*/
//
//			// If carbuncle test is used -> check whether cells i and i-1 are
//			// flagged as lying within the vicinity of a strong shock
//			if(gdata.use_carbuncleFlag) {
//				// if both cells are flagged -> use hll for all variables instead
//				if(pfM.carbuncle_flag*pfP.carbuncle_flag > 0) {
//					for(int q=0; q<pfM.get_num(); ++q) {
//						REAL fac = 1./(f_num.v_ch_p+f_num.v_ch_m+veps);
//
//						f_num.flux_num[q] = (f_num.v_ch_m*pfP.flux_phys[q] +
//								f_num.v_ch_p*pfM.flux_phys[q] -
//								f_num.v_ch_m*f_num.v_ch_p*(pfP.uCon[q] - pfM.uCon[q]))*fac;
//					}
//					return;
//				}
//			}
//
//
//			//REAL uConL[n_omInt], uConR[n_omInt];
//			std::vector<REAL> uConL(n_omInt);
//			std::vector<REAL> uConR(n_omInt);
//			//REAL uPriL[n_omInt], uPriR[n_omInt];
//			std::vector<REAL> uPriL(n_omInt);
//			std::vector<REAL> uPriR(n_omInt);
//			//REAL uConSL[n_omInt], uPriSL[n_omInt];
//			std::vector<REAL> uConSL(n_omInt);
//			std::vector<REAL> uPriSL(n_omInt);
//			//REAL uConSR[n_omInt], uPriSR[n_omInt];
//			std::vector<REAL> uConSR(n_omInt);
//			std::vector<REAL> uPriSR(n_omInt);
//
//			// Saving array values:
//
//			for(int q=0; q<n_omInt; ++q) {
//				/*
//					In this context "L" and "M" mean (L)eft from the cell-face
//					or on the (M)inus side
//					The same holds for "R" and "P"
//				 */
//				uConR[q] = pfP.uCon[q];
//				uConL[q] = pfM.uCon[q];
//				uPriR[q] = pfP.uPri[q];
//				uPriL[q] = pfM.uPri[q];
//			}
//
//#if (USE_COROTATION == CRONOS_ON)
//			// Overwrite settings for co-rotation case:
//			// uConR[qvPar] = uConL[q_rho]*uPriL[qvPar];
//			// uConR[qvP1]  = uConL[q_rho]*uPriL[qvP1];
//			// uConR[qvP2]  = uConL[q_rho]*uPriL[qvP2];
//			uConR[qvPar] = uConR[q_rho]*uPriR[qvPar];
//			uConR[qvP1]  = uConR[q_rho]*uPriR[qvP1];
//			uConR[qvP2]  = uConR[q_rho]*uPriR[qvP2];
//			uConL[qvPar] = uConL[q_rho]*uPriL[qvPar];
//			uConL[qvP1]  = uConL[q_rho]*uPriL[qvP1];
//			uConL[qvP2]  = uConL[q_rho]*uPriL[qvP2];
//#endif
//
//			REAL rhoL = uConL[q_rho];
//			REAL rhoR = uConR[q_rho];
//
//			REAL vParL = uPriL[qvPar];
//			REAL vParR = uPriR[qvPar];
//
//			REAL sParL = uConL[qvPar];
//			REAL sParR = uConR[qvPar];
//
//			REAL BPar = uConL[qBPar];
//
//			REAL pBL = sqr(uConL[qBPar]) + sqr(uConL[qBP1 ]) + sqr(uConL[qBP2 ]);
//			REAL pBR = sqr(uConR[qBPar]) + sqr(uConR[qBP1 ]) + sqr(uConR[qBP2 ]);
//
//			ptotalL = pfM.ptherm + 0.5*pBL;
//			ptotalR = pfP.ptherm + 0.5*pBR;
//
//			REAL vCharL = -f_num.v_ch_m;
//			REAL vCharR =  f_num.v_ch_p;
//
//			REAL vSigL = vCharL - vParL;
//			REAL vSigR = vCharR - vParR;
//
//
//			REAL idenom = 1./(vSigR*rhoR - vSigL*rhoL + veps);
//
//
//			// Eq. (38)
//			// S_M according to eq 38 in hlld paper
//			REAL vCharM = (vSigR*sParR - vSigL*sParL -
//					ptotalR + ptotalL)*idenom;
//
//
//			/*----------------------------------------------------------
//				  Start computation of plasma state in Regions L* and R*
//				  --------------------------------------------------------*/
//
//			// Total pressure in Regions L* and R* - see Eq. (41)
//			REAL ptotalS = (vSigR*rhoR*ptotalL - vSigL*rhoL*ptotalR +
//					rhoL*rhoR*vSigL*vSigR*(vParR - vParL))*idenom;
//
//
//			// Mass denisty in Regions L* and R* - Eq. (43)
//			uConSL[q_rho] = uConL[q_rho]*vSigL/(vCharL - vCharM);
//			uConSR[q_rho] = uConR[q_rho]*vSigR/(vCharR - vCharM);
//
//			// Needed lateron
//			REAL sqrtRhoSL = sqrt(uConSL[q_rho]);
//			REAL sqrtRhoSR = sqrt(uConSR[q_rho]);
//
//			// Characteristic velocities of rim of regions L* and R* - see Eq. (51)
//			REAL vCharSL = vCharM - std::abs(BPar)/sqrtRhoSL;
//			REAL vCharSR = vCharM + std::abs(BPar)/sqrtRhoSR;
//
//
//			/*
//
//				  According to Mignone a destinction has to be made in the
//				  following cases:
//
//				  Here we just compute intermediate values for the *-values of
//				  the individual components of the magnetic induction
//
//			 */
//
//			bool use_hllc = false;
//			if((vCharSL - vCharL) <  loceps*(vCharM - vCharL)) {
//				use_hllc = true;
//			}
//
//			if((vCharSR - vCharR) > -loceps*(vCharR - vCharM)) {
//				use_hllc = true;
//			}
//
//
//			if(use_hllc) {
//
//				idenom = 1./(vCharR - vCharL);
//
//				for(int q=qBMin; q<=qBMax; ++q) {
//
//					REAL value = (vCharR*uConR[q] - vCharL*uConL[q] +
//							pfM.flux_phys[q] - pfP.flux_phys[q])*idenom;
//					uConSL[q] = value;
//					uConSR[q] = value;
//
//				}
//
//				// By this a computation of Regions ** is avoided
//				vCharSL = vCharSR = vCharM;
//
//			} else {
//
//				// Using hlld
//
//				REAL facL = (rhoL*sqr(vSigL) - sqr(BPar)) /
//						(rhoL*vSigL*(vCharL - vCharM) - sqr(BPar));
//
//
//				// Compute components of mag induction
//				// see also Eqs. (45) and (47)
//				uConSL[qBPar] = BPar;
//
//				if(std::abs(vCharM) < 1.e-120 && std::abs(vCharL - vSigL) < 1.e-120) {
//					uConSL[qBP1 ] = uConL[qBP1 ];
//					uConSL[qBP2 ] = uConL[qBP2 ];
//				} else {
//					uConSL[qBP1 ] = uConL[qBP1 ]*facL;
//					uConSL[qBP2 ] = uConL[qBP2 ]*facL;
//				}
//
//
//				REAL facR = (rhoR*sqr(vSigR) - sqr(BPar)) /
//						(rhoR*vSigR*(vCharR - vCharM) - sqr(BPar));
//
//				uConSR[qBPar] = BPar;
//
//				if(std::abs(vCharM) < 1.e-120 && std::abs(vCharR - vSigR) < 1.e-120) {
//					uConSR[qBP1 ] = uConR[qBP1 ];
//					uConSR[qBP2 ] = uConR[qBP2 ];
//				} else {
//					uConSR[qBP1 ] = uConR[qBP1 ]*facR;
//					uConSR[qBP2 ] = uConR[qBP2 ]*facR;
//				}
//
//			}
//
//
//
//			// Compute perp components of velocity in regions L* and R*
//			// (this is simplified by using the above) see Eqs. (44) &
//			// (46)
//
//			REAL facL = BPar/(rhoL*vSigL);
//			REAL facR = BPar/(rhoR*vSigR);
//
//			uPriSL[qvP1] = uPriL[qvP1] - (uConSL[qBP1] - uConL[qBP1])*facL;
//			uPriSR[qvP1] = uPriR[qvP1] - (uConSR[qBP1] - uConR[qBP1])*facR;
//
//			uPriSL[qvP2] = uPriL[qvP2] - (uConSL[qBP2] - uConL[qBP2])*facL;
//			uPriSR[qvP2] = uPriR[qvP2] - (uConSR[qBP2] - uConR[qBP2])*facR;
//
//
//			// Get corresponding values for the momentum:
//			uConSL[qvPar] = uConSL[q_rho]*vCharM;
//			uConSL[qvP1 ] = uConSL[q_rho]*uPriSL[qvP1];
//			uConSL[qvP2 ] = uConSL[q_rho]*uPriSL[qvP2];
//
//			uConSR[qvPar] = uConSR[q_rho]*vCharM;
//			uConSR[qvP1 ] = uConSR[q_rho]*uPriSR[qvP1];
//			uConSR[qvP2 ] = uConSR[q_rho]*uPriSR[qvP2];
//
//
//			// Compute Energy in stared state (see Eq. (48))
//			REAL egesL = uConL[q_Eges];
//			REAL vScalBL  = (uPriL[qvPar]*BPar +
//					uPriL[qvP1 ]*uConL[qBP1 ] +
//					uPriL[qvP2 ]*uConL[qBP2 ]);
//			REAL vScalBLS = (vCharM*BPar +
//					uPriSL[qvP1 ]*uConSL[qBP1 ] +
//					uPriSL[qvP2 ]*uConSL[qBP2 ]);
//
//
//			uConSL[q_Eges] = (vSigL*egesL - ptotalL*vParL + ptotalS*vCharM +
//					BPar*(vScalBL - vScalBLS))/(vCharL - vCharM);
//
//
//			REAL egesR = uConR[q_Eges];
//			REAL vScalBR  = (uPriR[qvPar]*BPar +
//					uPriR[qvP1 ]*uConR[qBP1 ] +
//					uPriR[qvP2 ]*uConR[qBP2 ]);
//			REAL vScalBRS = (vCharM*BPar +
//					uPriSR[qvP1 ]*uConSR[qBP1 ] +
//					uPriSR[qvP2 ]*uConSR[qBP2 ]);
//
//
//			uConSR[q_Eges] = (vSigR*egesR - ptotalR*vParR + ptotalS*vCharM +
//					BPar*(vScalBR - vScalBRS))/(vCharR - vCharM);
//
//				// In case of dual energy approach we need to compute
//				// the entropy at the different positions - for this
//				// we just use the fact that entropy s is given as s =
//				// p/rho^{gamma-1}
//#if(CRSWITCH_DUAL_ENERGY == CRONOS_ON)
//
//			// First compute thermal pressure from total pressure
//			// (constant in starred region -- see Eq (40)
//
//			REAL BsqSL = (sqr(uConSL[qBPar]) + sqr(uConSL[qBP1]) +
//					sqr(uConSL[qBP2]));
//			REAL BsqSR = (sqr(uConSR[qBPar]) + sqr(uConSR[qBP1]) +
//					sqr(uConSR[qBP2]));
//
//			REAL pthSL = ptotalS - 0.5*BsqSL;
//			REAL pthSR = ptotalS - 0.5*BsqSR;
//
//			// Compute entropy:
//			uConSL[q_Eadd] = pthSL/pow(uConSL[q_rho], gamma-1.);
//			uConSR[q_Eadd] = pthSR/pow(uConSR[q_rho], gamma-1.);
//
//#endif
//
//			/*-------------------------------------------------------
//				  With these the computation for the stared regions is
//				  completed. Therefore, we can now obtain the corresponding
//				  fluxes:
//				  -----------------------------------------------------*/
//
//			if(vCharSL >= 0.) {
//
//				for(int q=0; q<n_omInt; ++q) {
//					f_num.flux_num[q] = pfM.flux_phys[q] + vCharL*(uConSL[q] -
//							uConL[q]);
//				}
//
//				f_num.ptotal_num = ptotalL;
//
//			} else if (vCharSR <= 0.) {
//
//				for(int q=0; q<n_omInt; ++q) {
//					f_num.flux_num[q] = pfP.flux_phys[q] + vCharR*(uConSR[q] -
//							uConR[q]);
//				}
//				f_num.ptotal_num = ptotalR;
//
//			} else {
//
//				/*------------------------------------------------------- I
//					  not inside the above regions we now have to compute state
//					  u**
//
//					  This is only computed if necessary. For Bx = 0 this state
//					  does not exist
//
//					  -----------------------------------------------------*/
//
//				//REAL uConSSL[n_omInt], uConSSR[n_omInt];
//				std::vector<REAL> uConSSL(n_omInt);
//				std::vector<REAL> uConSSR(n_omInt);
//				//REAL uPriSS[n_omInt];
//				std::vector<REAL> uPriSS(n_omInt);
//
//
//				// Density does not change - see Eq. (49):
//				uConSSL[q_rho] = uConSL[q_rho];
//				uConSSR[q_rho] = uConSR[q_rho];
//
//
//				REAL idenom = 1./(sqrtRhoSL + sqrtRhoSR);
//
//				// Computation for velocity (only one ** state) - see
//				// Eqs. (59) and (60)
//				uPriSS[qvP1] = (sqrtRhoSL*uPriSL[qvP1] + sqrtRhoSR*uPriSR[qvP1] +
//						(uConSR[qBP1] - uConSL[qBP1])*sign(BPar))*idenom;
//
//				uPriSS[qvP2] = (sqrtRhoSL*uPriSL[qvP2] + sqrtRhoSR*uPriSR[qvP2] +
//						(uConSR[qBP2] - uConSL[qBP2])*sign(BPar))*idenom;
//
//				// Corresponding Computation for Momentum:
//				uConSSL[qvPar] = uConSSL[q_rho]*vCharM;
//				uConSSR[qvPar] = uConSSR[q_rho]*vCharM;
//
//				uConSSL[qvP1 ] = uConSSL[q_rho]*uPriSS[qvP1];
//				uConSSR[qvP1 ] = uConSSR[q_rho]*uPriSS[qvP1];
//
//				uConSSL[qvP2 ] = uConSSL[q_rho]*uPriSS[qvP2];
//				uConSSR[qvP2 ] = uConSSR[q_rho]*uPriSS[qvP2];
//
//
//				// Computation for magnetic field (only one state **) - see
//				// Eqs. (61) and (62)
//				uPriSS[qBP1] = (sqrtRhoSL*uConSR[qBP1] + sqrtRhoSR*uConSL[qBP1] +
//						sqrtRhoSL*sqrtRhoSR*(uPriSR[qvP1] -
//								uPriSL[qvP1])*sign(BPar))*idenom;
//
//				uPriSS[qBP2] = (sqrtRhoSL*uConSR[qBP2] + sqrtRhoSR*uConSL[qBP2] +
//						sqrtRhoSL*sqrtRhoSR*(uPriSR[qvP2] -
//								uPriSL[qvP2])*sign(BPar))*idenom;
//
//
//				// Set relevant values
//				uConSSL[qBPar] = uConSSR[qBPar] = BPar;
//
//				uConSSL[qBP1 ] = uConSSR[qBP1 ] = uPriSS[qBP1 ];
//
//				uConSSL[qBP2 ] = uConSSR[qBP2 ] = uPriSS[qBP2 ];
//
//
//
//				// Computation for overall energy -- see Eq. (63)
//
//				REAL vScalBSS = (vCharM*BPar + uPriSS[qvP1]*uPriSS[qBP1] +
//						uPriSS[qvP2]*uPriSS[qBP2]);
//
//				uConSSL[q_Eges] = uConSL[q_Eges] - sqrtRhoSL*(vScalBLS -
//						vScalBSS)*sign(BPar);
//
//				uConSSR[q_Eges] = uConSR[q_Eges] + sqrtRhoSR*(vScalBRS -
//						vScalBSS)*sign(BPar);
//
//#if(CRSWITCH_DUAL_ENERGY == CRONOS_ON)
//				// Compute local entropy for dual energy approach:
//				// Get local mag field (square):
//				REAL BsqSS = (sqr(uConSSL[qBPar]) + sqr(uConSSL[qBP1 ]) +
//						sqr(uConSSL[qBP2 ]));
//				// Get local thermal pressure:
//				REAL pthSS = ptotalS - 0.5*BsqSS;
//
//				// Compute resulting entropy
//				REAL entSS = pthSS/(pow(uConSSL[q_rho], gamma-1.));
//
//				// Set values:
//				uConSSL[q_Eadd] = uConSSR[q_Eadd] = entSS;
//
//#endif
//
//
//				// Now compute the corresponding fluxes:
//				if(vCharM >= 0.) {
//					for(int q=0; q<n_omInt; ++q) {
//
//						f_num.flux_num[q] = pfM.flux_phys[q] + (vCharSL*(uConSSL[q] -
//								uConSL[q]) +
//								vCharL *(uConSL[q] -
//										uConL[q]));
//
//					}
//					f_num.ptotal_num = ptotalL;
//
//				} else {
//					for(int q=0; q<n_omInt; ++q) {
//
//						f_num.flux_num[q] = pfP.flux_phys[q] + (vCharSR*(uConSSR[q] -
//								uConSR[q]) +
//								vCharR *(uConSR[q] -
//										uConR[q]));
//
//					}
//					f_num.ptotal_num = ptotalR;
//
//				}
//
//			}
//
//		}
//
//
//	} else {
//
//
//		// Derivation for isothermal case.
//
//		/*-------------------------------------------------------
//			  Compute fluxes for regions outside the Riemann fan:
//			  -----------------------------------------------------*/
//
//		if(f_num.v_ch_m <= 0.) { // Riemann fan going right
//
//			for(int q=0; q<n_omInt; ++q) {
//				f_num.flux_num[q] = pfM.flux_phys[q];
//			}
//			f_num.ptotal_num = pfM.ptotal;
//
//		} else if (f_num.v_ch_p <= 0.) { // Riemann fan going left
//
//			for(int q=0; q<n_omInt; ++q) {
//				f_num.flux_num[q] = pfP.flux_phys[q];
//			}
//			f_num.ptotal_num = pfP.ptotal;
//
//		} else {
//
//
//			/*-------------------------------------------------------
//				  The rest is taking place inside the Riemann fan
//				  -------------------------------------------------------*/
//
//			//REAL uConL[n_omInt], uConR[n_omInt];
//			std::vector<REAL> uConL(n_omInt);
//			std::vector<REAL> uConR(n_omInt);
//			//REAL uPriL[n_omInt], uPriR[n_omInt];
//			std::vector<REAL> uPriL(n_omInt);
//			std::vector<REAL> uPriR(n_omInt);
//			//REAL uConSL[n_omInt], uConSR[n_omInt];
//			std::vector<REAL> uConSL(n_omInt);
//			std::vector<REAL> uConSR(n_omInt);
//
//			// Saving array values:
//
//			for(int q=0; q<n_omInt; ++q) {
//				/*
//									In this context "L" and "M" mean (L)eft from the cell-face
//									or on the (M)inus side
//									The same holds for "R" and "P"
//				 */
//				uConR[q] = pfP.uCon[q];
//				uConL[q] = pfM.uCon[q];
//				uPriR[q] = pfP.uPri[q];
//				uPriL[q] = pfM.uPri[q];
//			}
//
//			REAL rhoL = uConL[q_rho];
//			REAL rhoR = uConR[q_rho];
//
//			REAL vParL = uPriL[qvPar];
//			REAL vParR = uPriR[qvPar];
//
//			REAL BPar = uConL[qBPar];
//
//			REAL vCharL = -f_num.v_ch_m;
//			REAL vCharR =  f_num.v_ch_p;
//
//			REAL vSigL = vCharL - vParL;
//			REAL vSigR = vCharR - vParR;
//
//			REAL idenom = 1./(vCharR - vCharL);
//			REAL rhoHll = (rhoR*vSigR - rhoL*vSigL)*idenom;
//			REAL rhoS = rhoHll;
//			REAL fluxRhoHll = (vCharL*rhoR*vSigR - vCharR*rhoL*vSigL)*idenom;
//
//
//			// See text below Eq. (23) in Mignone (2007)
//			REAL uS = fluxRhoHll/rhoHll;
//			REAL sqrtRho = sqrt(rhoHll);
//
//			REAL vCharSL = uS - std::abs(BPar)/sqrtRho;
//			REAL vCharSR = uS + std::abs(BPar)/sqrtRho;
//
//
//			/*------------------------------------------------------
//				  Whenever a degeneracy occurrs (if vCharSL -> vCharL or
//				  vCharSR -> vCharR) we revert to hll
//				  -----------------------------------------------------*/
//
//			bool use_hll = false;
//
//			if((vCharSL - vCharL) <  loceps*(vCharR - vCharL)) {
//				use_hll = true;
//			}
//
//			if((vCharSR - vCharR) > -loceps*(vCharR - vCharL)) {
//				use_hll = true;
//			}
//
//
//
//			idenom = 1./(vCharR - vCharL + veps);
//
//#if EXTRACT_PRESSURE == TRUE
//			f_num.ptotal_num  = (vCharR*pfP.ptotal -
//					vCharL*pfM.ptotal)*idenom;
//#endif
//
//
//			if(use_hll) {
//
//				for(int q=0; q<n_omInt; ++q) {
//
//					f_num.flux_num[q] = (f_num.v_ch_m*pfP.flux_phys[q] +
//							f_num.v_ch_p*pfM.flux_phys[q] -
//							f_num.v_ch_m*f_num.v_ch_p*(pfP.uCon[q] - pfM.uCon[q]))*idenom;
//
//				}
//
//			} else {
//
//				REAL idenom = 1./(vCharR - vCharL + veps);
//
//				/*----------------------------------------------------
//					  Set fluxes, which are constant throughout the fan:
//					  these are: density, parallel velocity and parallel
//					  magnetic field
//					  --------------------------------------------------*/
//
//				f_num.flux_num[q_rho] = fluxRhoHll;
//
//				f_num.flux_num[qvPar] = (vCharR*pfM.flux_phys[qvPar] -
//						vCharL*pfP.flux_phys[qvPar] +
//						vCharR*vCharL*(uConR[qvPar] - uConL[qvPar]))*idenom;
//
//				f_num.flux_num[qBPar] = 0.;
//
//				/*----------------------------------------------------
//					  Compute state u*
//					  ---------------------------------------------------*/
//
//				REAL idenomL = 1./((vCharL - vCharSL)*(vCharL - vCharSR));
//
//				// See Eq. (30) in Mignone (2007)
//				uConSL[qvP1] = rhoS*uPriL[qvP1] -
//						BPar*uConL[qBP1]*(uS - uPriL[qvPar])*idenomL;
//				// See Eq. (31) in Mignone (2007)
//				uConSL[qvP2] = rhoS*uPriL[qvP2] -
//						BPar*uConL[qBP2]*(uS - uPriL[qvPar])*idenomL;
//
//
//				REAL idenomR = 1./((vCharR - vCharSL)*(vCharR - vCharSR));
//
//				uConSR[qvP1] = rhoS*uPriR[qvP1] -
//						BPar*uConR[qBP1]*(uS - uPriR[qvPar])*idenomR;
//				uConSR[qvP2] = rhoS*uPriR[qvP2] -
//						BPar*uConR[qBP2]*(uS - uPriR[qvPar])*idenomR;
//
//
//				REAL rhoSInv = 1./rhoS;
//				// See Eq. (32) in Mignone (2007)
//				uConSL[qBP1] = uConL[qBP1]*rhoSInv*(uConL[q_rho]*sqr(vSigL) -
//						sqr(BPar))*idenomL;
//				// See Eq. (33) in Mignone (2007)
//				uConSL[qBP2] = uConL[qBP2]*rhoSInv*(uConL[q_rho]*sqr(vSigL) -
//						sqr(BPar))*idenomL;
//
//				uConSR[qBP1] = uConR[qBP1]*rhoSInv*(uConR[q_rho]*sqr(vSigR) -
//						sqr(BPar))*idenomR;
//				uConSR[qBP2] = uConR[qBP2]*rhoSInv*(uConR[q_rho]*sqr(vSigR) -
//						sqr(BPar))*idenomR;
//
//
//				// Set velocity & magnetic fluxes (see Eq. (38) in Mignone (2007))
//				if(vCharSL >= 0.) {
//
//					f_num.flux_num[qvP1] = pfM.flux_phys[qvP1] + vCharL*(uConSL[qvP1] -
//							uConL[qvP1]);
//					f_num.flux_num[qvP2] = pfM.flux_phys[qvP2] + vCharL*(uConSL[qvP2] -
//							uConL[qvP2]);
//
//					f_num.flux_num[qBP1] = pfM.flux_phys[qBP1] + vCharL*(uConSL[qBP1] -
//							uConL[qBP1]);
//					f_num.flux_num[qBP2] = pfM.flux_phys[qBP2] + vCharL*(uConSL[qBP2] -
//							uConL[qBP2]);
//				} else if (vCharSR <= 0.) {
//
//					f_num.flux_num[qvP1] = pfP.flux_phys[qvP1] + vCharR*(uConSR[qvP1] -
//							uConR[qvP1]);
//					f_num.flux_num[qvP2] = pfP.flux_phys[qvP2] + vCharR*(uConSR[qvP2] -
//							uConR[qvP2]);
//
//					f_num.flux_num[qBP1] = pfP.flux_phys[qBP1] + vCharR*(uConSR[qBP1] -
//							uConR[qBP1]);
//					f_num.flux_num[qBP2] = pfP.flux_phys[qBP2] + vCharR*(uConSR[qBP2] -
//							uConR[qBP2]);
//				} else {
//
//					/*----------------------------------------------------------
//						  Doing computation for central fan (only if necessary):
//						  --------------------------------------------------------*/
//
//					int sBPar = sign(BPar);
//					REAL fac = sBPar*sqrtRho;
//
//					//REAL uConSC[n_omInt];
//					std::vector<REAL> uConSC(n_omInt);
//
//					// See Eq. (34) in Mignone (2007)
//					uConSC[qvP1] = 0.5*((uConSL[qvP1] + uConSR[qvP1]) +
//							(uConSR[qBP1] - uConSL[qBP1])*fac);
//
//					// See Eq. (35) in Mignone (2007)
//					uConSC[qvP2] = 0.5*((uConSL[qvP2] + uConSR[qvP2]) +
//							(uConSR[qBP2] - uConSL[qBP2])*fac);
//
//					// See Eq. (36) in Mignone (2007)
//					uConSC[qBP1] = 0.5*((uConSL[qBP1] + uConSR[qBP1]) +
//							(uConSR[qvP1] - uConSL[qvP1])/fac);
//
//					// See Eq. (37) in Mignone (2007)
//					uConSC[qBP2] = 0.5*((uConSL[qBP2] + uConSR[qBP2]) +
//							(uConSR[qvP2] - uConSL[qvP2])/fac);
//
//					// Compute numerical fluxes from physical flux prescription
//					// (see Eq. (24) in Mignone (2007))
//					f_num.flux_num[qvP1] = uConSC[qvP1]*uS - BPar*uConSC[qBP1];
//					f_num.flux_num[qvP2] = uConSC[qvP2]*uS - BPar*uConSC[qBP2];
//
//					f_num.flux_num[qBP1] = uConSC[qBP1]*uS - BPar*uConSC[qvP1]/rhoS;
//					f_num.flux_num[qBP2] = uConSC[qBP2]*uS - BPar*uConSC[qvP2]/rhoS;
//
//				}
//			}
//		}
//	}
//}


