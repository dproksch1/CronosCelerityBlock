#include "RiemannSolverHD.H"

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

/*-------------------------------------------------------------------
  HLLC Rieman solver according to the book "Riemann Solvers and
  Numerical Methods for Fluid Dynamics" by Toro (3rd. edition from
  2009)

  Specifics of the scheme are given on pages 331, 332
  -----------------------------------------------------------------*/

void HLLCSolver_Hydro::get_NumFlux(Queue queue, Data &gdata,
                             phys_fields_1D &pfL,
                             phys_fields_1D &pfR,
                             fields_1D &fl, const int &dir, int iFluid) const {

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

	int n_omInt(fl.get_num());

	//assert_sycl_eq(fl.v_ch_mORIG, fl.v_ch_mSYCL);
	//assert_sycl_eq(pfL.ptotalORIG, pfL.ptotalSYCL);
	//assert_sycl_eq(pfR.ptotalORIG, pfR.ptotalSYCL);
	//for (int q = 0; q < n_omInt; ++q) {
	//	assert_sycl_eq(pfL.fluxORIG[q], pfL.fluxSYCL[q]);
	//	assert_sycl_eq(pfR.fluxORIG[q], pfR.fluxSYCL[q]);
	//}

	//for (int q = 0; q < n_omInt; ++q) {
	for (int i = -1; i <= gdata.mx[dir]+1; ++i){

		/*-------------------------------------------------------
		  Compute fluxes for regions outside the Riemann fan:
		  -----------------------------------------------------*/
		
		if(fl.v_ch_mORIG(i) <= 0.) { // Riemann fan going right
			
			for(int q=0; q<n_omInt; ++q) {
				fl.fluxORIG[q](i) = pfR.fluxORIG[q](i-1);
			}
//#if(FLUID_TYPE==CRONOS_MULTIFLUID)
//			fl.ptotal[iFluid](i) = pfR.ptotal(i-1);
//#else
			fl.ptotalORIG(i) = pfR.ptotalORIG(i-1);
//#endif
			
		} else if (fl.v_ch_pORIG(i) <= 0.) { // Riemann fan going left
			
			for(int q=0; q<n_omInt; ++q) {
				fl.fluxORIG[q](i) = pfL.fluxORIG[q](i);
			}
//#if(FLUID_TYPE==CRONOS_MULTIFLUID)
//			fl.ptotal[iFluid](i) = pfL.ptotal(i);
//#else
			fl.ptotalORIG(i) = pfL.ptotalORIG(i);
//#endif
			
		} else {
		
			/*-------------------------------------------------------
			  The rest is taking place inside the Riemann fan
			  -------------------------------------------------------*/

			// If carbuncle test is used -> check whether cells i and i-1 are
			// flagged as lying within the vicinity of a strong shock
			if(gdata.use_carbuncleFlag) {
				// if both cells are flagged -> use hll for all variables instead
				if(fl.carbuncle_flag(i)*fl.carbuncle_flag(i-1) > 0) {
					for(int q=0; q<fl.get_num(); ++q) {
						REAL fac = 1./(fl.v_ch_pORIG(i)+fl.v_ch_mORIG(i)+veps);

						fl.fluxORIG[q](i) = (fl.v_ch_mORIG(i)*pfL.fluxORIG[q](i) +
								fl.v_ch_pORIG(i)*pfR.fluxORIG[q](i-1) -
								fl.v_ch_mORIG(i)*fl.v_ch_pORIG(i)*(pfL.uConORIG[q](i  )-
										pfR.uConORIG[q](i-1)))*fac;
					}
					continue;
				}
			}


			//REAL uConL[n_omInt], uConR[n_omInt];
			std::vector<REAL> uConL(n_omInt);
			std::vector<REAL> uConR(n_omInt);
			//REAL uPriL[n_omInt], uPriR[n_omInt];
			std::vector<REAL> uPriL(n_omInt);
			std::vector<REAL> uPriR(n_omInt);
//			REAL uConSL[n_omInt], uPriSL[n_omInt];
			//REAL uConSL[n_omInt];
			std::vector<REAL> uConSL(n_omInt);
//			REAL uConSR[n_omInt], uPriSR[n_omInt];
			//REAL uConSR[n_omInt];
			std::vector<REAL> uConSR(n_omInt);

			// Saving array values:
				
			for(int q=0; q<n_omInt; ++q) {
				/* 
				   Be aware of the two differnt meanings of "L" For
				   hllc "L" means "on the Left hand side of the cell
				   interface" For the reconstruction it means on the
				   left side of the cell
				*/
				uConR[q] = pfL.uConORIG[q](i  );
				uConL[q] = pfR.uConORIG[q](i-1);
				uPriR[q] = pfL.uPriORIG[q](i  );
				uPriL[q] = pfR.uPriORIG[q](i-1);
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

			REAL pThermL(pfR.pthermORIG(i-1));
			REAL pThermR(pfL.pthermORIG(i  ));

			REAL vParL = uPriL[qvPar];
			REAL vParR = uPriR[qvPar];

			REAL sParL = uConL[qvPar];
			REAL sParR = uConR[qvPar];
			
			// Characteristic velocities like in Toro book
			REAL vCharL = -fl.v_ch_mORIG(i);
			REAL vCharR =  fl.v_ch_pORIG(i);

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
					fl.fluxORIG[q](i) = pfR.fluxORIG[q](i-1) + vCharL*(uConSL[q] -
					                                           uConL[q]);
				}
				
				// fl.ptotal(i) = ptotalL(i);
				fl.ptotalORIG(i) = pfR.pthermORIG(i-1);
				
			} else {
				
				for(int q=0; q<n_omInt; ++q) {
					fl.fluxORIG[q](i) = pfL.fluxORIG[q](i  ) + vCharR*(uConSR[q] -
					                                           uConR[q]);
				}
				// fl.ptotal(i) = ptotalR(i);
				fl.ptotalORIG(i) = pfL.pthermORIG(i  );
				
			}
	
		}
	}

	// SYCL
#ifdef USE_SYCL
	size_t size = gdata.mx[dir];

	for (int q = 0; q < n_omInt; ++q) {
		queue.submit([&](sycl::handler& cgh) {
			auto fluxL_acc = pfL.fluxSYCL[q].get_access<cl::sycl::access::mode::read>(cgh);
			auto fluxR_acc = pfR.fluxSYCL[q].get_access<cl::sycl::access::mode::read>(cgh);
			auto fluxF_acc = fl.fluxSYCL[q].get_access<cl::sycl::access::mode::write>(cgh);
			auto v_ch_m_acc = fl.v_ch_mSYCL.get_access<cl::sycl::access::mode::read>(cgh);
			auto v_ch_p_acc = fl.v_ch_pSYCL.get_access<cl::sycl::access::mode::read>(cgh);
			auto ptotalL_acc = pfL.ptotalSYCL.get_access<cl::sycl::access::mode::read>(cgh);
			auto ptotalR_acc = pfR.ptotalSYCL.get_access<cl::sycl::access::mode::read>(cgh);
			auto ptotalF_acc = fl.ptotalSYCL.get_access<cl::sycl::access::mode::write>(cgh);
			cgh.parallel_for<class FluxOutsideRiemann>(Range<1>(size + 3), Id<1>(1), [=](Item<1> item) {
				size_t i = item.get(0);
				/*-------------------------------------------------------
				  Compute fluxes for regions outside the Riemann fan:
				  -----------------------------------------------------*/

				if (v_ch_m_acc[i] <= 0.0) { // Riemann fan going right

					//for (int q = 0; q < n_omInt; ++q) {
						fluxF_acc[i] = fluxR_acc[i - 1];
					//}
					ptotalF_acc[i] = ptotalR_acc[i - 1];

				} else if (v_ch_p_acc[i] <= 0.0) { // Riemann fan going left

					//for (int q = 0; q < n_omInt; ++q) {
						fluxF_acc[i] = fluxL_acc[i];
					//}
					ptotalF_acc[i] = ptotalL_acc[i];

				} else {
					// else branch will be handled in a follow-up kernel
				}
			});
		});
	}

	using BufferType = Buffer<REAL, 2>;
	int buffSize = gdata.mx[dir] + 1 + 2 + 1;

	BufferType uConR = BufferType(Range<2>(buffSize, n_omInt));
	BufferType uConL = BufferType(Range<2>(buffSize, n_omInt));
	BufferType uConSR = BufferType(Range<2>(buffSize, n_omInt));
	BufferType uConSL = BufferType(Range<2>(buffSize, n_omInt));
	BufferType uPriR = BufferType(Range<2>(buffSize, n_omInt));
	BufferType uPriL = BufferType(Range<2>(buffSize, n_omInt));
	Buffer<REAL, 1> vCharS = Buffer<REAL, 1>(Range<1>(buffSize));

	for (int q = 0; q < n_omInt; ++q) {
		queue.submit([&](sycl::handler& cgh) {
			auto uConL_acc = pfL.uConSYCL[q].get_access<cl::sycl::access::mode::read>(cgh);
			auto uConR_acc = pfR.uConSYCL[q].get_access<cl::sycl::access::mode::read>(cgh);
			auto uPriL_acc = pfL.uPriSYCL[q].get_access<cl::sycl::access::mode::read>(cgh);
			auto uPriR_acc = pfR.uPriSYCL[q].get_access<cl::sycl::access::mode::read>(cgh);
			auto uConL_loc_acc = uConL.get_access<cl::sycl::access::mode::write>(cgh);
			auto uConR_loc_acc = uConR.get_access<cl::sycl::access::mode::write>(cgh);
			auto uPriL_loc_acc = uPriL.get_access<cl::sycl::access::mode::write>(cgh);
			auto uPriR_loc_acc = uPriR.get_access<cl::sycl::access::mode::write>(cgh);
			cgh.parallel_for<class SaveLocalState>(Range<1>(size + 3), Id<1>(1), [=](Item<1> item) {
				size_t i = item.get(0);
				uConR_loc_acc[i][q] = uConL_acc[i];
				uConL_loc_acc[i][q] = uConR_acc[i-1];
				uPriR_loc_acc[i][q] = uPriL_acc[i];
				uPriL_loc_acc[i][q] = uPriR_acc[i-1];
			});
		});
	}

	queue.submit([&](sycl::handler& cgh) {
		auto v_ch_m_acc = fl.v_ch_mSYCL.get_access<cl::sycl::access::mode::read>(cgh);
		auto v_ch_p_acc = fl.v_ch_pSYCL.get_access<cl::sycl::access::mode::read>(cgh);
		auto uConL_loc_acc = uConL.get_access<cl::sycl::access::mode::read>(cgh);
		auto uConR_loc_acc = uConR.get_access<cl::sycl::access::mode::read>(cgh);
		auto uConSL_loc_acc = uConSL.get_access<cl::sycl::access::mode::read_write>(cgh);
		auto uConSR_loc_acc = uConSR.get_access<cl::sycl::access::mode::read_write>(cgh);
		auto uPriL_loc_acc = uPriL.get_access<cl::sycl::access::mode::read>(cgh);
		auto uPriR_loc_acc = uPriR.get_access<cl::sycl::access::mode::read>(cgh);
		auto pthermL_acc = pfL.pthermSYCL.get_access<cl::sycl::access::mode::read>(cgh);
		auto pthermR_acc = pfR.pthermSYCL.get_access<cl::sycl::access::mode::read>(cgh);
		auto vCharS_acc = vCharS.get_access<cl::sycl::access::mode::write>(cgh);

		cgh.parallel_for<class FluxInsideRiemannCompute>(Range<1>(size + 3), Id<1>(1), [=](Item<1> item) {
			size_t i = item.get(0);
			if (v_ch_m_acc[i] <= 0.0) { // Riemann fan going right
				// branch handled in a previous kernel
			} else if (v_ch_p_acc[i] <= 0.0) { // Riemann fan going left
				// branch handled in a previous kernel
			} else {
				REAL rhoL = uConL_loc_acc[i][q_rho];
				REAL rhoR = uConR_loc_acc[i][q_rho];

				REAL pThermL(pthermR_acc[i-1]);
				REAL pThermR(pthermL_acc[i]);

				REAL vParL = uPriL_loc_acc[i][qvPar];
				REAL vParR = uPriR_loc_acc[i][qvPar];

				REAL sParL = uConL_loc_acc[i][qvPar];
				REAL sParR = uConR_loc_acc[i][qvPar];

				// Characteristic velocities like in Toro book
				REAL vCharL = -v_ch_m_acc[i];
				REAL vCharR = v_ch_p_acc[i];

				// Just the signal velocities relative to background flow
				REAL vSigL = vCharL - vParL;
				REAL vSigR = vCharR - vParR;

				// Velocity of contact discontinuity S*
				REAL idenom = 1. / (rhoL * vSigL - rhoR * vSigR + veps);
				// Eq. (10.70)
				REAL vCharS = (pThermR - pThermL + sParL * vSigL - sParR * vSigR) * idenom;

				// save for use in follow-up kernel
				vCharS_acc[i] = vCharS;

				// Starred quantities according to Eq. (10.73)
				uConSL_loc_acc[i][q_rho] = uConL_loc_acc[i][q_rho] * vSigL / (vCharL - vCharS);
				uConSR_loc_acc[i][q_rho] = uConR_loc_acc[i][q_rho] * vSigR / (vCharR - vCharS);

				// Velocities in middle region
				uConSL_loc_acc[i][qvPar] = uConSL_loc_acc[i][q_rho] * vCharS;
				uConSL_loc_acc[i][qvP1] = uConSL_loc_acc[i][q_rho] * uPriL_loc_acc[i][qvP1];
				uConSL_loc_acc[i][qvP2] = uConSL_loc_acc[i][q_rho] * uPriL_loc_acc[i][qvP2];

				uConSR_loc_acc[i][qvPar] = uConSR_loc_acc[i][q_rho] * vCharS;
				uConSR_loc_acc[i][qvP1] = uConSR_loc_acc[i][q_rho] * uPriR_loc_acc[i][qvP1];
				uConSR_loc_acc[i][qvP2] = uConSR_loc_acc[i][q_rho] * uPriR_loc_acc[i][qvP2];

				// Overall energy in middle region
				REAL egesL = uConL_loc_acc[i][q_Eges];
				uConSL_loc_acc[i][q_Eges] = uConSL_loc_acc[i][q_rho] * (egesL / rhoL + (vCharS - vParL) * (vCharS + pThermL / (rhoL * vSigL)));
				REAL egesR = uConR_loc_acc[i][q_Eges];
				uConSR_loc_acc[i][q_Eges] = uConSR_loc_acc[i][q_rho] * (egesR / rhoR + (vCharS - vParR) * (vCharS + pThermR / (rhoR * vSigR)));
			}
		});
	});

	for (int q = 0; q < n_omInt; ++q) {
		queue.submit([&](sycl::handler& cgh) {
			auto fluxL_acc = pfL.fluxSYCL[q].get_access<cl::sycl::access::mode::read>(cgh);
			auto fluxR_acc = pfR.fluxSYCL[q].get_access<cl::sycl::access::mode::read>(cgh);
			auto fluxF_acc = fl.fluxSYCL[q].get_access<cl::sycl::access::mode::write>(cgh);
			auto v_ch_m_acc = fl.v_ch_mSYCL.get_access<cl::sycl::access::mode::read>(cgh);
			auto v_ch_p_acc = fl.v_ch_pSYCL.get_access<cl::sycl::access::mode::read>(cgh);
			auto pthermL_acc = pfL.pthermSYCL.get_access<cl::sycl::access::mode::read>(cgh);
			auto pthermR_acc = pfR.pthermSYCL.get_access<cl::sycl::access::mode::read>(cgh);
			auto ptotalF_acc = fl.ptotalSYCL.get_access<cl::sycl::access::mode::write>(cgh);
			auto uConSL_loc_acc = uConSL.get_access<cl::sycl::access::mode::read>(cgh);
			auto uConSR_loc_acc = uConSR.get_access<cl::sycl::access::mode::read>(cgh);
			auto uConL_loc_acc = uConL.get_access<cl::sycl::access::mode::read>(cgh);
			auto uConR_loc_acc = uConR.get_access<cl::sycl::access::mode::read>(cgh);
			auto vCharS_acc = vCharS.get_access<cl::sycl::access::mode::read>(cgh);
			cgh.parallel_for<class FluxInsideRiemannSave>(Range<1>(size + 3), Id<1>(1), [=](Item<1> item) {
				size_t i = item.get(0);
				if (v_ch_m_acc[i] <= 0.0) { // Riemann fan going right
					// branch handled in a previous kernel
				} else if (v_ch_p_acc[i] <= 0.0) { // Riemann fan going left
					// branch handled in a previous kernel
				} else {
					// Fluxes from Eq (10.71)
					if (vCharS_acc[i] >= 0.0) {
						REAL vCharL = -v_ch_m_acc[i];

						//for (int q = 0; q < n_omInt; ++q) {
							fluxF_acc[i] = fluxR_acc[i - 1] + vCharL * (uConSL_loc_acc[i][q] - uConL_loc_acc[i][q]);
						//}
						ptotalF_acc[i] = pthermR_acc[i - 1];

					} else {
						REAL vCharR = v_ch_p_acc[i];

						//for (int q = 0; q < n_omInt; ++q) {
							fluxF_acc[i] = fluxL_acc[i] + vCharR * (uConSR_loc_acc[i][q] - uConR_loc_acc[i][q]);
						//}
						ptotalF_acc[i] = pthermL_acc[i];

					}
				}
			});
		});
	}

	//assert_sycl_eq(fl.ptotalORIG, fl.ptotalSYCL);
	//for (int q = 0; q < n_omInt; ++q) {
	//	assert_sycl_eq(fl.fluxORIG[q], fl.fluxSYCL[q]);
	//}

#endif // USE_SYCL

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

#endif
