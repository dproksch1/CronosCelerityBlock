/**
 * @file
 * @brief A collection of GPU-compatible functions for the main kernel
 *
 * All functions used for calculations in gpu::compute_step() aside from reconstruction and proprietary Celerity/SYCL math functions
 */

#ifndef CRONOS_SOLVER_FUNC_H
#define CRONOS_SOLVER_FUNC_H

#include "data.H"
#include "problem.H"
#include "vector.H"
#include "fields_local.H"
#include "utils.H"

#include <stdlib.h>

#define HLLCSOLVER_HYDRO_VEPS 1.e-120

namespace gpu {

	//! @brief Compute the physical flux based on primitive and conservation values
	template<typename PhysFluxAccessorType, typename UPriAccessorType, typename UConAccessorType, typename SinglePthermAccessorType, typename FluidConstType>
	void get_PhysFlux(PhysFluxAccessorType& physFlux, UPriAccessorType& uPri, UConAccessorType& uCon,
					SinglePthermAccessorType& ptherm, int face, FluidConstType fluidConst) {
		
		int dir = face/2;

		physFlux[fluidConst[fluidConst_q_sx]] = 0.;
		physFlux[fluidConst[fluidConst_q_sy]] = 0.;
		physFlux[fluidConst[fluidConst_q_sz]] = 0.;

		// flux for the density
		physFlux[fluidConst[fluidConst_q_rho]] = uPri[fluidConst[fluidConst_q_sx] + dir] * uPri[fluidConst[fluidConst_q_rho]];
		
		// flux for momentum vector 
#if EXTRACT_PRESSURE == TRUE
		for(int q = fluidConst[fluidConst_q_sx]; q <= fluidConst[fluidConst_q_sz]; ++q) {
			physFlux[q] = (uCon[q] * uPri[dir + fluidConst[fluidConst_q_sx]]);
		}
#else
		
		for(int q = fluidConst[fluidConst_q_sx]; q <= fluidConst[fluidConst_q_sz]; ++q) {
			if(q == dir+1) {
				physFlux[q] += ((uPri[q] * uCon[q]) + ptherm);
			} else {
				physFlux[q] = (uCon[q] * uPri[dir + fluidConst[fluidConst_q_sx]]);
			}
		}
		
#endif

		// flux for total energy
		if(ENERGETICS == FULL){
			physFlux[fluidConst[fluidConst_q_Eges]] = (uCon[fluidConst[fluidConst_q_Eges]] + ptherm)
						* uPri[fluidConst[fluidConst_q_sx] + dir];
		}

// #if(CRSWITCH_DUAL_ENERGY == CRONOS_ON)
// 		// Compute physical flux for entropy evolution -- flux is
// 		// Entropy*v_par
// 		physFlux[fluidConst[fluidConst_q_Eadd]] = uCon[fluidConst[fluidConst_q_Eadd]] * uPri[fluidConst[fluidConst_q_sx] + dir];
// #endif

	}

	//! @brief Compute cfl numbers and characteristic velocities
    template<typename UPriAccessorType, typename UConAccessorType, typename PThermAccessorType, typename NumValsAccessorType, typename FluidConstType>
	double get_vChar(UPriAccessorType &pfL_uPri, UConAccessorType &pfL_uCon, PThermAccessorType &pfL_ptherm,
			UPriAccessorType &pfR_uPri, UConAccessorType &pfR_uCon, PThermAccessorType &pfR_ptherm, NumValsAccessorType &f_num,
			int dir, double problem_gamma, const double idx[], FluidConstType fluidConst) {
		//! Compute characteristic velocities

		int shift_vec[3] = {0,0,0};
		shift_vec[dir] = -1;

		//int iPos[3] = {ix, iy, iz};
		int q_rho = fluidConst[fluidConst_q_rho];
		int q_sx = fluidConst[fluidConst_q_sx];

		double rhoinv_p = 1./pfL_uCon[q_rho];
		double rhoinv_m = 1./pfR_uCon[q_rho];

		// Flow velocity
		double u_p = pfL_uPri[dir+q_sx];
		double u_m = pfR_uPri[dir+q_sx];

		double pres_p = pfL_ptherm;
		double pres_m = pfR_ptherm;

		// Sound speed
		double cs_p   = cl::sycl::sqrt(problem_gamma*pres_p*rhoinv_p);
		double cs_m   = cl::sycl::sqrt(problem_gamma*pres_m*rhoinv_m);

		double v_ch_p = cl::sycl::fmax(cl::sycl::fmax(cs_p+u_p,cs_m+u_m),0.);
		double v_ch_m = cl::sycl::fmax(cl::sycl::fmax(cs_p-u_p,cs_m-u_m),0.);

		f_num[NumV_ch_p] = v_ch_p;
		f_num[NumV_ch_m] = v_ch_m;

		double vmax = cl::sycl::fmax(v_ch_p, v_ch_m);

		// Local computation of cfl number
	#if  (NON_LINEAR_GRID == CRONOS_OFF)
		double cfl_loc = vmax * idx[dir];
	#else
		//double cfl_loc = vmax*gdata.getCen_idx(dir, iPos[dir]);
		double cfl_loc = vmax*idx[dir];
	#endif

	#ifdef GEOM
	#if GEOM != CARTESIAN
		//Non-Cartesian Case not implemented
	#endif
	#endif
		
		return cfl_loc;
	}

	//! @brief Compute numerical flux from phys fluxes using hllc solver
	template<typename UPriAccessorType, typename UConAccessorType, typename PhysFluxAccessorType, typename SinglePthermAccessorType, typename SinglePtotalAccessorType, typename CFlagAccessorType, typename NumFluxAccessorType, typename NumValsAccessorType, typename FluidConstType>
	double get_NumFlux(const UPriAccessorType &pfM_uPri, const UConAccessorType &pfM_uCon, const PhysFluxAccessorType& pfM_flux,
				const SinglePthermAccessorType pfM_ptherm, const SinglePtotalAccessorType& pfM_ptotal, const CFlagAccessorType pfM_cflag,
				const UPriAccessorType &pfP_uPri, const UConAccessorType &pfP_uCon, const PhysFluxAccessorType& pfP_flux,
				const SinglePthermAccessorType pfP_ptherm, const SinglePtotalAccessorType& pfP_ptotal, const CFlagAccessorType pfP_cflag,
				NumFluxAccessorType& numFlux, const NumValsAccessorType numVals, bool use_carbuncle, int dir, const double gamma, FluidConstType fluidConst) {
		//! Compute numerical flux from phys fluxes using hllc solver
		/*
		* @PARAM pfM variables and fluxes on left-hand side of interface
		* @PARAM pfP variables and fluxes on right-hand side of interface
		*/

		double num_ptotal;
		int qvPar, qvP1, qvP2;
		int q_sx = fluidConst[fluidConst_q_sx];
		int q_sy = fluidConst[fluidConst_q_sy];
		int q_sz = fluidConst[fluidConst_q_sz];

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

		/*-------------------------------------------------------
		  Compute fluxes for regions outside the Riemann fan:
		  -----------------------------------------------------*/

		if(numVals[NumV_ch_m] <= 0.) { // Riemann fan going right

			for(int q=0; q<N_OMINT; ++q) {
				numFlux[q] = pfM_flux[q];
			}
			
			num_ptotal = pfM_ptotal;

		} else if (numVals[NumV_ch_p] <= 0.) { // Riemann fan going left

			for(int q=0; q<N_OMINT; ++q) {
				numFlux[q] = pfP_flux[q];
			}

			num_ptotal = pfP_ptotal;

		} else {

			/*-------------------------------------------------------
				The rest is taking place inside the Riemann fan
				-------------------------------------------------------*/

			// If carbuncle test is used -> check whether cells i and i-1 are
			// flagged as lying within the vicinity of a strong shock

			//TODO implement
			/*if(gdata.use_carbuncleFlag) {
				// if both cells are flagged -> use hll for all variables instead
				if(pfM.carbuncle_flag*pfP.carbuncle_flag > 0) {
					for(int q=0; q<pfM.get_num(); ++q) {
						double fac = 1./(f_num.v_ch_p + f_num.v_ch_m+ HLLCSOLVER_HYDRO_VEPS);

						f_num.flux_num[q] = (f_num.v_ch_m*pfP.flux_phys[q] +
								f_num.v_ch_p*pfM.flux_phys[q] -
								f_num.v_ch_m*f_num.v_ch_p*(pfP.uCon[q] - pfM.uCon[q]))*fac;
					}
					return;
				}
			}*/



			if(use_carbuncle) {
				// if both cells are flagged -> use hll for all variables instead
				if(pfM_cflag*pfP_cflag > 0) {
					for(int q=0; q < N_OMINT; ++q) {
						double fac = 1./(numVals[NumV_ch_p] + numVals[NumV_ch_m] + HLLCSOLVER_HYDRO_VEPS);

						numFlux[q] = (numVals[NumV_ch_m]*pfP_flux[q] +
								numVals[NumV_ch_p]*pfM_flux[q] -
								numVals[NumV_ch_m]*numVals[NumV_ch_p]*(pfP_uCon[q] - pfM_uCon[q]))*fac;
					}
					return 0;
				}
			}

			double uConL[N_OMINT], uConR[N_OMINT];
			//std::vector<double> uConL(n_omInt);
			//std::vector<double> uConR(n_omInt);
			double uPriL[N_OMINT], uPriR[N_OMINT];
			//std::vector<double> uPriL(n_omInt);
			//std::vector<double> uPriR(n_omInt);
			double uConSL[N_OMINT];
			//std::vector<double> uConSL(n_omInt);
			double uConSR[N_OMINT];
			//std::vector<double> uConSR(n_omInt);

			// Saving array values:

			for(int q=0; q<N_OMINT; ++q) {
				/*
					In this context "L" and "M" mean (L)eft from the cell-face
					or on the (M)inus side
					The same holds for "R" and "P"
				*/
				uConR[q] = pfP_uCon[q];
				uConL[q] = pfM_uCon[q];
				uPriR[q] = pfP_uPri[q];
				uPriL[q] = pfM_uPri[q];
			}

			int q_rho = fluidConst[fluidConst_q_rho];

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


			double rhoL = uConL[q_rho];
			double rhoR = uConR[q_rho];

			double pThermL(pfM_ptherm);
			double pThermR(pfP_ptherm);

			double vParL = uPriL[qvPar];
			double vParR = uPriR[qvPar];

			double sParL = uConL[qvPar];
			double sParR = uConR[qvPar];

			// Characteristic velocities like in Toro book
			double vCharL = -numVals[NumV_ch_m];
			double vCharR =  numVals[NumV_ch_p];

			// Just the signal velocities relative to background flow
			double vSigL = vCharL - vParL;
			double vSigR = vCharR - vParR;

			// Velocity of contact discontinuity S*
			double idenom = 1./(rhoL*vSigL - rhoR*vSigR + HLLCSOLVER_HYDRO_VEPS);
			// Eq. (10.70)
			double vCharS = (pThermR - pThermL +
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
			int q_Eges = fluidConst[fluidConst_q_Eges];
			double egesL = uConL[q_Eges];
			uConSL[q_Eges] = uConSL[q_rho]*(egesL/rhoL + (vCharS - vParL)*
					(vCharS + pThermL/(rhoL*vSigL)));
			double egesR = uConR[q_Eges];
			uConSR[q_Eges] = uConSR[q_rho]*(egesR/rhoR + (vCharS - vParR)*
					(vCharS + pThermR/(rhoR*vSigR)));

#if(CRSWITCH_DUAL_ENERGY == CRONOS_ON)
		// Compute the entropy in the starred region for dual
			// energy description

			// Pressure in starred region (according to Eq. (10.36) in
			// Toro book) - should be identical for L and R values
			double pThermS = pThermL + rhoL*(vSigL - vParL)*(vCharS - vParL);

			// With this compute entropy:
			uConSL[q_Eadd] = pThermS/cl::sycl::pow(uConSL[q_rho], gamma-1.);
			uConSR[q_Eadd] = pThermS/cl::sycl::pow(uConSR[q_rho], gamma-1.);
#endif


			// Fluxes from Eq (10.71)
			if(vCharS >= 0.) {

				for(int q=0; q<N_OMINT; ++q) {
					numFlux[q] = pfM_flux[q] + vCharL*(uConSL[q] -
							uConL[q]);
				}

				// fl.ptotal(i) = ptotalL(i);
				
				num_ptotal = pfM_ptherm;

			} else {

				for(int q=0; q<N_OMINT; ++q) {
					numFlux[q] = pfP_flux[q] + vCharR*(uConSR[q] -
							uConR[q]);
				}
				// fl.ptotal(i) = ptotalR(i);
				
				num_ptotal = pfP_ptherm;
			}
		}

		return num_ptotal;
	}

}

void compute_carbuncleFlag(Queue &, Data &, double alpha_c);

#endif