#ifndef CRONOS_TRANSFORMATIONS_BLOCK_H
#define CRONOS_TRANSFORMATIONS_BLOCK_H

#include "CException.H"
#include "data.H"
#include "gridfunc.H"
#include "problem.H"
#include "eos.H"
#include "vector.H"
#include "fields_local.H"
#include "utils.H"

namespace gpu {

	inline double pressure(double dens, const double problem_cs2, const double gamma, const double denominator, const double half_beta) {
#ifdef SOUND
		if(gamma != 1.){
			return problem_cs2*cl::sycl::pow(dens, gamma)*denominator;
		} else {
			return problem_cs2*dens;
		}
#else
		if(gamma != 1.){
			return problem_cs2*cl::sycl::pow(dens, gamma)*denominator*half_beta;
		} else {
			return problem_cs2*dens*half_beta;
		}
#endif
	}

	inline double TransEth2E_gpu(double rhoinv, double psq, double Bsq, double ETherm, const int fluidType) {
		
		#if (FLUID_TYPE == CRONOS_MHD)
			double Energy = ETherm + 0.5*psq*rhoinv + 0.5*Bsq;
		#elif (FLUID_TYPE == CRONOS_HYDRO)
			double Energy = ETherm + 0.5*psq*rhoinv + 0.5*Bsq;
		#else
			double Energy;
			if(fluidType==CRONOS_MHD) {
				Energy = ETherm + 0.5*psq*rhoinv + 0.5*Bsq;
			} else {
				Energy = ETherm + 0.5*psq*rhoinv;
			}
		#endif
		
		return Energy;
	}

	inline double TransT2Eth_gpu(double problem_gamma, double rhoinv, double ETherm) {

		return 1./(rhoinv*(problem_gamma-1.))*ETherm;
	}

	template<typename UPriAccessorType, typename UConAccessorType, typename SinglePthermAccessorType, typename SinglePtotalAccessorType, typename FluidConstType>
	void get_Cons(UPriAccessorType& uPri, UConAccessorType& uCon, SinglePthermAccessorType& ptherm,
				 SinglePtotalAccessorType& ptotal, int ix, int iy, int iz, int face, int thermal, double problem_gamma,
				 const double problem_cs2, const double denominator, const double half_beta, const int fluidType, FluidConstType fluidConst) {
		
		int dir = face / 2;
		
		// Compute thermal pressure:
		if(ENERGETICS == FULL) {
			if(thermal) { // Compute thermal pressure from Etherm
				ptherm = (problem_gamma - 1.)*uPri[fluidConst[fluidConst_q_Eges]];
			} else { // Compute thermal pressure from Temperature
				ptherm = uPri[fluidConst[fluidConst_q_rho]]*uPri[fluidConst[fluidConst_q_Eges]];
			}
		} else {
			ptherm = pressure(uPri[fluidConst[fluidConst_q_rho]], problem_cs2, problem_gamma, denominator, half_beta);
		}

		// Set / compute conservative variables:

		uCon[fluidConst[fluidConst_q_rho]] = uPri[fluidConst[fluidConst_q_rho]]; // Density

		uCon[fluidConst[fluidConst_q_sx]] = uPri[fluidConst[fluidConst_q_sx]]*uPri[fluidConst[fluidConst_q_rho]]; // v_x -> s_x
		uCon[fluidConst[fluidConst_q_sy]] = uPri[fluidConst[fluidConst_q_sy]]*uPri[fluidConst[fluidConst_q_rho]]; // v_y -> s_y
		uCon[fluidConst[fluidConst_q_sz]] = uPri[fluidConst[fluidConst_q_sz]]*uPri[fluidConst[fluidConst_q_rho]]; // v_z -> s_z


		double Bsq(0.);
	// Thermal energy / temperature / overall energy
		if(ENERGETICS == FULL) {
			double rhoinv = 1./uPri[fluidConst[fluidConst_q_rho]];

			double psq = (sqr(uPri[fluidConst[fluidConst_q_sx]]) + sqr(uPri[fluidConst[fluidConst_q_sy]]) +
					sqr(uPri[fluidConst[fluidConst_q_sz]]))*sqr(uPri[fluidConst[fluidConst_q_rho]]);

			if(thermal) {
				uCon[fluidConst[fluidConst_q_Eges]] = TransEth2E_gpu(rhoinv,psq,Bsq,
						uPri[fluidConst[fluidConst_q_Eges]], fluidType);
			} else {
				double Energy = TransT2Eth_gpu(problem_gamma, rhoinv, uPri[fluidConst[fluidConst_q_Eges]]);
				uCon[fluidConst[fluidConst_q_Eges]] = TransEth2E_gpu(rhoinv,psq,Bsq, Energy, fluidType);
			}

#if(CRSWITCH_DUAL_ENERGY == CRONOS_ON)
			// Computation of entropy for dual energy corrections
			double Etherm = ptherm/(problem_gamma - 1.);
			uPri[fluidConst[fluidConst_q_Eadd]] = Etherm*cl::sycl::pow(uCon[fluidConst[fluidConst_q_rho]],
					1.-problem_gamma);
			uCon[fluidConst[fluidConst_q_Eadd]] = uPri[fluidConst[fluidConst_q_Eadd]];
#endif

		}

		// Total pressure:
		ptotal = ptherm;

	}

}

class Transformations_Block {
public:
	Transformations_Block(const CronosFluid&, ProblemType &, bool, int iFluid=0);

    void set_thermal(bool thermal);
	bool get_thermal();

	void TransPrim2Cons(Queue &, Data &gdata, gridFunc &gfunc, ProblemType &Problem);
	void TransCons2Prim(Queue &, Data &gdata, gridFunc &gfunc, ProblemType &Problem);
	void TransVel2Momen(Queue &, Data &, gridFunc &, ProblemType &);
	void TransMomen2Vel(Queue &, Data &, gridFunc &, ProblemType &);
	void TransE2Eth(Queue &, Data &, gridFunc &, ProblemType &, int, bool);
	void TransEth2E(Queue &, const Data &gdata, gridFunc &gfunc, ProblemType &Problem) const;
	void TransT2Eth(Queue &, const Data &, gridFunc &, ProblemType &) const;
	void TransT2E(Queue &, const Data &, gridFunc &, ProblemType &) const;

private:
    double TNorm;
	double kBoverMeanMolWeight_num;
#if (USE_COROTATION == CRONOS_ON)
	double omegaZ;
#endif
	int TempErr, thermal;
	int fluidType, iFluid;
	int q_rho, q_sx, q_sy, q_sz, q_Bx, q_By, q_Bz, q_Eges, q_Eadd;
	int q_rho_loc, q_sx_loc, q_sy_loc, q_sz_loc;
	int q_Bx_loc, q_By_loc, q_Bz_loc, q_Eges_loc, q_Eadd_loc;
	bool TPhys;
	bool magFluid;

};

#endif