#ifndef CRONOS_RIEMANNSOLVERHD_H
#define CRONOS_RIEMANNSOLVERHD_H

#include "RiemannSolver.H"

#include "utils.H"

#define HLLCSOLVER_HYDRO_VEPS 1.e-120

namespace gpu {

	template<typename AccessorType, typename UPriAccessorType, typename UConAccessorType, typename PThermAccessorType, typename NumValsAccessorType, typename IdxConstType, typename FluidConstType>
	double get_vChar(AccessorType& omAcc, UPriAccessorType &pfL_uPri, UConAccessorType &pfL_uCon, PThermAccessorType &pfL_ptherm,
			UPriAccessorType &pfR_uPri, UConAccessorType &pfR_uCon, PThermAccessorType &pfR_ptherm, NumValsAccessorType &f_num,
			int dir, double problem_gamma, IdxConstType idx, FluidConstType fluidConst) {
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
		double cfl_loc = vmax*idx[dir];
	#else
		//dproksch TODO implement
		//double cfl_loc = vmax*gdata.getCen_idx(dir, iPos[dir]);
		double cfl_loc = vmax*idx[dir];
	#endif
//dproksch TODO implement
	#ifdef GEOM
	#if GEOM != CARTESIAN
		/*if(dir==0) {
			cfl_loc /= gdata.h0(ix, iy, iz, shift_vec[0],shift_vec[1],shift_vec[2]);
		} else if (dir==1) {
			cfl_loc /= gdata.h1(ix, iy, iz, shift_vec[0],shift_vec[1],shift_vec[2]);
		} else {
			cfl_loc /= gdata.h2(ix, iy, iz, shift_vec[0],shift_vec[1],shift_vec[2]);
		}*/
	#endif
	#endif
		
		return cfl_loc;
	}


}

class RiemannSolverHD: public RiemannSolver {
public:
	RiemannSolverHD(){}
	RiemannSolverHD(const Data &gdata, int, int);
	virtual void get_vChar(const Data &, const ProblemType &, const phys_fields_0D &, const phys_fields_0D &,
				num_fields_0D &, int ix, int iy, int iz, int dir, double &cfl_lin) const;
	void reset_Indices(const CronosFluid &fluid) {
		//! Reset field indices (only necessary for multifluid simulations)
		this->q_rho = fluid.get_q_rho();
		this->q_sx = fluid.get_q_sx();
		this->q_sy = fluid.get_q_sy();
		this->q_sz = fluid.get_q_sz();
		this->q_Eges = fluid.get_q_Eges();
		this->q_Eadd = fluid.get_q_Eadd();
	};
protected:
	int q_rho, q_sx, q_sy, q_sz, q_Eges, q_Eadd;
};

// HLLC Solver for Hydrodynamics
class HLLCSolver_Hydro: public RiemannSolverHD {
public:
	HLLCSolver_Hydro(const Data &, int, int);
	virtual void get_NumFlux(const Data &, const phys_fields_0D &,
		const phys_fields_0D &, num_fields_0D &, int dir, int iFluid=0) const;
	~HLLCSolver_Hydro(){};
private:
	int qvPar, qvP1, qvP2;
	double veps;
	double gamma;
};

void get_vChar2(const Data &, const ProblemType &, const phys_fields_0D &, const phys_fields_0D &,
			num_fields_0D &, int dir, double &);

double get_vChar3(const Data &, const ProblemType &, const phys_fields_0D &, const phys_fields_0D &,
			num_fields_0D &, int dir);

void get_NumFlux2(const Data &, const phys_fields_0D &, const phys_fields_0D &,
			num_fields_0D &, int dir, int iFluid=0);

#endif
