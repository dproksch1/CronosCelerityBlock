#ifndef CRONOS_RIEMANNSOLVERHD_H
#define CRONOS_RIEMANNSOLVERHD_H

#include "RiemannSolver.H"

#define HLLCSOLVER_HYDRO_VEPS 1.e-120

class RiemannSolverHD: public RiemannSolver {
public:
	RiemannSolverHD(){}
	RiemannSolverHD(const Data &gdata, int, int);
	virtual void get_vChar(const Data &, const ProblemType &, const phys_fields_0D &, const phys_fields_0D &,
				num_fields_0D &, int ix, int iy, int iz, int dir, REAL &cfl_lin) const;
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
	REAL veps;
	REAL gamma;
};

void get_vChar2(const Data &, const ProblemType &, const phys_fields_0D &, const phys_fields_0D &,
			num_fields_0D &, int dir, REAL &);

void get_NumFlux2(const Data &, const phys_fields_0D &, const phys_fields_0D &,
			num_fields_0D &, int dir, int iFluid=0);

#endif
