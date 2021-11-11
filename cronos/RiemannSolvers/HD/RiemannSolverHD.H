#ifndef CRONOS_RIEMANNSOLVERHD_H
#define CRONOS_RIEMANNSOLVERHD_H

#include "RiemannSolver.H"

#define HLLCSOLVER_HYDRO_VEPS 1.e-120

class RiemannSolverHD: public RiemannSolver {
public:
	RiemannSolverHD(){}
	RiemannSolverHD(const Data &gdata, int, int);
	void get_vChar(Queue /*queue*/, Data &, ProblemType &, cronos::vector<double> &,
			phys_fields_1D &, phys_fields_1D &, NumMatrix<REAL,1> &, NumMatrix<REAL,1> &,
			Buffer<REAL, 1> /*v_ch_mSYCL*/, Buffer<REAL, 1> /*v_ch_pSYCL*/,
//	               phys_fields_1D &, phys_fields_1D &, fields_1D &,
	               const int &, REAL &) const;
	virtual void get_vChar(const Data &, const ProblemType &, const phys_fields_0D &, const phys_fields_0D &,
				num_fields_0D &, int ix, int iy, int iz, int dir, REAL &cfl_lin) const;
//	void reset_Indices(ProblemType &Problem) {
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


// Generic HLL solver
//class HLLSolver_gen: public RiemannSolverHD {
//public:
//	HLLSolver_gen(const Data &, int, int);
//	void get_NumFlux(Queue /*queue*/, Data &, phys_fields_1D &,
//	                 phys_fields_1D &, fields_1D &, const int &, int iFluid=0) const;
//	virtual void get_NumFlux(Data &, phys_fields_0D &,
//			phys_fields_0D &, num_fields_0D &, int dir, int iFluid=0) const;
//	~HLLSolver_gen(){};
//private:
//	REAL veps;
//};


// HLLC Solver for Hydrodynamics
class HLLCSolver_Hydro: public RiemannSolverHD {
public:
	HLLCSolver_Hydro(const Data &, int, int);
	void get_NumFlux(Queue /*queue*/, Data &, phys_fields_1D &,
	                 phys_fields_1D &, fields_1D &, const int &, int iFluid=0) const;
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
