#ifndef CRONOS_RIEMANNSOLVERMHD_H
#define CRONOS_RIEMANNSOLVERMHD_H

#include "RiemannSolver.H"

class RiemannSolverMHD: public RiemannSolver {
public:
	RiemannSolverMHD(){}
	RiemannSolverMHD(const Data &gdata, int, int);
	void get_vChar(Queue /*queue*/, Data &, ProblemType &, cronos::vector<double> &,
			phys_fields_1D &, phys_fields_1D &, NumMatrix<REAL,1> &, NumMatrix<REAL,1> &,
			Buffer<REAL, 1> /*v_ch_mSYCL*/, Buffer<REAL, 1> /*v_ch_pSYCL*/,
			//	               phys_fields_1D &, phys_fields_1D &, fields_1D &,
	               const int &, REAL &) const;
	virtual void get_vChar(const Data &, const ProblemType &, const phys_fields_0D &, const phys_fields_0D &,
			num_fields_0D &, int ix, int iy, int iz, int dir, REAL &cfl_lin) const;
	void get_vChar2D(Data &, ProblemType &, cronos::vector<double> &,
	                 phys_fields_2D &, phys_fields_2D &,
	                 phys_fields_2D &, phys_fields_2D &,
	                 fields_2D &, const int &, REAL &);
	virtual void get_NumEmf2D(Data &, phys_fields_2D &, phys_fields_2D &,
	                          phys_fields_2D &, phys_fields_2D &, fields_2D &,
	                          const int &){};
	virtual void GetEx2D(Data &, fields_2D &, NumMatrix<REAL,3> &,
	                     const int &, const int &){};
	void reset_Indices(const CronosFluid &fluid) {
		//! Reset field indices (only necessary for multifluid simulations)
		this->q_rho = fluid.get_q_rho();
		this->q_sx = fluid.get_q_sx();
		this->q_sy = fluid.get_q_sy();
		this->q_sz = fluid.get_q_sz();
		this->q_Bx = fluid.get_q_Bx();
		this->q_By = fluid.get_q_By();
		this->q_Bz = fluid.get_q_Bz();
		this->q_Eges = fluid.get_q_Eges();
		this->q_Eadd = fluid.get_q_Eadd();
	};
protected:
	int q_rho, q_sx, q_sy, q_sz, q_Bx, q_By, q_Bz, q_Eges, q_Eadd;
};


class HLLSolver: public RiemannSolverMHD {
public:
	HLLSolver(const Data &, int, int);
	HLLSolver(){}
	void get_NumFlux(Queue /*queue*/, Data &, phys_fields_1D &,
	                 phys_fields_1D &, fields_1D &, const int &, int iFluid=0) const;
	virtual void get_NumFlux(const Data &, const phys_fields_0D &, const phys_fields_0D &,
			num_fields_0D &, int dir, int iFluid=0) const;
	void get_NumEmf2D(Data &, phys_fields_2D &, phys_fields_2D &,
	                  phys_fields_2D &, phys_fields_2D &, fields_2D &,
	                  const int &);
	void GetEx2D(Data &, fields_2D &, NumMatrix<REAL,3> &,
	             const int &, const int &){};
	~HLLSolver(){};
private:
	REAL veps;
};


//class HLLDSolver: public RiemannSolverMHD {
//public:
//	HLLDSolver(const Data &, int, int);
//	void get_NumFlux(Queue /*queue*/, Data &, phys_fields_1D &,
//	                 phys_fields_1D &, fields_1D &, const int &, int iFluid=0) const;
//	virtual void get_NumFlux(Data &, phys_fields_0D &, phys_fields_0D &,
//			num_fields_0D &, int dir, int iFluid=0) const;
//	~HLLDSolver(){};
//private:
//	int qvPar, qvP1, qvP2, qBPar, qBP1, qBP2;
////	int RHO, qEges;
//	int qBMin, qBMax;
//	REAL veps, loceps;
//	REAL gamma;
//};

#endif
