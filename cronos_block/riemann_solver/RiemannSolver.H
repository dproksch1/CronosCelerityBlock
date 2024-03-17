/**
 * @file
 * @brief The solver configuration functionality
 *
 * This file provides solver configuration functionality via object-oriented solver objects.
 */

#ifndef CRONOS_RIEMANNSOLVER_H
#define CRONOS_RIEMANNSOLVER_H = 1

#include "data.H"
#include "problem.H"
#include "vector.H"
#include "fields_local.H"
#include "utils.H"
#include "solver_func.H"
#include <stdlib.h>


/*!
 * @brief Base class for Riemann solvers supplying base variables and interface
 */
class RiemannSolver {
public:
	RiemannSolver(){
		if(value_exists("use_carbuncleFlag")) {
			alpha_carbuncle = value((char*)"alpha_carbuncle");
		}
		verbosity = 0;
	};
	RiemannSolver(const Data &gdata, int /*dir*/, int Fluid_Type) {
		this->Fluid_Type = Fluid_Type;
		verbosity = 0;

		if(value_exists("use_carbuncleFlag")) {
			alpha_carbuncle = value((char*)"alpha_carbuncle");
		}
	};
	virtual ~RiemannSolver(){};
	virtual void reset_Indices(const CronosFluid &) = 0;
	int get_Fluid_Type() {
		return Fluid_Type;
	}
	double getAlphaCarbuncle();
	void set_verbosity(int _verb);
	
protected:
	double norm, alpha_carbuncle;
	int Fluid_Type;
	int verbosity;
	bool use_carbuncleFlag;
};

class RiemannSolverHD: public RiemannSolver {
public:
	RiemannSolverHD(){}
	RiemannSolverHD(const Data &gdata, int, int);
	
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


/*!
 * @brief HLLC Solver for Hydrodynamics
 */
class HLLCSolver_Hydro: public RiemannSolverHD {
public:
	HLLCSolver_Hydro(const Data &, int, int);
	~HLLCSolver_Hydro(){};
private:
	int qvPar, qvP1, qvP2;
	double veps;
	double gamma;
};


#endif
