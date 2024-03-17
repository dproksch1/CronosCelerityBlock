/**
 * @file
 * @brief Solver execution functions
 *
 * Provides the functionality for solver execution.
 */

#ifndef CRONOS_SOLVER_H
#define CRONOS_SOLVER_H = 1

#include "problem.H"
#include "Hdf5File_cbase.H"
#include "data.H"
#include "gridfunc.H"
#include "CException.H"
#include "queue.H"
#include "lazy_assertion.H"


/*!
 * @brief Generic Runge-Kutta solver class
 */
class RKSolver {
public:
	virtual double singlestep(Data &, GridFunc &, ProblemType &, int, Queue&) {return 0.;}
	virtual void phystest(Data &, GridFunc &, ProblemType &, int){}
	virtual void init_constants(Data &) {}
	virtual void addConstants_toH5(Data &, Hdf5Stream &h5out) {}
	virtual void getConstants_fromH5(Data &, Hdf5iStream &h5in) {}
	virtual void init(Queue &, Data &, GridFunc &, ProblemType &) {}
	virtual void restart(Queue &, Data &, GridFunc &, ProblemType &) {}
	virtual void Abort(CException) {}
	virtual ~RKSolver(){}
};

#endif