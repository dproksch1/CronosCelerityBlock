#ifndef CRONOS_SOLVER_H
#define CRONOS_SOLVER_H = 1

#include "problem.H"
#include "Hdf5File_cbase.H"
#include "data.H"
#include "gridfunc.H"
#include "CException.H"
#include "queue.H"

class rksolver {
public:
	virtual REAL singlestep(Data &, gridFunc &, ProblemType &, int, Queue&) {return 0.;}
	virtual void phystest(Data &, gridFunc &, ProblemType &, int){}
	virtual void init_constants(Data &) {}
	virtual void addConstants_toH5(Data &, Hdf5Stream &h5out) {}
	virtual void getConstants_fromH5(Data &, Hdf5iStream &h5in) {}
	virtual void init(Data &, gridFunc &, ProblemType &) {}
	virtual void restart(Data &, gridFunc &, ProblemType &) {}
	virtual void Abort(CException) {}
	virtual ~rksolver(){}
	REAL dt_idl, dt_mov;
};

#endif