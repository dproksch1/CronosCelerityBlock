#ifndef CRONOS_PROJ_H_INCLUDED
#define CRONOS_PROJ_H_INCLUDED

#ifdef _CRAY
typedef long Long;
#else
typedef long long Long;
#endif

#include <time.h>
#include <stdlib.h>
#include <fstream>
#include <sys/types.h>
#include <sys/stat.h>
#include "matrix.H"
#include "data.H"
#include "gridfunc.H"
#include "Hdf5File_cbase.H"
#if (CRONOS_MOVIE == CRONOS_ON)
#include "movie.H"
#endif
#include "constants.H"


#include "specific.H"
#include "CException.H"
#include "problem.H"
#include "solver.H"

class Environment
{
public:
	Environment(Data &);
	~Environment();
	void setType(Data &);
	int integrate(Data &); //used
	void singlestep(Data &);
	void compute_dt(Data &, double scale=1.); //used
	void LoadData(Data &); //not used regularily, still leave in
#if(FLUID_TYPE != CRONOS_MULTIFLUID)
	void LoadData_flt(Data &, int); //not used atm
#endif
	void InitOutput(Data &); //used
	void Abort(Data &, CException); //not used regularily, still leave in
	int Finalize(Data &, string); //used
	unique_ptr<rksolver> RKSolver;
//	rksolver *RKSolver;
	/**vector of Euler solvers*/
	//vector< unique_ptr<eulersolver> > EulerSolvers;
	//	eulersolver *EulerSolver;
	int outputflag;
	REAL restart_time;
	std::unique_ptr<ProblemType> Problem;
	std::unique_ptr<gridFunc> gfunc;
private:
	void setup(Data &); //used
	void init_solvers(Data &);
	//   void set_MovieOut(Data &, int);
#if (CRONOS_MOVIE == CRONOS_ON)
	void WriteMovies(Data &);
	Movie mov;
#endif
	void pdestep(Data &);
	bool CheckEnd(Data &); //used
	bool CheckEnd_User(Data &); //used
	void CheckOut(Data &); //used
	//	void Finalize(Data &, string);
	void WriteDivB(Data &gdata); //not used atm
	void Output(Data &, bool, bool); //used
	string MakeFilename(const int [DIM], int, bool, bool); //used
	void AddParallelFileNameStuff(string &filename, const int [DIM]); //not used
	void getCoords(const Data &, string, int [DIM], int [DIM]); //part of loadData
	timeval tick,tock;
	int numdbl_out, numdbl_done, numdbl_pass;
	int numflt_out, numflt_done, numflt_pass;
	int numascii_out, numascii_done, numascii_pass;
	int numinfo_out, numinfo_done, numinfo_pass;
	int failnum, dt_debug;
	int restart_step;
	int nummov_done;
	int N_MOV;
	REAL dt_dbl, dt_flt, dt_ascii, dt_info, dt_mov;
	REAL t_last_dbl, t_last_flt;
	REAL cfl_set, cfl_min, cfl_max;
	string OldName;
};



#endif
