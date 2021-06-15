#ifndef CRONOS_PROJ_H_INCLUDED
#define CRONOS_PROJ_H_INCLUDED

#ifdef _CRAY
typedef long Long;
#else
typedef long long Long;
#endif

#ifdef parallel
#include "mpi.h"
#ifdef FFT
#include <rfftw_mpi.h>
#endif
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
#include "problem.H"

#include "queue.H"

class Environment
{
public:
	Environment(Data &);
	~Environment();
	void setType(Data &);
	int integrate(Data &, Queue& queue);
	void singlestep(Data &);
	void compute_dt(Data &, double scale=1.);
	void LoadData(Data &);
#if(FLUID_TYPE != CRONOS_MULTIFLUID)
	void LoadData_flt(Data &, int);
#endif
	void InitOutput(Data &);
	void Abort(Data &, CException);
	int Finalize(Data &, string);
	unique_ptr<rksolver> RKSolver;
//	rksolver *RKSolver;
	/**vector of Euler solvers*/
	vector< unique_ptr<eulersolver> > EulerSolvers;
	//	eulersolver *EulerSolver;
	int outputflag;
	REAL restart_time;
	std::unique_ptr<ProblemType> Problem;
	std::unique_ptr<gridFunc> gfunc;
private:
	void setup(Data &);
	void init_solvers(Data &);
	//   void set_MovieOut(Data &, int);
#if (CRONOS_MOVIE == CRONOS_ON)
	void WriteMovies(Data &);
	Movie mov;
#endif
	void pdestep(Data &, Queue& queue);
	bool CheckEnd(Data &);
	bool CheckEnd_User(Data &);
	void CheckOut(Data &);
	//	void Finalize(Data &, string);
	void WriteDivB(Data &gdata);
	void Output(Data &, bool, bool);
	string MakeFilename(const int [DIM], int, bool, bool);
	void AddParallelFileNameStuff(string &filename, const int [DIM]);
	void getCoords(const Data &, string, int [DIM], int [DIM]);
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
