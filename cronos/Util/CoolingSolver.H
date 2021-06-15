#ifndef COOLING_SOLVER_H
#define COOLING_SOLVER_H


//#include "constants.H"
#include "CoolingFunction.H"
#include "solver.H"
#include "normalisation.H"
//#include "data.H"
//#include "problem.H"
//#include "normalisation.H"



class CoolingSolver: public eulersolver {
public:
	CoolingSolver(Data &gdata);
	~CoolingSolver();
	REAL singlestep(Data &, gridFunc &, ProblemType &);
	void init(Data &, gridFunc &, ProblemType &);
	int test_me();
private:
	/** Implementation of Euler scheme with adaptive time-step size
	 * */
	REAL singlestep_adaptiveEuler(Data &, gridFunc &, ProblemType &);
	/** Implementation of exact integration scheme from Townsend (2009)
	 * */
	REAL singlestep_exactIntScheme(Data &, gridFunc &, ProblemType &);

	REAL get_Tnext_adaptiveEuler(REAL rho_num, REAL T_old_num, REAL del_t);

	REAL get_Tnext_exactIntScheme(REAL rho_num, REAL T_old_num, REAL del_t);

	/** Implementation of function Y given in Eq. (A5) in Townsend (2009)
	 * */
	REAL get_Y(REAL Temp);

	/** Implementation of inverse of function Y given in Eq. (A7) in Townsend (2009)
	 * */
	REAL get_Yinv(REAL Temp);
	CoolingFunction *fCool;

	/** Prepare exact solution scheme by Townsend (2009)
	 * */
	void prepare_exactSolver();

	/** Find discrete interval where current temperature is located -- needed for exact solver
	 * */
	int get_tempInterval(double Temp_num);
	int find_Interval(NumArray<double> &val_arr, double val);

	NumArray<double> powerLawLambda, Yk_bin, discTemp_num;
	NumArray<double> discTemp_phys, discLam_phys;

	/** Reference temperature for implementation by Townsend (2009)
	 * */
	double T_ref;
	/** Lower boundary for temperature
	 * */
	double T_min_num;
	int q_rho, q_Eges, q_wind1, q_wind2;
	/** Switch for type of solver
	 * 1 -> adaptive Euler
	 * 2 -> exact integration scheme by Townsend (2009)
	 * */
	int solver_type;
	int ionisation_index;
};

#endif
