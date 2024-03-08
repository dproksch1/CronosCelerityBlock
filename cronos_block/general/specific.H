#ifndef CRONOS_SPECIFIC_H
#define CRONOS_SPECIFIC_H = 1

#include "Hdf5File_cbase.H"
#include "solver.H"
#ifdef _MSC_VER
#include "timewrapper.H"
#else
#include <time.h>
#include <sys/time.h>
#endif
#include <string>
#include <vector>

#include "vector.H"
#include "CException.H"
#include "problem.H"
#include "transformations_block.H"
#include "timestepping.H"
#include "RiemannSolver.H"
#include "reconst_block.H"
#include "queue.H"
#include "lazy_assertion.H"


class HyperbolicSolver: public RKSolver {
public:
	HyperbolicSolver(Data &, ProblemType &);
	~HyperbolicSolver();
	virtual void init_constants(Data &);
	virtual void addConstants_toH5(Data &, Hdf5Stream &h5out);
	virtual void getConstants_fromH5(Data &, Hdf5iStream &h5in);
	virtual void init(Queue &, Data &, GridFunc &, ProblemType &);
	virtual void restart(Queue &, Data &, GridFunc &, ProblemType &);
	virtual double singlestep(Data &, GridFunc &, ProblemType &, int, Queue&);
	// virtual void phystest(Data &, GridFunc &, ProblemType &, int, int iFluid=0);

private:
	void init_general(Queue &, Data &, GridFunc &, ProblemType &);
	double compute_cfl(Data &, ProblemType &, double, double, int);
	void gen_MagBoundValues(Data &, GridFunc &, ProblemType &);
	/**Compute magnetic field from vector potential
	 * @PARAM initial_setup indicated if computation is done for the first time ->
	 * to be used e.g. by TransPot2Mag
	 * */
	double compute_divB(Data &, GridFunc &, ProblemType &);
	double compute_divB(Data &, GridFunc &, ProblemType &, NumMatrix<double,3> &);
	double power(double, int);
	inline double pressure(const double &);
	inline double pressure(double, double, double, double, double);
	inline double pressure(Data &, ProblemType &, const double &,
	                     const double &, const double &, const double &);
	void UserEquations(const Data &);
	void set_UserPde(const Data &);
	void set_TimeIntegrator(Queue &, const Data &, GridFunc &);

	std::vector<std::unique_ptr<TimeIntegrator>> TimeIntegratorGeneric;
	std::vector<std::unique_ptr<TimeIntegrator>> TimeIntegratorUser;

	//	int Fluid_Type;

	double yama;
	double rho0;
	double c2_iso, denomPres;
	double rhomin, Ethmin, vmax;

	int thermal;
	int TempNeg;
	int debug;
	int q_rho, q_sx, q_sy, q_sz, q_Bx, q_By, q_Bz, q_Eges, q_Eadd;

	bool IntegrateA, bcVecPotResized;
	bool is_OldRestart;

	clock_t cstart, cstep, cend;
	timeval tick,tock,tock2;
	timeval tstep;
	timeval tin1, tin2, tin3, tin4, tin5, tin6, tin7, tin8, tin9, tin10;
	timeval tin11, tin12, tin13, tin14;
	std::unique_ptr<Transformations_Block> Trafo;
	// For the new solver:
	std::vector<std::unique_ptr<RiemannSolver>> Riemann;
	
	//std::vector<std::unique_ptr<Reconstruction>> ReconstX, ReconstY, ReconstZ;
	// std::unique_ptr<PhysFluxes> PhysFlux;
	//Saves *Save;

	/** Buffer for HyperbolicSolver::get_Changes */
	// NumMatrix<double,1> gflux[DIM];

	int n_om, n_omInt, n_omIntAll, n_Omega;
	int n_omUser, n_omIntUser, n_OmegaUser;
	bool with_mag;	
};

#endif