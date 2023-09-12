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
#include "limiter.H"
#include "problem.H"
#include "eos.H"
#include "transformations.H"
#include "sources.H"
#include "timestepping.H"
#include "RiemannSolver.H"
#include "RiemannSolverHD.H"
#include "PhysFluxes.H"
#include "PhysFluxesHD.H"
#include "reconst_block.H"
#include "queue.H"


class HyperbolicSolver: public rksolver {
public:
	HyperbolicSolver(Data &, ProblemType &);
	~HyperbolicSolver();
	virtual void init_constants(Data &);
	virtual void addConstants_toH5(Data &, Hdf5Stream &h5out);
	virtual void getConstants_fromH5(Data &, Hdf5iStream &h5in);
	virtual void init(Queue &, Data &, gridFunc &, ProblemType &);
	virtual void restart(Queue &, Data &, gridFunc &, ProblemType &);
	virtual double singlestep(Data &, gridFunc &, ProblemType &, int, Queue&);
	virtual void phystest(Data &, gridFunc &, ProblemType &, int, int iFluid=0);

private:
	void init_general(Queue &, Data &, gridFunc &, ProblemType &);
	double compute_cfl(Data &, ProblemType &, double, double, int);
	void gen_MagBoundValues(Data &, gridFunc &, ProblemType &);
	/**Compute magnetic field from vector potential
	 * @PARAM initial_setup indicated if computation is done for the first time ->
	 * to be used e.g. by TransPot2Mag
	 * */
	void compute_B(Data &, gridFunc &, ProblemType &, bool initial_setup=false);
	double compute_divB(Data &, gridFunc &, ProblemType &);
	double compute_divB(Data &, gridFunc &, ProblemType &, NumMatrix<double,3> &);
	void CheckNan(NumMatrix<double,3> &, int, int, int, string);
	void TransPot2Mag(Data &, gridFunc &, ProblemType &);
	void TransMag2Pot(Data &, gridFunc &, ProblemType &);
	void ComputeChecksumA(Data &, double []);
	double power(double, int);
	inline double pressure(const double &);
	inline double pressure(double, double, double, double, double);
	inline double pressure(Data &, ProblemType &, const double &,
	                     const double &, const double &, const double &);
	void UserEquations(const Data &);
	void set_UserPde(const Data &);
	void set_TimeIntegrator(Queue &, const Data &, gridFunc &);

	std::vector<std::unique_ptr<TimeIntegrator>> TimeIntegratorGeneric;
	std::vector<std::unique_ptr<TimeIntegrator>> TimeIntegratorUser;

	//	int Fluid_Type;

	double yama;
	double rho0;
	double c2_iso, denomPres;
	double Mass0, Eges0, mv0[3], BAve0[3];
	double ChecksumA[3];
	double rhomin, Ethmin, vmax;
	double Ekin_tave, Emag_tave, Etherm_tave, Eges_tave;
	double Ekfluc_tave, Ebfluc_tave, Vorticity_tave, Enstrophy_tave;
	double cs_tave, v2ave_tave;
	double time_ave;

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
	std::unique_ptr<EquationOfState> eos;
	std::unique_ptr<SourceTerms> sources;
	std::unique_ptr<Transformations> Trafo;
	NumMatrix<double,3> *v_max_p;
	NumMatrix<double,3> *v_max_m;
	// For the new solver:
	std::vector<std::unique_ptr<RiemannSolver>> Riemann;
	
	//std::vector<std::unique_ptr<Reconstruction>> ReconstX, ReconstY, ReconstZ;
	std::unique_ptr<PhysFluxes> PhysFlux;
	//Saves *Save;

	/** Buffer for HyperbolicSolver::get_Changes */
	NumMatrix<double,1> gflux[DIM];

	int n_om, n_omInt, n_omIntAll, n_Omega;
	int n_omUser, n_omIntUser, n_OmegaUser;
	bool with_mag;	
};

#endif