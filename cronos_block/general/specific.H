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
	HyperbolicSolver(Data &, ProblemType &);//used
	~HyperbolicSolver();
	virtual void init_constants(Data &);//used
	virtual void addConstants_toH5(Data &, Hdf5Stream &h5out);//not_used (keep)
	virtual void getConstants_fromH5(Data &, Hdf5iStream &h5in);//not_used (keep)
	virtual void init(Data &, gridFunc &, ProblemType &);//used
	virtual void restart(Data &, gridFunc &, ProblemType &);//not used (could be eliminated, maybe keep)
	virtual REAL singlestep(Data &, gridFunc &, ProblemType &, int, Queue&);
	virtual void phystest(Data &, gridFunc &, ProblemType &, int, int iFluid=0);//used

private:
	void init_general(Data &, gridFunc &, ProblemType &);//used
	REAL compute_cfl(Data &, ProblemType &, REAL, REAL, int);
	void gen_MagBoundValues(Data &, gridFunc &, ProblemType &);//not used (keep if you keep restart)
	/**Compute magnetic field from vector potential
	 * @PARAM initial_setup indicated if computation is done for the first time ->
	 * to be used e.g. by TransPot2Mag
	 * */
	void compute_B(Data &, gridFunc &, ProblemType &, bool initial_setup=false);//not used
	REAL compute_divB(Data &, gridFunc &, ProblemType &);//not used
	REAL compute_divB(Data &, gridFunc &, ProblemType &, NumMatrix<double,3> &);//not used
	void CheckNan(NumMatrix<double,3> &, int, int, int, string);//used
	void TransPot2Mag(Data &, gridFunc &, ProblemType &); //not used
	void TransMag2Pot(Data &, gridFunc &, ProblemType &); //not used
	void ComputeChecksumA(Data &, REAL []); //not used
	REAL power(double, int); //not used
	inline REAL pressure(const REAL &);
	inline REAL pressure(REAL, REAL, REAL, REAL, REAL);
	inline REAL pressure(Data &, ProblemType &, const REAL &,
	                     const REAL &, const REAL &, const REAL &);
	void UserEquations(const Data &);
	void set_UserPde(const Data &);
	void set_TimeIntegrator(const Data &, gridFunc &); //used

	std::vector<std::unique_ptr<TimeIntegrator>> TimeIntegratorGeneric;
	std::vector<std::unique_ptr<TimeIntegrator>> TimeIntegratorUser;

	//	int Fluid_Type;

	REAL yama;
	REAL rho0;
	REAL c2_iso, denomPres;
	REAL Mass0, Eges0, mv0[3], BAve0[3];
	REAL ChecksumA[3];
	REAL rhomin, Ethmin, vmax;
	REAL Ekin_tave, Emag_tave, Etherm_tave, Eges_tave;
	REAL Ekfluc_tave, Ebfluc_tave, Vorticity_tave, Enstrophy_tave;
	REAL cs_tave, v2ave_tave;
	REAL time_ave;

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
	NumMatrix<REAL,3> *v_max_p;
	NumMatrix<REAL,3> *v_max_m;
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