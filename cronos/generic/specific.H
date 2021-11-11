#ifndef CRONOS_SPECIFIC_H
#define CRONOS_SPECIFIC_H = 1


#include "Hdf5File_cbase.H"
#include "solver.H"
#ifdef _MSC_VER
#include "timewrapper.H"
#else
#include <time.h>
#endif
#include <string>
#include <vector>
//#include <sys/time.h>
#ifdef _MSC_VER
#include "timewrapper.H"
#else
#include <time.h>
#endif
#include "vector.H"
#include "CException.H"
#include "limiter.H"
#include "problem.H"
#include "eos.H"
#include "transformations.H"
#include "sources.H"
// #include "H5Cpp.h"

// for the new solver
#include "RiemannSolver.H"
#include "RiemannSolverHD.H"
#include "RiemannSolverMHD.H"
#include "reconst.H"
#include "reconst_2nd.H"
#include "reconst_WENO.H"
#include "timestepping.H"
#include "save_data.H"
#include "PhysFluxes.H"
#include "ConstrainedTransport.H"
#include "queue.H"

// #ifndef H5_NO_NAMESPACE
// using namespace H5;
// #endif

class HyperbolicSolver: public rksolver {
public:
	HyperbolicSolver(Data &, ProblemType &);
	~HyperbolicSolver();
	virtual void init_constants(Data &);
	virtual void addConstants_toH5(Data &, Hdf5Stream &h5out);
	virtual void getConstants_fromH5(Data &, Hdf5iStream &h5in);
	virtual void init(Data &, gridFunc &, ProblemType &);
	virtual void restart(Data &, gridFunc &, ProblemType &);
	virtual REAL singlestep(Data &, gridFunc &, ProblemType &, int, Queue&);
	virtual void phystest(Data &, gridFunc &, ProblemType &, int, int iFluid=0);

private:
	void init_general(Data &, gridFunc &, ProblemType &);
	void get_Changes(Data &, NumMatrix<REAL, 1> [], NumMatrix<REAL,3> [],
	                 NumMatrix<REAL,1> &, int, int, int, int, int, int iFluid=0);
	void get_ChangesEmf(Data &, NumMatrix<REAL, 3> [],
	                    NumMatrix<REAL, 2> &, int, int);
	void get_Changes(Data &, NumMatrix<REAL, 1> &, NumMatrix<REAL, 1> &, NumMatrix<REAL,3> [],
	                 NumMatrix<REAL,1> &, int, int, int, int, int iFluid=0);
	void get_Changes(Data &, NumArray<REAL> &, NumArray<REAL> &,
			NumArray<REAL> &, NumArray<REAL> &, NumArray<REAL> &, NumArray<REAL> &,
			NumMatrix<REAL,3> [], double, int, int, int, int, int iFluid=0);
	void get_Changes(const Data &, const num_fields_0D &, const num_fields_0D &,
			NumMatrix<REAL,3> [], int, int, int, int dir, int, int iFluid=0) const;
//	void get_Changes(Data &, NumArray<REAL> &, NumArray<REAL> &,
//			NumMatrix<REAL,3> [], double, double, int, int, int, int dir, int, int iFluid=0);
	void get_ChangesX(Data &, NumMatrix<REAL, 1> &, NumMatrix<REAL, 1> &, NumMatrix<REAL,3> [],
			NumMatrix<REAL,1> &, int, int, int, int, int iFluid=0);
	REAL get_Changes1DCart(const Data &, const NumArray<REAL> &, const NumArray<REAL> &,
			int, int, int) const;
	REAL compute_cfl(Data &, ProblemType &, REAL, REAL, int);
	void gen_MagBoundValues(Data &, gridFunc &, ProblemType &);
	/**Compute magnetic field from vector potential
	 * @PARAM initial_setup indicated if computation is done for the first time ->
	 * to be used e.g. by TransPot2Mag
	 * */
	void compute_B(Data &, gridFunc &, ProblemType &, bool initial_setup=false);
	REAL compute_divB(Data &, gridFunc &, ProblemType &);
	REAL compute_divB(Data &, gridFunc &, ProblemType &, NumMatrix<double,3> &);
	void phystest(Data &, gridFunc &, ProblemType &);
	void CheckNan(NumMatrix<double,3> &, int, int, int, string);
	void CheckNeg(NumMatrix<double,3> &, int, int);
	// inline REAL TransE2Eth(REAL &, REAL &, REAL &, REAL &);
	// inline REAL TransEth2E(REAL &, REAL &, REAL &, REAL &);
	// inline REAL TransT2E(REAL , REAL , REAL , REAL );
	void TransPot2Mag(Data &, gridFunc &, ProblemType &);
	void TransMag2Pot(Data &, gridFunc &, ProblemType &);
	void ComputeChecksumA(Data &, REAL []);
	REAL power(double, int);
	inline REAL pressure(const REAL &);
	inline REAL pressure(REAL, REAL, REAL, REAL, REAL);
	inline REAL pressure(Data &, ProblemType &, const REAL &,
	                     const REAL &, const REAL &, const REAL &);
	void UserEquations(const Data &);
	void set_UserPde(const Data &);
	void set_TimeIntegrator(const Data &, gridFunc &);

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
	//std::unique_ptr<RiemannSolver> RiemannX, RiemannY, RiemannZ;
	std::unique_ptr<fields_1D> fieldsX, fieldsY, fieldsZ;
	//std::shared_ptr<fields_1D> fields[DIR_NUM];
	std::unique_ptr<phys_fields_1D> physValxL, physValxR, physValyL, physValyR;
	std::unique_ptr<phys_fields_1D> physValzL, physValzR;
	std::vector<std::unique_ptr<Reconstruction>> ReconstX, ReconstY, ReconstZ;
	std::unique_ptr<PhysFluxes> PhysFlux;
	Saves *Save;
	ConstrainedTransport *CTSolveX, *CTSolveY, *CTSolveZ;
	// Declarations for the User fields:
	Reconstruction **ReconstXUser, **ReconstYUser, **ReconstZUser;
	phys_fields_1D *physValxLUser, *physValxRUser;
	phys_fields_1D *physValyLUser, *physValyRUser;
	phys_fields_1D *physValzLUser, *physValzRUser;
	PhysFluxesUser *PhysFluxUser;
	RiemannSolver *RiemannXUser, *RiemannYUser, *RiemannZUser;

	/** Buffer for HyperbolicSolver::get_Changes */
	NumMatrix<double,1> gflux[DIM];

	int n_om, n_omInt, n_omIntAll, n_Omega;
	int n_omUser, n_omIntUser, n_OmegaUser;
	bool with_mag;
#if (FLUID_TYPE == CRONOS_MULTIFLUID)
//	RiemannSolver *RiemannSolversX, *RiemannSolversY, *RiemannSolversZ;
//	std::vector<HLLSolver>  RiemannSolversX, RiemannSolversY, RiemannSolversZ;
	std::vector<RiemannSolver*>  RiemannSolversX, RiemannSolversY, RiemannSolversZ;
//	RiemannSolver *RiemannSolversY, *RiemannSolversZ;
//	fields_1D *fieldsMultiX, *fieldsMultiY, *fieldsMultiZ;
	std::vector<fields_1D*> fieldsMultiX, fieldsMultiY, fieldsMultiZ;
//	phys_fields_1D *physValMultixL, *physValMultiyL, *physValMultizL;
//	phys_fields_1D *physValMultixR, *physValMultiyR, *physValMultizR;
	std::vector<phys_fields_1D*> physValMultixL, physValMultiyL, physValMultizL;
	std::vector<phys_fields_1D*> physValMultixR, physValMultiyR, physValMultizR;
//	Reconstruction *ReconstMultiX, *ReconstMultiY, *ReconstMultiZ;
	std::vector<Reconstruction**> ReconstMultiX, ReconstMultiY, ReconstMultiZ;
//	PhysFluxes *PhysFluxMulti;
	std::vector<PhysFluxes*> PhysFluxMulti;
	std::vector<Transformations*> TrafoMulti;
	int numFluids;
//	int N_OM, N_OMINT, N_OMINT_ALL;
//	int N_OM_USER, N_OMINT_USER;
#else
#endif
	fields_1D *fieldsXUser, *fieldsYUser, *fieldsZUser;
	
};

void get_Changes2(const Data &, const num_fields_0D &, const num_fields_0D &,
			NumMatrix<REAL,3> [], int, int, int, int dir, int, int iFluid=0);

REAL get_Changes1DCart2(const Data &, const NumArray<REAL> &, const NumArray<REAL> &,
			int, int, int);


inline REAL HyperbolicSolver::pressure(const REAL &dens)
{
#ifdef SOUND
	if(yama != 1.){
		return c2_iso*pow(dens,yama)*denomPres;
	} else {
		return c2_iso*dens;
	}
#else
	if(yama != 1.){
		return c2_iso*pow(dens,yama)*denomPres*half_beta;
	} else {
		return c2_iso*dens*half_beta;
	}
#endif
}


inline REAL HyperbolicSolver::pressure(Data &gdata, ProblemType &Problem,
                                       const REAL &dens, const REAL &i,
                                       const REAL &j, const REAL &k)
{
#ifdef SOUND
	if(yama != 1.){
		return Problem.c2_iso(gdata,i,j,k)*pow(dens,yama)*denomPres;
	} else {
		return Problem.c2_iso(gdata,i,j,k)*dens;
	}
#else
	if(yama != 1.){
		return Problem.c2_iso(gdata,i,j,k)*pow(dens,yama)*denomPres*half_beta;
	} else {
		return Problem.c2_iso(gdata,i,j,k)*dens*half_beta;
	}
#endif
}


inline REAL HyperbolicSolver::pressure(REAL rhoinv, REAL psq,
                                       REAL Bsq, REAL Energy, REAL yama)
{
	REAL PresVal;
	if (yama == 1){
		PresVal = c2_iso/rhoinv;
	}else{
		PresVal= (yama-1)*(Energy - 0.5*psq*rhoinv - 0.5*Bsq);
	}
	if (PresVal < 0.){
		TempNeg+=1;
		PresVal = 1.e-12;}
	return PresVal;
}


// inline REAL HyperbolicSolver::TransE2Eth(REAL &rhoinv, REAL &psq, REAL &Bsq, REAL &Energy)
// {
// 	REAL ETherm;
// 	ETherm = (Energy - 0.5*psq*rhoinv - 0.5*Bsq);
// 	if (ETherm < 0.){
// 		TempNeg+=1;
// 		ETherm = 1.e-12;
// 	}
// 	return ETherm;
// }


// inline REAL HyperbolicSolver::TransEth2E(REAL &rhoinv, REAL &psq, REAL &Bsq, REAL &ETherm)
// {
// 	REAL Energy;
// 	Energy = (ETherm + 0.5*psq*rhoinv + 0.5*Bsq);
// 	return Energy;
// }


// inline REAL HyperbolicSolver::TransT2E(REAL rhoinv, REAL psq, REAL Bsq, REAL ETherm)
// {
// 	REAL Energy;
// 	ETherm *= 1./(rhoinv*(yama-1.));
// 	Energy = (ETherm + 0.5*psq*rhoinv + 0.5*Bsq);
// 	return Energy;
// }

#endif
