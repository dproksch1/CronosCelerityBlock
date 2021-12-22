#ifndef CRONOS_PROBLEM_H
#define CRONOS_PROBLEM_H = 1

#include <string>
#include <map>
#include "matrix.H"
#include "util.H"
#include "data.H"
#include "constants.H"
#include "movie.H"

#include "vector.H"
#include "fields_local.H"
#include "normalisation.H"
#include "units.H"
#include "physical_constants.H"

using namespace std;


class ProblemType {
public:
	ProblemType(const Data &gdata);
	virtual void init_fields(Data &, int [], int []) = 0;
	// virtual void init_fields(Data &, int q, int [], int []);
	// virtual REAL init_fields(Data &, int q, REAL xPos, REAL yPos, REAL zPos);
	void init_fields(Data &);
	void init_fields(Data &, int);
	virtual void init_problem(Data &) {};
	virtual string get_Name();   
	int get_Type();
	virtual void bc_User(Data &, NumMatrix<REAL,3> &, int, int, int, int);
	virtual void bc_periodic(Data &,
	                         NumMatrix<double,3> &, NumMatrix<double,3> &,
	                         int, int, bool){};
	virtual void bc_Flux(Data &, NumMatrix<REAL,3> [],
	                     NumMatrix<REAL,3> [], NumMatrix<REAL,3> []) {}
	virtual void bc_Flux(Data &, NumMatrix<REAL,3> [], int) {}
	virtual void bc_emf(Data &, NumMatrix<double,3> &,
	                    NumMatrix<double,3> &, NumMatrix<double,3> &) {}
	virtual REAL c2_iso(const Data &, const REAL &, const REAL &, const REAL &) const;
	virtual void src_User(Data &, NumMatrix<REAL,3> [], NumMatrix<REAL,3> []) {}
	virtual void src_Add(Data &, bool [N_OMINT]) {}
	virtual REAL eta(Data &, const REAL &, const REAL &, const REAL &);
	virtual REAL nu(Data &, const REAL &, const REAL &, const REAL &);
	virtual void computePhystest(Data &){}
	virtual void writePhystestInfo(Data &){}
	virtual void writePhystest(Data &, ofstream &){}
	virtual void computeFluct(Data &, double &, double &);
	virtual int checkConvergence(Data &);
	virtual void writeMovie(const Data &gdata, Movie & mov) const ;
	virtual bool force_max(const int &);
	virtual REAL max_Val(const int &);
	virtual bool force_min(const int &);
	virtual REAL min_Val(const int &);
	virtual REAL grid_user_x(REAL ratio);
	virtual REAL grid_user_y(REAL ratio);
	virtual REAL grid_user_z(REAL ratio);
	virtual void ReadFromH5(Hdf5iStream &) {}
	virtual void WriteToH5(Hdf5Stream &) {}
	// Allow a change of the User fields between conservative and
	// characteristic form
	virtual void get_Cons(Data &, cronos::vector<REAL> &Pos,
	                      phys_fields_1D &physField, 
	                      phys_fields_1D &physFieldsUser,
	                      int dir, REAL shift);
	// Routines replaced by TransPrim2Cons et al...
	// virtual void TransConsCharUser(Data &) {}
	// virtual void TransCharConsUser(Data &) {}

	virtual void TransPrim2Cons(Data &) {}
	virtual void TransCons2Prim(Data &) {}


	virtual void get_PhysFluxUser(Data &gdata, cronos::vector<REAL> &iPos,
	                              phys_fields_1D &pf, phys_fields_1D &pfUser,
	                              int dir, REAL shift) {}
	virtual void get_PhysFlux(Data &gdata, cronos::vector<REAL> &iPos,
			phys_fields_1D &pf, int dir, REAL shift) {}

	void name_User(Data &);
	/**Routine to access some internal variables from anywhere (e.g. from an Euler solver)*/
	double get_internalValue(int){return 0.;}
	int get_q(NumMatrix<REAL,3> &, string fluidName="fluid00");
	void set_q(int q, NumMatrix<REAL,3> &);
	void set_Info(bool);
	void set_AsciiOut(bool);
	void set_InfoData(Data &);
	bool checkout(int);
	bool get_Info();
	bool get_AsciiOut();
	virtual ~ProblemType(){}

	normalisation * TrafoNorm;

	/** Metallicity -> mass fractions of hydrogen and helium
	 * */
	double massFraction_X, massFraction_Y;
	/**Default mean particle mass for pure, fully-ionised electron-proton plasma
	 * */
	Quantity meanParticleMass=0.5*CRONOS_CONSTANTS::HydrogenMass;
	REAL gamma;
	REAL rho0;
	bool mag, info, asciiOut, infodata[9];
	int q_rho, q_sx, q_sy, q_sz, q_Bx, q_By, q_Bz, q_Eges, q_Eadd;
	int q_Ax, q_Ay, q_Az;
	int n_omInt, n_omIntAll;
protected:
	REAL nu0, eta0;
	REAL cs2;
	string name;
	int type;
#if(FLUID_TYPE==CRONOS_MULTIFLUID)
	int numFluids;
#endif
	map<std::string, int> fieldNames;
};  

#endif
