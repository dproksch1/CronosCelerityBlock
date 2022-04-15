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
	virtual void init_fields(Data &, int [], int []) {} /*=*/;
	void init_fields(Data &);
	virtual void init_problem(Data &) {};
	virtual string get_Name(); 
	int get_Type();
	virtual void bc_User(Data &, NumMatrix<double,3> &, int, int, int, int);
	virtual void bc_periodic(Data &,
	                         NumMatrix<double,3> &, NumMatrix<double,3> &,
	                         int, int, bool){};
	virtual void bc_Flux(Data &, NumMatrix<double,3> [],
	                     NumMatrix<double,3> [], NumMatrix<double,3> []) {}
	virtual void bc_Flux(Data &, NumMatrix<double,3> [], int) {}
	virtual void bc_emf(Data &, NumMatrix<double,3> &,
	                    NumMatrix<double,3> &, NumMatrix<double,3> &) {}
	virtual double c2_iso(const Data &, const double &, const double &, const double &) const; //not used
	virtual void src_User(Data &, NumMatrix<double,3> [], NumMatrix<double,3> []) {}
	virtual void src_Add(Data &, bool [N_OMINT]) {}
	virtual double eta(Data &, const double &, const double &, const double &); //not used
	virtual double nu(Data &, const double &, const double &, const double &); //not used
	virtual void computePhystest(Data &){}
	virtual void writePhystestInfo(Data &){}
	virtual void writePhystest(Data &, ofstream &){}
	virtual void computeFluct(Data &, double &, double &);
	virtual int checkConvergence(Data &);
	virtual void writeMovie(const Data &gdata, Movie & mov) const ;
	virtual bool force_max(const int &);
	virtual double max_Val(const int &);
	virtual bool force_min(const int &);
	virtual double min_Val(const int &);
	virtual double grid_user_x(double ratio);
	virtual double grid_user_y(double ratio);
	virtual double grid_user_z(double ratio);
	virtual void ReadFromH5(Hdf5iStream &) {}
	virtual void WriteToH5(Hdf5Stream &) {}

/*	virtual void get_PhysFluxUser(Data &gdata, cronos::vector<double> &iPos,
	                              phys_fields_1D &pf, phys_fields_1D &pfUser,
	                              int dir, double shift) {}
	virtual void get_PhysFlux(Data &gdata, cronos::vector<double> &iPos,
			phys_fields_1D &pf, int dir, double shift) {}*/

	virtual void TransPrim2Cons(Data &) {}
	virtual void TransCons2Prim(Data &) {}

	void name_User(Data &);
	/**Routine to access some internal variables from anywhere (e.g. from an Euler solver)*/
	double get_internalValue(int){return 0.;}
	int get_q(NumMatrix<double,3> &, string fluidName="fluid00");
	void set_q(int q, NumMatrix<double,3> &);
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
	double gamma;
	double rho0;
	bool mag, info, asciiOut, infodata[9];
	int q_rho, q_sx, q_sy, q_sz, q_Bx, q_By, q_Bz, q_Eges, q_Eadd;
	int q_Ax, q_Ay, q_Az;
	int n_omInt, n_omIntAll;
protected:
	double nu0, eta0;
	double cs2;
	string name;
	int type;
	map<std::string, int> fieldNames;
};  

#endif
