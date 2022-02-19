#ifndef CRONOS_SOURCES_H
#define CRONOS_SOURCES_H = 1

#include "data.H"
#include "problem.H"
#include "CException.H"
#include "eos.H"
#include "limiter.H"

class SourceTerms {
public:
	SourceTerms(const Data &, const ProblemType &, bool IntegrateA,
	            int thermal);
	~SourceTerms();
	void src_rhs(Data &, ProblemType &, NumMatrix<double,3> nom[]);
	void src_Geom(Data &, ProblemType &, NumMatrix<double,3> nom[]);
	void mod_Geom_Visc_WE(Data &, const ProblemType &, 
	                      double [], double [],
	                      const double &, const double &, const double &, const int &);
	void mod_Geom_Visc_SN(Data &, const ProblemType &, 
	                      double [], double [],
	                      const double &, const double &, const double &, const int &);
	void mod_Geom_Visc_BT(Data &, const ProblemType &, 
	                      double [], double [],
	                      const double &, const double &, const double &, const int &);
	double shift_Geom_WE(Data &, int, int, int);
	double shift_Geom_SN(Data &, int, int, int);
	double shift_Geom_BT(Data &, int, int, int);

	double shift_Geom_WE(Data &, int);
	double shift_Geom_SN(Data &, int);
	double shift_Geom_BT(Data &, int);
private:
	void src_Geom_Hydro(Data &, ProblemType &, NumMatrix<double,3> nom[], const CronosFluid &);
	void src_Geom_Viscosity(Data &, ProblemType &, NumMatrix<double,3> nom[], const CronosFluid &);
	void src_rhs_Ideal(Data &, ProblemType &, NumMatrix<double,3> nom[], const CronosFluid &fluid);
	void src_rhs_Viscosity(Data &, ProblemType &, NumMatrix<double,3> nom[], const CronosFluid &fluid);
	void src_Axis(Data &, NumMatrix<double,3> nom[], const CronosFluid &fluid);

	string geom;
	EquationOfState *eos;
	double facShift;
	bool IntegrateA, thermal;
#ifdef GEOM_SOURCE_CORRECTION
	limiter* Limiter;
#endif
#if (USE_COROTATION == CRONOS_ON)
	double omegaZ;
#endif
};




#endif
