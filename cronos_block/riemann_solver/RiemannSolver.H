#ifndef CRONOS_RIEMANNSOVLER_H
#define CRONOS_RIEMANNSOVLER_H = 1

#include "data.H"
#include "problem.H"
#include "vector.H"
#include "fields_local.H"
#include <stdlib.h>


/*-----------------------------------------------------

  Base class for Riemann solvers supplying base variables and
  interface

  -----------------------------------------------------*/

class RiemannSolver {
public:
	RiemannSolver(){
		if(value_exists("use_carbuncleFlag")) {
			alpha_carbuncle = value((char*)"alpha_carbuncle");
		}
		verbosity = 0;
	};
	RiemannSolver(const Data &gdata, int /*dir*/, int Fluid_Type) {
		this->Fluid_Type = Fluid_Type;
		verbosity = 0;

		if(value_exists("use_carbuncleFlag")) {
			alpha_carbuncle = value((char*)"alpha_carbuncle");
		}
	};
	virtual ~RiemannSolver(){};
	virtual void get_vChar(const Data &, const ProblemType &, const phys_fields_0D &, const phys_fields_0D &,
			num_fields_0D &,
			int ix, int iy, int iz, int dir, double &cfl_lin) const = 0;
	virtual void get_NumFlux(const Data &, const phys_fields_0D &,
		const phys_fields_0D &, num_fields_0D &, int dir, int iFluid=0) const = 0;
	virtual void reset_Indices(const CronosFluid &) = 0;
	int get_Fluid_Type() {
		return Fluid_Type;
	}
	void compute_carbuncleFlag(Data &gdata);
	void compute_carbuncleFlag(Queue &queue, Data &gdata);
	void set_verbosity(int _verb);
	
protected:
	double norm, alpha_carbuncle;
	int Fluid_Type;
	int verbosity;
	bool use_carbuncleFlag;
	//NumMatrix<double, 1> ptotalL, ptotalR;
};


#endif
