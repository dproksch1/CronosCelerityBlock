#ifndef CRONOS_PHYS_FLUXES_H
#define CRONOS_PHYS_FLUXES_H = 1

#include "data.H"
#include "fields_local.H"
#include "problem.H"
#include "vector.H"

class PhysFluxes {
public:
	PhysFluxes(const Data &, const CronosFluid &fluid){
		q_rho = fluid.get_q_rho();
		q_sx = fluid.get_q_sx();
		q_sy = fluid.get_q_sy();
		q_sz = fluid.get_q_sz();
		q_Eges = fluid.get_q_Eges();
		q_Eadd = fluid.get_q_Eadd();
	};
  virtual ~PhysFluxes(){};
  virtual void get_PhysFlux(const Data &gdata, const ProblemType &Problem, phys_fields_0D &fields,
	        int ix, int iy, int iz, int face, int iFluid=0) = 0;
protected:
  int q_rho, q_sx, q_sy, q_sz, q_Eges, q_Eadd;
};


class PhysFluxesUser {
public:
	PhysFluxesUser(){};
	virtual void get_PhysFlux(Data &gdata, ProblemType &Problem, phys_fields_0D &fields,
	        int ix, int iy, int iz, int face)=0;
};

#endif
