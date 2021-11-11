#ifndef CRONOS_PHYS_FLUXES_HD_H
#define CRONOS_PHYS_FLUXES_HD_H = 1

#include "PhysFluxes.H"

class PhysFluxesHD: public PhysFluxes {
public:
	PhysFluxesHD(const Data &, const CronosFluid &fluid);
	~PhysFluxesHD(){};
	void get_PhysFlux(const Data &gdata, const ProblemType &Problem, phys_fields_0D &fields,
		        int ix, int iy, int iz, int face, int iFluid=0);
	void get_PhysEmfs(Data &, ProblemType &, cronos::vector<REAL> &,
	                  phys_fields_2D &, const int &) {}
};

#endif
