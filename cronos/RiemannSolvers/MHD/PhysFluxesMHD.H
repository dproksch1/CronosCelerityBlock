#ifndef CRONOS_PHYS_FLUXES_MHD_H
#define CRONOS_PHYS_FLUXES_MHD_H = 1

#include "PhysFluxes.H"

class PhysFluxesMHD: public PhysFluxes {
public:
	PhysFluxesMHD(const Data &, const CronosFluid &);
	~PhysFluxesMHD(){};
	void get_PhysFlux(Queue /*queue*/, Data &, ProblemType &, cronos::vector<REAL> &,
	                  phys_fields_1D &, const int &, REAL, int iFluid=0);
	void get_PhysFlux(Queue /*queue*/, Data &gdata, ProblemType &Problem, phys_fields_0D &fields,
			int ix, int iy, int iz, int face, int iFluid=0);
	void get_PhysEmfs(Data &, ProblemType &, cronos::vector<REAL> &,
	                  phys_fields_2D &, const int &);
private:
	int DQ_SB;
	int q_Bx, q_By, q_Bz;
};

#endif
