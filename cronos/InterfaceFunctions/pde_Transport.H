#ifndef CRONOS_PDE_TRANSPORT_H
#define CRONOS_PDE_TRANSPORT_H

class PhysFluxUserTransport: public PhysFluxesUser {
public:
	PhysFluxUserTransport(){};
	~PhysFluxUserTransport(){};
	void get_PhysFlux(Data &, ProblemType &, cronos::vector<REAL> &,
	                  phys_fields_1D &, phys_fields_1D &, int, REAL);
	void get_PhysFlux(Data &gdata, ProblemType &Problem,
			phys_fields_0D &fieldsUser, int ix, int iy, int iz, int face);
};


void PhysFluxUserTransport::get_PhysFlux(Data &gdata,
                                         ProblemType &Problem,
                                         cronos::vector<REAL> &iPos,
                                         phys_fields_1D &pf,
                                         phys_fields_1D &pfUser,
                                         int dir, REAL shift){
	pfUser.fluxORIG[0].clear();
	
	// Example: transport of a passive tracer-field
	for (int i = -1; i <= gdata.mx[dir]+1; ++i){
		
		pfUser.fluxORIG[0](i) = pfUser.uPriORIG[0](i)*pf.uPriORIG[1+dir](i);

	}
}


void PhysFluxUserTransport::get_PhysFlux(Data &gdata, ProblemType &Problem,
		phys_fields_0D &fieldsUser, int ix, int iy, int iz, int face) {

	int dir = face/2;
	fieldsUser.flux_phys[0] = 0.;

	// Example: transport of a passive tracer-field

	fieldsUser.flux_phys[0] = fieldsUser.uPri[0]*fieldsUser.uPri[1+dir];

}

#endif
