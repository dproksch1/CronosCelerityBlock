#include "PhysFluxesHD.H"

PhysFluxesHD::PhysFluxesHD(const Data &gdata,
		const CronosFluid &fluid) : PhysFluxes(gdata, fluid) {

}

void PhysFluxesHD::get_PhysFlux(const Data &gdata,
	const ProblemType &Problem, phys_fields_0D &fields,
        int /*ix*/, int /*iy*/, int /*iz*/, int face, int iFluid) {

	int dir = face/2;

	fields.flux_phys[q_sx] = 0.;
	fields.flux_phys[q_sy] = 0.;
	fields.flux_phys[q_sz] = 0.;

	// flux for the density
	fields.flux_phys[q_rho] = fields.uPri[q_sx+dir]*fields.uPri[q_rho] ;


	// flux for momentum vector
#if EXTRACT_PRESSURE == TRUE
	for(int q=q_sx; q<=q_sz; ++q) {
		fields.flux_phys[q] = (fields.uCon[q  ]*fields.uPri[dir+q_sx]);
	}
#else
	for(int q=q_sx; q<=q_sz; ++q) {
		if(q == dir+1) {
			fields.flux_phys[q] += (fields.uPri[q]*fields.uCon[q] + fields.ptherm);
		} else {
			fields.flux_phys[q] = (fields.uCon[q]*fields.uPri[dir+q_sx]);
		}
	}
#endif


	// flux for total energy
	if(ENERGETICS == FULL){

#if (USE_COROTATION == CRONOS_ON)
		fields.flux_phys[q_Eges] = fields.uCon[q_Eges]*fields.uPri[q_sx+dir] +
				fields.ptherm*fields.uInertial[dir];
#else
		fields.flux_phys[q_Eges] = (fields.uCon[q_Eges] +
				fields.ptherm)*fields.uPri[q_sx+dir];
#endif

	}

#if(CRSWITCH_DUAL_ENERGY == CRONOS_ON)
	// Compute physical flux for entropy evolution -- flux is
	// Entropy*v_par
	fields.flux_phys[q_Eadd] = fields.uCon[q_Eadd]*fields.uPri[q_sx+dir];
#endif

#if (USE_ANGULAR_MOMENTUM == TRUE)


	int ishift = 0;
	if(face % 2 > 0) ishift = 1;
	int shift_vec[3] = {0,0,0};
	shift_vec[dir] = ishift;

	fields.flux_phys[2] *= gdata.h1(ix, iy, iz,
			shift_vec[0],shift_vec[1],shift_vec[2]);


#endif

}


