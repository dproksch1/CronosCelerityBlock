#ifndef CRONOS_PHYS_FLUXES_HD_H
#define CRONOS_PHYS_FLUXES_HD_H = 1

#include "PhysFluxes.H"
#include "utils.H"

namespace gpu {

	template<typename PhysFluxAccessorType, typename UPriAccessorType, typename UConAccessorType, typename SinglePthermAccessorType, typename FluidConstType>
	void get_PhysFlux(PhysFluxAccessorType& physFlux, UPriAccessorType& uPri, UConAccessorType& uCon,
					SinglePthermAccessorType& ptherm, int face, FluidConstType fluidConst) {
		
		int dir = face/2;

		physFlux[fluidConst[fluidConst_q_sx]] = 0.;
		physFlux[fluidConst[fluidConst_q_sy]] = 0.;
		physFlux[fluidConst[fluidConst_q_sz]] = 0.;

		// flux for the density
		physFlux[fluidConst[fluidConst_q_rho]] = uPri[fluidConst[fluidConst_q_sx] + dir] * uPri[fluidConst[fluidConst_q_rho]];
		
		// flux for momentum vector 
#if EXTRACT_PRESSURE == TRUE
		for(int q = fluidConst[fluidConst_q_sx]; q <= fluidConst[fluidConst_q_sz]; ++q) {
			physFlux[q] = (uCon[q] * uPri[dir + fluidConst[fluidConst_q_sx]]);
		}
#else
		
		for(int q = fluidConst[fluidConst_q_sx]; q <= fluidConst[fluidConst_q_sz]; ++q) {
			if(q == dir+1) {
				physFlux[q] += ((uPri[q] * uCon[q]) + ptherm);
			} else {
				physFlux[q] = (uCon[q] * uPri[dir + fluidConst[fluidConst_q_sx]]);
			}
		}
		
#endif

		// flux for total energy
		if(ENERGETICS == FULL){

#if (USE_COROTATION == CRONOS_ON)
			// not implemented atm (dproksch TODO)
			//fields.flux_phys[q_Eges] = fields.uCon[q_Eges]*fields.uPri[q_sx+dir] +
			//		fields.ptherm*fields.uInertial[dir];
#else
			physFlux[fluidConst[fluidConst_q_Eges]] = (uCon[fluidConst[fluidConst_q_Eges]] + ptherm)
						* uPri[fluidConst[fluidConst_q_sx] + dir];
#endif
		}

#if(CRSWITCH_DUAL_ENERGY == CRONOS_ON)
		// Compute physical flux for entropy evolution -- flux is
		// Entropy*v_par
		physFlux[fluidConst[fluidConst_q_Eadd]] = uCon[fluidConst[fluidConst_q_Eadd]] * uPri[fluidConst[fluidConst_q_sx] + dir];
#endif

#if (USE_ANGULAR_MOMENTUM == TRUE)
		// not yet implemented
#endif

	}


}

class PhysFluxesHD: public PhysFluxes {
public:
	PhysFluxesHD(const Data &, const CronosFluid &fluid);
	~PhysFluxesHD(){};
	void get_PhysFlux(const Data &gdata, const ProblemType &Problem, phys_fields_0D &fields,
		        int ix, int iy, int iz, int face, int iFluid=0);
};

#endif
