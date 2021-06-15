#include "multifluid.H"

void CronosMultifluid::setup() {
	//! Here the user is supposed to suppy all information needed for the multifluid case

	// Here we assume to be working with two fields
	fluidTypes(0) = CRONOS_MHD;
	fluidTypes(1) = CRONOS_HYDRO;

	// Default: no user fields to be included
	for(int iFluid=0; iFluid<numFluids; ++iFluid) {
		userFields(iFluid) = 0;
	}

}
//const int Mulfifluid_NumFluids = 2;
//const int Multifluid_Type[Mulfifluid_NumFluids] = {CRONOS_MHD, CRONOS_HYDRO};
//#define OMS_USER FALSE
//const int Multifluid_N_OMINT_USER[Mulfifluid_NumFluids] = {0,0};
