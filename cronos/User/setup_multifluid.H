#ifndef CRONOS_MULTIFLUID_H_INCLUDED
#define CRONOS_MULTIFLUID_H_INCLUDED

#include "util.H"
#include "cronos.H"

// This is an example configuration for a Cronos multifluid project
// Beware: configuration is very different from standard Cronos

const int Mulfifluid_NumFluids = 2;
const int Multifluid_Type[Mulfifluid_NumFluids] = {CRONOS_MHD, CRONOS_HYDRO};
#define OMS_USER FALSE
const int Multifluid_N_OMINT_USER[Mulfifluid_NumFluids] = {0,0};


#endif
