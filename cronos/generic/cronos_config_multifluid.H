// Case without dual energy
#if (CRSWITCH_DUAL_ENERGY == CRONOS_OFF)

#define AUX_ENERGY NONE

#if (ENERGETICS == ISOTHERMAL)
const int N_OMINT = 7;
#elif (ENERGETICS == FULL)
const int N_OMINT = 8;
#endif // ENERGETICS

// Case with dual energy
#else // CRSWITCH_DUAL_ENERGY -> From here: CRONOS_ON

#define AUX_ENERGY ENTROPY

#if (ENERGETICS == ISOTHERMAL)
const int N_OMINT = 8;

#elif (ENERGETICS == FULL)
const int N_OMINT = 9;
#endif

#endif // End test for dual energy

const int N_SUBS = 3;     // To hold alternative quantities (e.g. vec pot)



const int N_OMEGA = N_OMINT+N_ADD+N_SUBS;
const int N_OM = N_OMINT+N_ADD+N_SUBS;


// Default value for OMS_USER
#ifndef OMS_USER
#define OMS_USER FALSE
const int N_OMINT_USER = 0;
#endif

// If user fields are used determine number of arrays
#if (OMS_USER == TRUE)

const int N_OMEGA_USER = N_OMINT_USER;
const int N_OM_USER = N_OMINT_USER;

#else

const int N_OMINT_USER = 0;

#endif

// Default value for inclusion of movies
#ifndef CRONOS_MOVIE
#define CRONOS_MOVIE CRONOS_ON
#endif

const int N_OMINT_ALL = N_OMINT + N_OMINT_USER;

