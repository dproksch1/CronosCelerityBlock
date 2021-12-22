// Case without dual energy
#if (CRSWITCH_DUAL_ENERGY == CRONOS_OFF)

#define AUX_ENERGY NONE

#if (ENERGETICS == ISOTHERMAL)
#if (FLUID_TYPE == CRONOS_MHD)
const int N_OMINT = 7;
#elif (FLUID_TYPE == CRONOS_HYDRO)
const int N_OMINT = 4;
#endif

#elif (ENERGETICS == FULL)
#if (FLUID_TYPE == CRONOS_MHD)
const int N_OMINT = 8;
#elif (FLUID_TYPE == CRONOS_HYDRO)
const int N_OMINT = 5;
#endif
#endif // ENERGETICS

// Case with dual energy
#else // CRSWITCH_DUAL_ENERGY -> From here: CRONOS_ON

#define AUX_ENERGY ENTROPY

#if (ENERGETICS == ISOTHERMAL)
#if (FLUID_TYPE == CRONOS_MHD)
const int N_OMINT = 8;
#elif (FLUID_TYPE == CRONOS_HYDRO)
const int N_OMINT = 5;
#endif

#elif (ENERGETICS == FULL)
#if (FLUID_TYPE == CRONOS_MHD)
const int N_OMINT = 9;
#elif (FLUID_TYPE == CRONOS_HYDRO)
const int N_OMINT = 6;
#endif
#endif

#endif // End test for dual energy

//// q-values for different fields
//const int Q_RHO = 0;
//const int Q_Sx = 1;
//const int Q_Sy = 2;
//const int Q_Sz = 3;
//#if (FLUID_TYPE == CRONOS_MHD)
//const int Q_Bx = 4;
//const int Q_By = 5;
//const int Q_Bz = 6;
//const int Q_Eges = 7;
//const int Q_Eadd = 8;
//const int N_SUBS = 3;     // To hold alternative quantities (e.g. vec pot)
//const int N_SUBS = 3;     // To hold alternative quantities (e.g. vec pot)
//#elif (FLUID_TYPE == CRONOS_HYDRO)
//const int Q_Eges = 4;
//const int Q_Eadd = 5;
//const int Q_Bx = -42;
//const int Q_By = -42;
//const int Q_Bz = -42;
//const int N_SUBS = 0;
//#endif


#if (FLUID_TYPE == CRONOS_MHD)
const int N_SUBS = 3;     // To hold alternative quantities (e.g. vec pot)
#elif (FLUID_TYPE == CRONOS_HYDRO)
const int N_SUBS = 0;
#endif


const int N_OMEGA = N_OMINT+N_ADD+N_SUBS;
const int N_OM = N_OMINT+N_ADD+N_SUBS;



// Default value for OMS_USER
#ifndef OMS_USER
#define OMS_USER CRONOS_OFF
//const int N_OMINT_USER = 0;
#endif

//#if (OMS_USER == CRONOS_OFF)
//const int N_OMINT_USER = 0;
//const int N_OM_USER = 0;
//#endif

// If user fields are used determine number of arrays
#if (OMS_USER == CRONOS_ON)

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

//// Finally initial values for all
//const int q_rho = Q_RHO;
//const int q_sx = Q_Sx;
//const int q_sy = Q_Sy;
//const int q_sz = Q_Sz;
//const int q_Bx = Q_Bx;
//const int q_By = Q_By;
//const int q_Bz = Q_Bz;
//const int q_Eges = Q_Eges;
//const int q_Eadd = Q_Eadd;
