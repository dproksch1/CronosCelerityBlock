#ifndef CRONOS_CONFIG_H
#define CRONOS_CONFIG_H = 1

// Set default value for fluid type
#ifndef FLUID_TYPE
#define FLUID_TYPE CRONOS_MHD
#endif

// Set default value for energetics
#ifndef ENERGETICS
#define ENERGETICS ISOTHERMAL
#endif

// Set default value for use of angular momentum instead of normal momentum
#ifndef USE_ANGULAR_MOMENTUM
#define USE_ANGULAR_MOMENTUM FALSE
#endif

// Set number of ghost cells
#if (FLUID_TYPE == CRONOS_MHD)
const int B = 3; // Default value for number of ghostcells
#elif (FLUID_TYPE == CRONOS_HYDRO)
const int B = 3; // Default value for number of ghostcells
#else // Indicating multifluid type
const int B = 3; // Default value for number of ghostcells
#endif

// # of oms to integrate (formerly N_OMEGA)


// the definition RK_steps implies the usage of a runge kutta integrator
#ifdef RK_STEPS
#ifndef TIME_INTEGRATOR
#define TIME_INTEGRATOR RUNGEKUTTA
#endif
#endif

// Set the number of time integration substeps for the respective scheme
#if (TIME_INTEGRATOR == RUNGEKUTTA)
#define TIME_SUBSTEPS RK_STEPS
#elif (TIME_INTEGRATOR == VANLEER)
#define TIME_SUBSTEPS 2
#endif

#include "cronos_config_singlefluid.H"


// Setting several default values:

// Set variable for constrained transport according to Gardiner & Stone
#if (CT_TYPE == STONE)
#define FLUX_CT = 1
#endif

// Default value for constrained transport
#ifndef STONE_TYPE
#define STONE_TYPE STONE_SIMPLE
#endif

// Default value for pressure extractio
#ifndef EXTRACT_PRESSURE
#define EXTRACT_PRESSURE FALSE
#endif

// Default value for Geometry
#ifndef GEOM
#define GEOM CARTESIAN
#endif

// Default valze for VEC_POT_BCS
#ifndef VEC_POT_BCS
#define VEC_POT_BCS FALSE
#endif


#endif // End inclusion of file

#define COUT_LOCK do {} while (0)
#define COUT_UNLOCK do {} while (0)

