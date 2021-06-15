#ifndef CRONOS_CONSTANTS_H_INCLUDED
#define CRONOS_CONSTANTS_H_INCLUDED

#include "util.H"
#include "cronos.H"

#define FLUID_TYPE CRONOS_MHD
#define MHD
//#define CT_TYPE STONE
#define STONE_TYPE STONE_CENTRE
#define CT_TYPE CONSISTENT

#define EXTRACT_PRESSURE FALSE 

#define RK_STEPS 2
#define GEOM 1
#define DIM 3
#undef SAVEMEM

#define CRONOS_MOVIE CRONOS_OFF

#define USE_ANGULAR_MOMENTUM FALSE 

#define NON_LINEAR_GRID CRONOS_OFF

#define REAL double 
#define SOUND
#undef PHYSDISS
//const int B = 3;
const int BOUT_FLT = 0;
const int BOUT_DBL = 1;
#define VEC_POT_BCS FALSE 

#define CRSWITCH_DUAL_ENERGY CRONOS_OFF

#define ENERGETICS FULL 
const int N_ADD = 2;      // Used for div(B)


#define OMS_USER FALSE 
//const int N_OMINT_USER = 0;

const int N_P = 0;

#include "cronos_config.H"

#endif


