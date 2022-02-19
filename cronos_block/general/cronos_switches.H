#ifndef CRONOS_SWITCHES_H_INCLUDED
#define CRONOS_SWITCHES_H_INCLUDED

// true and false:
#define TRUE	1
#define FALSE	0

// on and off
#define CRONOS_ON 1
#define CRONOS_OFF 0

// Available Geometries:
#define CARTESIAN	1
#define CYLINDRICAL	2
#define SPHERICAL	3

// Available Riemann Solvers
#define HLL		1
#define HLLD		2

// Available Fluid Types
#define CRONOS_HYDRO   1
#define CRONOS_MHD     2
#define CRONOS_MULTIFLUID 3
#define CRONOS_USER   42 

// Available CT-Versions
#define CONSISTENT 	1 // Londrillo & Del Zanna (2000)
#define STONE		2 // Gardiner & Stone

// Available Stone-CT versions
#define STONE_SIMPLE    1
#define STONE_CENTRE    2

// Available Versions of Energy Treatment:
#define ISOTHERMAL  1
#define FULL        2

// Available types of auxiliary energy (for dual energy treatment)
#define ENTROPY 2
//#define NONE -1

// Available Time integration schemes:
#define RUNGEKUTTA 1
#define VANLEER 2

#endif
