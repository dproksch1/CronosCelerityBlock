#ifndef CRONOS_H_INCLUDED
#define CRONOS_H_INCLUDED

#include "cronos_switches.H"

// Use linear correction for geometrical sources
#undef GEOM_SOURCE_CORRECTION

// Use angular momentum
#undef ANGULAR_MOMENTUM

// Physical constants:
#include "physical_constants.H"

// Mathematical constants
const double onethird = 1./3.;
const double twothirds = 2./3.;

const char BEEP = '\a';

#endif
