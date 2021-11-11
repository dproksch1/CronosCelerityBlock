#ifndef CRONOS_RIEMANNSOLVERHD_H
#define CRONOS_RIEMANNSOLVERHD_H

#include "data.H"
#include "problem.H"
#include "vector.H"
#include "fields_local.H"
#include <stdlib.h>

#define HLLCSOLVER_HYDRO_VEPS 1.e-120

void get_vChar2(const Data &, const ProblemType &, const phys_fields_0D &, const phys_fields_0D &,
			num_fields_0D &, int dir, REAL &);

void get_NumFlux2(const Data &, const phys_fields_0D &, const phys_fields_0D &,
			num_fields_0D &, int dir, int iFluid=0);

#endif
