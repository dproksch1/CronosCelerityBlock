#ifndef CHANGES_H
#define CHANGES_H = 1

//#include "save_data.H"
#include "fields_local.H"
#include "fluid.H"


void get_Changes(const Data &, const num_fields_0D &, const num_fields_0D &,
			NumMatrix<double,3> [], int, int, int, int dir, int, int iFluid=0);

double get_Changes1DCart(const Data &, const NumArray<double> &, const NumArray<double> &,
			int, int, int);

#endif