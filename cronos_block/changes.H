#ifndef CHANGES_H
#define CHANGES_H = 1

#include "save_data.H"

void get_Changes2(const Data &, const num_fields_0D &, const num_fields_0D &,
			NumMatrix<REAL,3> [], int, int, int, int dir, int, int iFluid=0);

REAL get_Changes1DCart2(const Data &, const NumArray<REAL> &, const NumArray<REAL> &,
			int, int, int);

#endif