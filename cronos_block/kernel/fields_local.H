/**
 * @file
 * @brief Field data structures for computations in Non-Celerity Cronos
 *
 * This file provides data structures for fields in Non-Celerity simulations.
 * It is unused for the standalone usage of this version, but has to remain for adequate proxy-app integration.
 */

#ifndef CRONOS_FIELDS_1D
#define CRONOS_FIELDS_1D = 1

#include <matrix.H>
#include <vector>

#include "data.H"
#include "constants.H"

//! @brief physical flux tracking for Non-Celerity Cronos
class phys_fields_0D {
public:
	phys_fields_0D(const Data &, int, const CronosFluid &fluid);
	phys_fields_0D() {}
	~phys_fields_0D() {};
	bool isgeneric();
	int get_num() const;

	NumArray<double> uPri;
	NumArray<double> uCon;
	NumArray<double> flux_phys;
#if (USE_COROTATION == CRONOS_ON)
	NumArray<double> uInertial;
#endif

	double ptotal, ptherm;
	int carbuncle_flag;
private:
	int face, num;
};

//! @brief numerical flux tracking for Non-Celerity Cronos
class num_fields_0D {
public:
	num_fields_0D(const Data &, const CronosFluid &fluid);
	num_fields_0D() {};
	~num_fields_0D() {};
	double v_ch_p;
	double v_ch_m;
	NumArray<double> flux_num;
	double ptotal_num;
	int carbuncle_flag;
private:
	int num;
};

#endif