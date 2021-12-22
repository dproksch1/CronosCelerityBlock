#ifndef CRONOS_FIELDS_1D
#define CRONOS_FIELDS_1D = 1

#include <matrix.H>
#include <vector>

#include "data.H"
#include "constants.H"

class phys_fields_0D {
public:
	phys_fields_0D(const Data &, int, const CronosFluid &fluid);
	phys_fields_0D() {}
	~phys_fields_0D() {};
	bool isgeneric();
	int get_num() const;

	NumArray<REAL> uPri;
	NumArray<REAL> uCon;
	NumArray<REAL> flux_phys;
#if (USE_COROTATION == CRONOS_ON)
	NumArray<REAL> uInertial;
#endif

	double ptotal, ptherm;
	int carbuncle_flag;
private:
	int face, num;
};


class num_fields_0D {
public:
	num_fields_0D(const Data &, const CronosFluid &fluid);
	num_fields_0D() {};
	~num_fields_0D() {};
	REAL v_ch_p;
	REAL v_ch_m;
	NumArray<REAL> flux_num;
	double ptotal_num;
	int carbuncle_flag;
private:
	int num;
};

#endif