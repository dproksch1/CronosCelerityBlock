#ifndef CRONOS_FIELDS_1D
#define CRONOS_FIELDS_1D = 1

#include "data.H"
#include "matrix.H"
#include "constants.H"
#include "queue.H"
#include <vector>

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

#endif