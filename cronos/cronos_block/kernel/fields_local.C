
#include "fields_local.H"
#include <iostream>

phys_fields_0D::phys_fields_0D(const Data &gdata, int face, const CronosFluid &fluid) {

#if (FLUID_TYPE == CRONOS_MULTIFLUID)
//	fieldLists List(fluid, dir);
	this->num = fluid.get_N_OMINT();
#else
	this->num = N_OMINT;
#endif

	this->face = face;

	uPri.resize(num);
	uCon.resize(num);
	flux_phys.resize(num);
#if (USE_COROTATION == CRONOS_ON)
	uInertial.resize(num);
#endif

	ptotal = 0.;
	ptherm = 0.;
	carbuncle_flag = 0;

}

int phys_fields_0D::get_num() const {
	return num;
}

num_fields_0D::num_fields_0D(const Data &, const CronosFluid &fluid) {
#if (FLUID_TYPE == CRONOS_MULTIFLUID)
//	fieldLists List(fluid, dir);
	this->num = fluid.get_N_OMINT();
#else
	this->num = N_OMINT;
#endif

	v_ch_p = 1.;
	v_ch_m = 1.;
	flux_num.resize(num);

}