#include "data.H"
#include "mod_Sod.H"

void Environment::setType(Data &gdata) {
	int type = static_cast<int>(value((char*)"type"));

	if (type == 99) {
		Problem = std::make_unique<ShockTubeTestHydro>(gdata);
	} else {
		cerr << "   Unknown problem type: " << type <<" -- exiting " << endl;
		exit(-22);
	}

	Problem->init_problem(gdata);
	Problem->set_InfoData(gdata);

}

// In case of multifluid - add setup file
#if(FLUID_TYPE == CRONOS_MULTIFLUID)
#include "multifluid_setup.C"
#endif
