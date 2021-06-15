#include "data.H"
#include "mod_MultTest.H"

void Environment::setType(Data &gdata) {
	int type = static_cast<int>(value((char*)"type"));

	if (type == 240) {
		Problem = new MultTest(gdata);
	} else {
		cerr << "   Unknown problem type -- exiting " << endl;
		exit(-22);
	}

	Problem->init_problem(gdata);
	Problem->set_InfoData(gdata);

}



#if (NON_LINEAR_GRID == CRONOS_ON)
void Grid::set_UserGridFunction(int dir, REAL beg, REAL end, REAL len) {
//	myNonlinGrid[dir] = new GridFunction_Sin();
	myNonlinGrid[dir] = new GridFunction_modEXP(beg,end,len);
  // myNonlinGrid[dir] = new GridFunction_Lin();

}
#endif

// In case of multifluid - add setup file
#if(FLUID_TYPE == CRONOS_MULTIFLUID)
#include "multifluid_setup.C"
#endif
