#include "data.H"
#include "gridgen.H"
#include "problem.H"

#include "mod_ParkerSpherical1D_iso.H"

void Environment::setType(Data &gdata) {
	int type = static_cast<int>(value((char*)"type"));

	if (type == 65) {
		Problem = new ParkerSph1DIso(gdata);
	} else {
		cerr << "   Unknown problem type -- exiting " << endl;
		exit(-22);
	}

	Problem->init_problem(gdata);
	Problem->set_InfoData(gdata);

}



#if (NON_LINEAR_GRID == CRONOS_ON)
void Grid::set_UserGridFunction(int dir) {
	myNonlinGrid[dir] = new GridFunction_Sin();
  // myNonlinGrid[dir] = new GridFunction_Lin();

}
#endif


// In case of multifluid - add setup file
#if(FLUID_TYPE == CRONOS_MULTIFLUID)
#include "multifluid_setup.C"
#endif
