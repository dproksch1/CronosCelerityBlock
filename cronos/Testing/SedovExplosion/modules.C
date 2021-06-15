#include "data.H"
#include "gridgen.H"
#include "problem.H"
#include "mod_Sedov2D.H"
#include "mod_Sedov3D.H"
#include "mod_Sedov3DSpherical.H"

void Environment::setType(Data &gdata) {
	int type = static_cast<int>(value((char*)"type"));

	if (type == 84) {
		Problem = new Sedov3D(gdata);
	} else if (type==83) {
		Problem = new Sedov2D(gdata);
	} else if (type==86) {
		Problem = new Sedov3DSpherical(gdata);
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


// In case of multifluid - also 
#if(FLUID_TYPE == CRONOS_MULTIFLUID)
#include "multifluid_setup.C"
#endif

