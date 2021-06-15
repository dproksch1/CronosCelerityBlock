#include "data.H"
#include "gridgen.H"
#include "problem.H"
// // #include "mod_teststar3D.H"
// #include "mod_CAK2.H"
// #include "mod_CWB.H"
#include "mod_linearWaves.H"
#include "mod_EigenWaves.H"

void Environment::setType(Data &gdata) {
	int type = static_cast<int>(value((char*)"type"));

	// if (type == 3) {
	// 	Problem = new teststar3D(gdata);
	// } else if (type == 26) {
	//   Problem = new mod_CAK2(gdata);
	// } else if (type == 222) {
	//   Problem = new CWBVar(gdata);
	// } else
	if (type == 153) {
	  Problem = new LinearWaveDecay(gdata);
	} else if (type == 163) {
	  Problem = new EigenWavesDecay(gdata);
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
