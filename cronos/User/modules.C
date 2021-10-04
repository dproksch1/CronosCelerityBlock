#include "data.H"
#include "mod_linHydroWaves.H"
#include "mod_BrioWu.H"
#include "mod_HydroDisc.H"
//#include "mod_Heliosphere.H"
#include "mod_MRICylDisc.H"
#include "mod_SOD.H"
#include "mod_Sedov3D.H"
#include "mod_TestNonlinGrid.H"
#include "mod_MHDBlastPolar.H"
#include "mod_MHDBlast2D.H"
#include "mod_MHDBlast3DPolar.H"
#include "mod_axistestHydro.H"
#include "mod_axistestMHD.H"
#include "mod_axistestVecPot.H"
#include "mod_RotorProblem2D.H"
#include "mod_MHDBlast3D.H"
// #if(FLUID_TYPE==CRONOS_MULTIFLUID)
// #include "mod_MultTest.H"
// #endif

void Environment::setType(Data &gdata) {
	int type = static_cast<int>(value((char*)"type"));

	if (type == 63) {
		Problem = new linHydroWaves(gdata);
	} else if (type == 23) {
		Problem = new MHDBlastPolar(gdata);
	} else if (type == 24) {
		Problem = new MHDBlast3DPolar(gdata);
	} else if (type == 26) {
		Problem = new MHDBlast2D(gdata);
#if(GEOM < 3)
	} else if (type == 27) {
		Problem = new MHDBlast3D(gdata);
#endif
	} else if (type == 33) {
		Problem = new AxistestHydro(gdata);
	} else if (type == 35) {
		Problem = new AxistestMHD(gdata);
#if(GEOM < 3)
	} else if (type == 51) {
		Problem = new RotorProblem2D(gdata);
#endif
	} else if (type == 36) {
		Problem = new AxistestVecPot(gdata);
//	}else if (type == 12){
//		Problem = new Heliosphere(gdata);
	} else if (type == 65) {
		Problem = new HydroDisc(gdata);
	} else if (type == 48) {
		Problem = new MRICylDisc(gdata);
	} else if (type == 92) {
		Problem = new BrioWuTest(gdata);
	} else if (type == 84) {
		Problem = new Sedov3D(gdata);
	} else if (type == 99) {
		Problem = new SOD(gdata);
	} else if (type == 120) {
		Problem = new TestNonlinGrid(gdata);
// #if(FLUID_TYPE==CRONOS_MULTIFLUID)
// 	} else if (type == 240) {
// 		Problem = new MultTest(gdata);
// #endif
	} else {
		cerr << "   Unknown problem type: " << type << " -- exiting " << endl;
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
