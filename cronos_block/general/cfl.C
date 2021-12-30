#include "gridgen.H"

using namespace std;

REAL HyperbolicSolver::compute_cfl(Data &gdata, ProblemType &Problem,
                                   REAL cfl_eta, REAL cfl_lin, int n) {

	if(cfl_eta > cfl_lin && n == TIME_SUBSTEPS-1 &&
	   Problem.get_Info() && gdata.rank == 0) {
		cout << " Timestep given by dissipation " << endl;
	}
	REAL cfl = std::max(cfl_eta,cfl_lin);

#ifdef DTCOMP_OLD
	cfl *= gdata.dt;
	if(gdata.rank==0) {
		if(n == TIME_SUBSTEPS-1 && Problem.checkout(1)) {
			cout << "  cfl = " << cfl;
#ifdef PHYSDISS
			cout << " " << cfl_eta;
			cout << " " << comp_max_glob;
#endif
			cout << endl;
		}
	}
	if (cfl > 0.8) {
		cout << " *** cfl > cfl_max (=0.8) *** " << endl;
	}
#endif

	return cfl;

}
