#include "RiemannSolverHD.H"
#include "specific.H"

void HyperbolicSolver::UserEquations(const Data &gdata) {

	if(gdata.rank == 0) {
		cout << " Adding solver for user fields " << endl;
	}

	// Making instances for 1D fields
	fieldsXUser = new fields_1D(gdata, 0, N_OMINT_USER);
	fieldsYUser = new fields_1D(gdata, 1, N_OMINT_USER);
	fieldsZUser = new fields_1D(gdata, 2, N_OMINT_USER);

	physValxLUser = new phys_fields_1D(gdata, 0, N_OMINT_USER);
	physValxRUser = new phys_fields_1D(gdata, 0, N_OMINT_USER);
	physValyLUser = new phys_fields_1D(gdata, 1, N_OMINT_USER);
	physValyRUser = new phys_fields_1D(gdata, 1, N_OMINT_USER);
	physValzLUser = new phys_fields_1D(gdata, 2, N_OMINT_USER);
	physValzRUser = new phys_fields_1D(gdata, 2, N_OMINT_USER);

	ReconstXUser = new Reconstruction*[TIME_SUBSTEPS];
	ReconstYUser = new Reconstruction*[TIME_SUBSTEPS];
	ReconstZUser = new Reconstruction*[TIME_SUBSTEPS];

	for (int i=0; i < TIME_SUBSTEPS; i++) {
		ReconstXUser[i] = new Reconstruction(gdata, 0, N_OMINT_USER, i);
		ReconstYUser[i] = new Reconstruction(gdata, 1, N_OMINT_USER, i);
		ReconstZUser[i] = new Reconstruction(gdata, 2, N_OMINT_USER, i);
	}

	//RiemannXUser = new HLLSolver_gen(gdata, 0, CRONOS_USER);
	//RiemannYUser = new HLLSolver_gen(gdata, 1, CRONOS_USER);
	//RiemannZUser = new HLLSolver_gen(gdata, 2, CRONOS_USER);

}


