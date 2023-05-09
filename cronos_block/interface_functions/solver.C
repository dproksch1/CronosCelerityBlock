#include "solver.H"

#include "gridgen.H"

#define CR_CELERITY CELERITY_OFF

void Environment::init_solvers(Data &gdata)
{
	// Specified user class
	RKSolver = std::make_unique<HyperbolicSolver>(gdata, *Problem);

	try{
		RKSolver->init_constants(gdata);
	} catch (CException exep) {
		Abort(gdata, exep);
	}

}


void Environment::pdestep(Data &gdata, Queue& queue)
{
	gdata.cfl = 0.;

	// double timeinit(gdata.time);

#ifdef TIME_INTEGRATOR

	for (int n = 0; n < TIME_SUBSTEPS; n++) {
		if(gdata.rank==0 && Problem->checkout(0)){
			cout << "  RKSTEP = " << n+1 << endl;
		}
		try {
			// gdata.cfl = (std::max(RKSolver->singlestep(gdata, *gfunc, *Problem, n, queue), gdata.cfl));
			RKSolver->singlestep(gdata, *gfunc, *Problem, n, queue);
		} catch (CException exep) {

			Abort(gdata, exep);
		}
	}

	// queue.slow_full_sync();

	// auto omRange = gdata.omSYCL[0].get_range();
	// size_t nom_max[3] = {omRange.get(0), omRange.get(1), omRange.get(2)};
	// for (int q = 0; q < N_OMINT; q++) {

	// 	queue.submit(celerity::allow_by_ref, [=, &gdata](celerity::handler& cgh) {
	// 		celerity::accessor omSYCL_acc{gdata.omSYCL[q], cgh, celerity::access::all{}, celerity::read_only_host_task};
	// 		cgh.host_task(celerity::on_master_node, [=, &gdata]{
	// 			for (int i = 0; i < nom_max[0]; i++) {
	// 				for (int j = 0; j < nom_max[1]; j++) {
	// 					for (int k = 0; k < nom_max[2]; k++) {
	// 						gdata.om[q](i-3,j-3,k-3) = omSYCL_acc[i][j][k];
	// 					}
	// 				}
	// 			}
	// 		});
	// 	});
	// }
	// queue.slow_full_sync();

#endif

}
