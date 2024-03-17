#include "gridgen.H"
#include <float.h>
#include <string.h>
#include "timewrapper.H"
#include <stdio.h>
#include <iostream>
#include <iomanip>
#include <sstream>
// #include "H5Cpp.h"
#include "CException.H"
// #ifndef H5_NO_NAMESPACE
// using namespace H5;
// #endif
#include <filesystem>

#include "queue.H"
#include "lazy_assertion.H"

// local functions:
void RemoveNullsFromString(string &input)
{
	string null = "0";
	while (input[input.length()-1] == null[0]){
		input.erase(input.length()-1);
	}
}


Environment::Environment(Data &gdata)
	: restart_time(0.), restart_step(0),
	  numdbl_out(static_cast<int>(value((char*)"num_bin"))),
	  numdbl_done(0),
	  numdbl_pass(0),
	  numflt_out(static_cast<int>(value((char*)"num_float"))),
	  numflt_done(0),
	  numflt_pass(0),
	  numascii_out(static_cast<int>(value((char*)"num_ascii"))),
	  numascii_done(1),
	  numascii_pass(0),
	  numinfo_out(static_cast<int>(value((char*)"num_info"))),
	  numinfo_done(0),
	  numinfo_pass(1),
	  failnum(0),
	  nummov_done(0),
	  dt_dbl(value((char*)"dt_bin") - 1.e-12),
	  dt_flt(value((char*)"dt_float") - 1.e-12),
	  dt_ascii(value((char*)"dt_ascii") - 1.e-12),
	  dt_info(value((char*)"dt_info") - 1.e-12),
	  dt_mov(value((char*)"dt_mov") - 1.e-12),
	  cfl_set(value((char*)"cfl_threshold")),
	  cfl_min(value((char*)"cfl_min")),
	  cfl_max(value((char*)"cfl_threshold"))
{

	gfunc = std::make_unique<GridFunc>(gdata);
	setup(gdata);

//	rksolver = new RKSolver();
//	EulerSolver = new eulersolver();
	init_solvers(gdata);

	// write out info about norms
	ofstream fout;
	char fname[255];
	//	if(Problem->TrafoNorm && gdata.rank==0) {
	if(Problem->TrafoNorm != NULL && gdata.rank==0) {
		Problem->TrafoNorm->print_norms(cout);
		sprintf(fname,"%s/%s_norms.info",getenv("poub"),getenv("pname"));
		fout.open(fname);
		Problem->TrafoNorm->print_norms(fout);
		fout.close();
	}

	gettimeofday(&tick, 0);

#ifdef DTCOMP_OLD
	dt_debug = static_cast<int>(value((char*)"dt_debug"));
#endif
};


Environment::~Environment() {
//	if(rksolver != NULL) {
//		delete rksolver;
//	}

	//if(Problem != NULL) {
	//	delete Problem;
	//}

	//if(gfunc != NULL) {
	//	delete gfunc;
	//}

	// Delete all new objects:
// 	delete gfunc;

}


//! @brief Performs the initial setup of the simulation environment
void Environment::setup(Data &gdata)
{
	// Set the problem type
	setType(gdata);
  
	// Naming of standard fields:
	Problem->name_User(gdata);
  
#if(FLUID_TYPE == CRONOS_MULTIFLUID)
	int n_omInt = gdata.fluids->get_N_OMINT();
	int n_Omega = gdata.fluids->get_N_OMEGA();
#else
	int n_omInt = gdata.fluid.get_N_OMINT();
	int n_Omega = gdata.fluid.get_N_OMEGA();
#endif

	for(int q=0; q<N_ADD; ++q) {
		std::ostringstream dummy;
		dummy << "Add" << q;
		gdata.om[q+n_omInt].rename(dummy.str());
	}

}


//! @brief Performs the initial writeouts
void Environment::InitOutput(Queue &queue, Data &gdata)
{
	// Initial output

	if(gdata.time == 0.) {
		// Output of movie:
#if (CRONOS_MOVIE == CRONOS_ON)
		WriteMovies(gdata);
#endif

	}

}


/*!
 * @brief Checks for completion using the time and timestep conditions and potentially
 * 			triggers finalization
 */
bool Environment::CheckEnd(Data &gdata, Queue &queue)
{
	// Check if runtime exceeded
	gettimeofday(&tock, 0);
	double delt = ((tock.tv_sec + tock.tv_usec/1.e6) -
	             (tick.tv_sec + tick.tv_usec/1.e6));
	double deltmax = value((char*)"runtime");

	bool RunTimeExit = false;
	if(delt > deltmax) {
		RunTimeExit = true;
	}

	if(RunTimeExit) {
		return Finalize(gdata, queue, static_cast<string>("Runtime exceeded"));
	}

	// Check if time-integration is completed
	bool Completed = false;
	if(gdata.time >= gdata.t_end) {
		Completed = true;
	}

	if(Completed) {
		return Finalize(gdata, queue, static_cast<string>("Time integration completed"));
	}
	return false;
}


/*!
 * @brief Checks for completion using the user-defined convergence rules and potentially
 * 			triggers finalization
 */
bool Environment::CheckEnd_User(Data &gdata, Queue &queue) {
	bool EndProgram_User = Problem->checkConvergence(gdata);

	if(EndProgram_User) {
		return Finalize(gdata, queue, static_cast<string>("User defined stopping criterion reached "));
	}

	return false;
}


//! Computes the timestep delta
void Environment::compute_dt(Data &gdata, double factor)
{

#ifdef DTCOMP_OLD
	// Decreasing timestep-size if cfl > cfl_threshold
	if (gdata.cfl > cfl_max){
		//    cout << " cfl: " << gdata.cfl << " " << cfl_max << endl;
		if(gdata.rank==0 && Problem->checkout(1)){
			cout << "performing half_dt" << endl;
		}
		gdata.dt *= 0.5;   
	}


	// Increasing timestep-size if cfl < clf_min && doubling fits
	// regular time-stepping

	double dtdoub = 2.*gdata.dt;
	double teps =1.e-8*gdata.dt;

#if (CRONOS_MOVIE == CRONOS_ON)
	double tnext_mov = (static_cast<int>(gdata.time/dt_mov)+1)*dt_mov;

	double dttest_mov = abs(fmod((tnext_mov - gdata.time)+teps,dtdoub));
#endif

	double trun = abs(gdata.time - restart_time);
	/*
	  Checking if:
	  (a) cfl < cfl_min
	  (b) time is multiple of dt_doub
	  (c) cfl is really initialised (only if trun > 0.)
	  (d) sufficiently far from restart step
	*/

	if(gdata.cfl < cfl_min && gdata.time > 0. && trun > 1.e-8) {
#if (CRONOS_MOVIE == CRONOS_ON)
		if(dttest_mov < 100*teps) {
#endif
			if(gdata.rank==0 && Problem->checkout(1)) {
				cout << "performing double_dt" << endl;
			}
			gdata.dt *= 2.;
			failnum = 0;
		} else {
			failnum++;
		}
#if (CRONOS_MOVIE == CRONOS_ON)
	}
#endif

	if(dt_debug && failnum > 10) {
		throw CException(" Tried to double dt 10 times without success ");
	}
#else
	// Here the timestep-size is determined by the cfl-number in the
	// cell with the maximum cfl-number being just cfl_value. The only
	// exceptions are outputs at specific times 

	gdata.dt = cfl_set/gdata.cfl;
	gdata.dt /= factor;

	// Decrease timestep if movie-output comes too soon
#if (CRONOS_MOVIE == CRONOS_ON)
	double tnext_mov = dt_mov*(nummov_done);
	double dtnext_mov = tnext_mov - gdata.time;

	if(gdata.dt > dtnext_mov) {
		gdata.dt = dtnext_mov;
	}
#endif

	// Decrease timestep if float-output comes too soon
	double tnext_flt = dt_flt*(numflt_done);
	double dtnext_flt = tnext_flt - gdata.time;

	if(gdata.dt > dtnext_flt && numflt_pass < numflt_out) {
		gdata.dt = dtnext_flt;
	}

	// Decrease timestep if ascii-output comes too soon
	double tnext_ascii = dt_ascii*(numascii_done+1);
	double dtnext_ascii = tnext_ascii - gdata.time;

	if(gdata.dt > dtnext_ascii) {
		gdata.dt = dtnext_ascii;
	}

#endif

}

//! Compares timestep-size against the output step-sizes and potentially performs writeouts
void Environment::CheckOut(Data &gdata, Queue &queue)
{
	int outputflag[5] = {0, 0, 0, 0, 0};
	bool buffer_fetch_necessary = false;

	numdbl_pass++;
	numflt_pass++;
	numascii_pass++;
	numinfo_pass++;

#ifdef DTCOMP_OLD

	double dt_eps = 1.e-4*gdata.dt;
  
	if(gdata.rank==0){


		/* If:
		   (a) More than "num_out" steps were done since last output
		   or
		   (b) More than "dt_out" time was integrated since last output
		   a HDF-output is triggered
		   But: Careful - condition (b) does not work properly if dt
		   changed in between.
		*/

		double tmod = abs(fmod(gdata.time+dt_eps,dt_dbl));
		double trun = gdata.time - restart_time;
		if ((trun <= 1e-6) || (dt_dbl <= 0) ||
		    (gdata.time - t_last_dbl + dt_eps) < dt_dbl) tmod =  1.;

		if(((dt_dbl > 0) && (tmod <= 10*dt_eps)) || numdbl_pass >= numdbl_out){
			outputflag[0] = 1;
			t_last_dbl = gdata.time;
			buffer_fetch_necessary = true;
		}


		tmod = abs(fmod(gdata.time+dt_eps,dt_flt));
		if ((trun <= 1e-6) || (dt_flt <= 0) ||
		    (gdata.time - t_last_flt + dt_eps) < dt_flt) tmod =  1.;

		if(((dt_flt > 0) && (tmod <= 10*dt_eps)) || numflt_done >= numflt_out){
			outputflag[1] = 1;
			t_last_flt = gdata.time;
			buffer_fetch_necessary = true;
		}


		tmod = abs(fmod(gdata.time+dt_eps,dt_ascii));
		if ((trun <= 1e-6) || (dt_ascii <= 0)) tmod =  1.;

		if(((dt_ascii > 0) && (tmod <= 10*dt_eps)) || numascii_done >= numascii_out){
			outputflag[2] = 1;
			buffer_fetch_necessary = true;
		}


		tmod = abs(fmod(gdata.time+dt_eps,dt_info));
		if ((trun <= 1e-6) || (dt_info <= 0)) tmod =  1.;

		if(((dt_info > 0) && (tmod <= 10*dt_eps)) || numinfo_done >= numinfo_out){
			outputflag[3] = 1;
			buffer_fetch_necessary = true;
		}


#if (CRONOS_MOVIE == CRONOS_ON)
		tmod = abs(fmod(gdata.time+dt_eps,dt_mov));
		if ((trun <= 1e-6) || (dt_mov <= 0)) tmod =  1.;

		if((dt_mov > 0) && (tmod <= 10*dt_eps)) {
			outputflag[4] = 1;
			buffer_fetch_necessary = true;
		}
	}
#endif
#else

	// if(gdata.rank == 0) {
		const double dt_eps = gdata.dt*1.e-2;

		// Double output
		const double tnext_dbl = dt_dbl*(numdbl_done);
		const double dtnext_dbl = tnext_dbl - gdata.time;

		if(dtnext_dbl < dt_eps || numdbl_pass >= numdbl_out) {
			outputflag[0] = 1;
			buffer_fetch_necessary = true;
		}

    
		// Float output
		const double tnext_flt = dt_flt*(numflt_done);
		const double dtnext_flt = tnext_flt - gdata.time;

		if(dtnext_flt < dt_eps || numflt_pass >= numflt_out) {
			outputflag[1] = 1;
			buffer_fetch_necessary = true;
		}


		// Ascii output
		const double tnext_ascii = dt_ascii*(numascii_done);
		const double dtnext_ascii = tnext_ascii - gdata.time;

		if(dtnext_ascii < dt_eps || numascii_pass >= numascii_out) {
			outputflag[2] = 1;
			buffer_fetch_necessary = true;
		}


		// Screen output
		const double tnext_info = dt_info*(numinfo_done);
		const double dtnext_info = tnext_info - gdata.time;

		if(dtnext_info< dt_eps || numinfo_pass >= numinfo_out) {
			outputflag[3] = 1;
			buffer_fetch_necessary = true;
		}


		// Movie output
#if (CRONOS_MOVIE == CRONOS_ON)
		double tnext_mov = dt_mov*(nummov_done);
		double dtnext_mov = tnext_mov - gdata.time;

		if(dtnext_mov < dt_eps) {
			outputflag[4] = 1;
			buffer_fetch_necessary = true;
		}
#endif


	// }
  
#endif

#if (CRONOS_DISTR_OUTPUT != CRONOS_ON)
	if (buffer_fetch_necessary) {
		FetchDataBuffer(gdata, queue);
	}
#endif

	// doing double output
	if(outputflag[0] == 1) {
#if (CRONOS_DISTR_OUTPUT == CRONOS_ON)
		Output_Distributed(queue, gdata, false);
#else
		if(outputflag[1] == 1) {
			Output_Master(queue, gdata, false, true);
		} else {
			Output_Master(queue, gdata, false, false);
		}
#endif
		
		// Output(gdata, false, false);
		numdbl_pass = 0;
	}

	// doing float output
	if(outputflag[1] == 1) {
#if (CRONOS_DISTR_OUTPUT == CRONOS_ON)
		Output_Distributed(queue, gdata, true);
#else
		Output_Master(queue, gdata, true, false);
#endif
		// Output(gdata, true, false);
		numflt_pass = 0;
	}

	// setting ascii output
	if(outputflag[2] == 1) {
		Problem->set_AsciiOut(true);
		numascii_done++;
		numascii_pass = 0;
	} else {
		Problem->set_AsciiOut(false);
	}

	// setting info output
	if(outputflag[3] == 1) {
		Problem->set_Info(true);
		numinfo_done++;
		numinfo_pass = 1;
	} else {
		Problem->set_Info(false);
	}


#if (CRONOS_MOVIE == CRONOS_ON)
	if(outputflag[4] == 1) {
		WriteMovies(gdata);
	}
#endif
}

//! @brief Performs a fetch of the Celerity buffer into the NumMatrix data structure and synchronizes
void Environment::FetchDataBuffer(Data &gdata, Queue &queue)
{
	auto omRange = gdata.omSYCL[0].get_range();
	size_t nom_max[3] = {omRange.get(0), omRange.get(1), omRange.get(2)};
	for (int q = 0; q < N_OMINT; q++) {

		queue.submit(celerity::allow_by_ref, [=, &gdata](celerity::handler& cgh) {
			celerity::accessor omSYCL_acc{gdata.omSYCL[q], cgh, celerity::access::all{}, celerity::read_only_host_task};
			cgh.host_task(celerity::on_master_node, [=, &gdata]{
				for (int i = 0; i < nom_max[0]; i++) {
					for (int j = 0; j < nom_max[1]; j++) {
						for (int k = 0; k < nom_max[2]; k++) {
							gdata.om[q](i-3,j-3,k-3) = omSYCL_acc[i][j][k];
						}
					}
				}
			});
		});
	}
	queue.slow_full_sync();
}


#if (CRONOS_MOVIE == CRONOS_ON)
//! @brief Performs the problems specified movie output
void Environment::WriteMovies(Data &gdata)
{
	if(gdata.rank == 0 && Problem->checkout(2)) {
		cout << "......................................................" << endl;
		cout << "    Movie Output triggered " << endl;
		cout << "......................................................" << endl;
	}
  
	Problem->writeMovie(gdata, mov);
	nummov_done++;
}
#endif


//! @brief Performs time step and checks for end conditions
int Environment::integrate(Data &gdata, Queue& queue)
{
	int EndProgram(0);
	outputflag = 0;
  
	if(gdata.rank==0){
		//    cout << "---------------------------------------------" << endl;
		cout << "======================================================" << endl;
		cout << "Integrating: time = " << gdata.time;
#ifndef DTCOMP_OLD
		cout << " with dt = " << gdata.dt;
#endif
		cout << " at step " << gdata.tstep;
		cout << endl;
	}

	try{
		pdestep(gdata, queue);
	} catch (CException exep) {
		Abort(gdata, exep);
	}

	gdata.tstep++;

	CheckOut(gdata, queue);

	try {
		compute_dt(gdata);
	} catch (CException exep) {
		Abort(gdata, exep);
	}

	// Check if possible user-defined ending condition has been reached
	EndProgram = CheckEnd_User(gdata, queue);

	// Check if ending condition is reached:
	if(EndProgram==0) {
		EndProgram = CheckEnd(gdata, queue);
	}

	// step++;
	// gdata.tstep++;

	//if(gdata.tstep==1) {
	//	EndProgram = Finalize(gdata, queue, "stopping for valgrind");
	//}

	return EndProgram;
}


//! @brief Performs final output and triggers the application's end condition
int Environment::Finalize(Data &gdata, Queue &queue, string message)
{
	if(gdata.rank == 0) {
		cout << " Program exiting normally: " << endl;
		cout << " ------ " << message << " ------ " << endl;
	}


#if (CRONOS_DISTR_OUTPUT == CRONOS_ON)
		Output_Distributed(queue, gdata, true);
#else	
		FetchDataBuffer(gdata, queue);
		Output_Master(queue, gdata, true, false);
#endif

	return 1; // Indicate normal end of program

}

//! @brief Aborts the application and outputs error
void Environment::Abort(Data &gdata, CException exep)
{
	cout << flush;
	cerr << " Program aborted due to this error: " << endl;
	cerr << exep.returnReport() << endl;

	int errorcode = -1;
	exit(errorcode);
}

//! @brief Performs the application's HDF5 writeout of its NumMatrix structure on the master node
void Environment::Output_Master(Queue &queue, Data &gdata, bool isfloat, bool sync)
{
	
	queue.submit(celerity::allow_by_ref, [=, &gdata](celerity::handler& cgh) {
		cgh.host_task(celerity::on_master_node, [=, &gdata] {
			Output(gdata, isfloat, false, false);
		});
	});

	//now adds one pre full completion. Probably not a problem
	if (isfloat) {
		numflt_done++;
	} else {
		numdbl_done++;
	}

	if (sync) {
		queue.slow_full_sync();
	}

}

//! @brief Performs a distributed version of the application's HDF5 writeout of its NumMatrix structure
void Environment::Output_Distributed(Queue &queue, Data &gdata, bool isfloat)
{
	/*
	  Routine to write a hdf5 output file for double data -- can be used
	  for restart
	*/

	if(gdata.rank == 0 && Problem->checkout(2)) {
		cout << "......................................................" << endl;
		if(isfloat) {
			cout << "   Output (flt) [distr_io] ";
		} else {
			cout << "   Output (dbl) [distr_io] ";
		}
		cout << " triggered at time: " << gdata.time << endl; 
	}

	string filename;
	char dirname[255];

	if(isfloat) {
		sprintf(dirname, "%s/%s_float", getenv("poub"),getenv("pname"));
	} else {
		sprintf(dirname, "%s/%s_double", getenv("poub"),getenv("pname"));
	}
	for (char *tmp=dirname+strlen(dirname)-1; *tmp=='0'; tmp--) *tmp=0;

	queue.submit([=](celerity::handler& cgh) {
		cgh.host_task(celerity::on_master_node, [=] {
			filesystem::create_directory(dirname);
			filesystem::permissions(dirname, filesystem::perms::owner_read | filesystem::perms::owner_exec | filesystem::perms::owner_write | filesystem::perms::group_read |
					 filesystem::perms::group_exec | filesystem::perms::group_write | filesystem::perms::others_exec);
			auto ftime = filesystem::last_write_time(dirname);
			filesystem::file_time_type currentTime(decltype(ftime)::clock::now());
			filesystem::last_write_time(dirname, currentTime);
		});
	});

	filename = dirname;
	filename += "/";
	filename += MakeFilename(gdata.coords, gdata.tstep, false, isfloat);

	int n_Omega = N_OMEGA;
	int n_omIntUser = N_OMINT_USER;
	int n_omIntAll = N_OMINT_ALL;

	int numout(n_omIntAll);

	if(!isfloat) {
		numout = n_omIntAll + N_SUBS;
	}

	Hdf5Stream *h5out;
	h5out = new Hdf5Stream(filename, numout, gdata.rank, true);

	if (!isfloat) {
		for (int q = 0; q < gdata.omSYCL.size(); q++) {
			queue.submit(celerity::allow_by_ref, [=, &gdata](celerity::handler& cgh) {
				celerity::accessor omSYCL_out_acc{gdata.omSYCL_out[q], cgh, celerity::access::one_to_one{}, celerity::write_only, celerity::no_init};
				celerity::accessor omSYCL_acc{gdata.omSYCL[q], cgh, celerity::access::neighborhood{6,6,6}, celerity::read_only};
				cgh.parallel_for<class BufferInitializationKernel>(gdata.omSYCL_out[q].get_range(), [=](celerity::item<3> item) {
					omSYCL_out_acc[item.get_id(0)][item.get_id(1)][item.get_id(2)] = omSYCL_acc[item.get_id(0) + 3][item.get_id(1) + 3][item.get_id(2) + 3];
				});
			});
		}

		queue.submit(celerity::allow_by_ref, [=, &gdata](celerity::handler& cgh) {
			celerity::accessor om_rho_acc{gdata.omSYCL_out[0], cgh, celerity::experimental::access::even_split<3>{}, celerity::read_only_host_task};
			celerity::accessor om_sx_acc{gdata.omSYCL_out[1], cgh, celerity::experimental::access::even_split<3>{}, celerity::read_only_host_task};
			celerity::accessor om_sy_acc{gdata.omSYCL_out[2], cgh, celerity::experimental::access::even_split<3>{}, celerity::read_only_host_task};
			celerity::accessor om_sz_acc{gdata.omSYCL_out[3], cgh, celerity::experimental::access::even_split<3>{}, celerity::read_only_host_task};
			celerity::accessor om_Eges_acc{gdata.omSYCL_out[4], cgh, celerity::experimental::access::even_split<3>{}, celerity::read_only_host_task};
			cgh.host_task(celerity::experimental::collective, [=, &gdata](celerity::experimental::collective_partition part) {

				auto comm = part.get_collective_mpi_comm();
				h5out->Initialize(comm);

				auto om_rho = om_rho_acc.get_allocation_window(part);
				auto om_sx = om_sx_acc.get_allocation_window(part);
				auto om_sy = om_sy_acc.get_allocation_window(part);
				auto om_sz = om_sz_acc.get_allocation_window(part);
				auto om_Eges = om_Eges_acc.get_allocation_window(part);

				distr_io::dataout(gdata, *h5out, *Problem, numout, isfloat, om_rho,
											om_sx, om_sy, om_sz, om_Eges);
				
				rksolver->addConstants_toH5(gdata, *h5out);

				delete h5out;
			});
		});

	} else {
		for (int q = 0; q < gdata.omSYCL.size(); q++) {
			queue.submit(celerity::allow_by_ref, [=, &gdata](celerity::handler& cgh) {
				celerity::accessor omSYCL_out_flt_acc{gdata.omSYCL_out_flt[q], cgh, celerity::access::one_to_one{}, celerity::write_only, celerity::no_init};
				celerity::accessor omSYCL_acc{gdata.omSYCL[q], cgh, celerity::access::neighborhood{6,6,6}, celerity::read_only};
				cgh.parallel_for<class BufferInitializationKernel>(gdata.omSYCL_out_flt[q].get_range(), [=](celerity::item<3> item) {
					omSYCL_out_flt_acc[item.get_id(0)][item.get_id(1)][item.get_id(2)] = (q < 0) ? 0. : static_cast<float>(omSYCL_acc[item.get_id(0) + 3][item.get_id(1) + 3][item.get_id(2) + 3]);
				});
			});
		}

		queue.submit(celerity::allow_by_ref, [=, &gdata](celerity::handler& cgh) {
			celerity::accessor om_rho_acc{gdata.omSYCL_out_flt[0], cgh, celerity::experimental::access::even_split<3>{}, celerity::read_only_host_task};
			celerity::accessor om_sx_acc{gdata.omSYCL_out_flt[1], cgh, celerity::experimental::access::even_split<3>{}, celerity::read_only_host_task};
			celerity::accessor om_sy_acc{gdata.omSYCL_out_flt[2], cgh, celerity::experimental::access::even_split<3>{}, celerity::read_only_host_task};
			celerity::accessor om_sz_acc{gdata.omSYCL_out_flt[3], cgh, celerity::experimental::access::even_split<3>{}, celerity::read_only_host_task};
			celerity::accessor om_Eges_acc{gdata.omSYCL_out_flt[4], cgh, celerity::experimental::access::even_split<3>{}, celerity::read_only_host_task};
			cgh.host_task(celerity::experimental::collective, [=, &gdata](celerity::experimental::collective_partition part) {

				auto comm = part.get_collective_mpi_comm();
				h5out->Initialize(comm);

				auto om_rho = om_rho_acc.get_allocation_window(part);
				auto om_sx = om_sx_acc.get_allocation_window(part);
				auto om_sy = om_sy_acc.get_allocation_window(part);
				auto om_sz = om_sz_acc.get_allocation_window(part);
				auto om_Eges = om_Eges_acc.get_allocation_window(part);

				distr_io::dataout(gdata, *h5out, *Problem, numout, isfloat, om_rho,
											om_sx, om_sy, om_sz, om_Eges);
				
				rksolver->addConstants_toH5(gdata, *h5out);

				delete h5out;
			});
		});
	}

	if(!isfloat) {
		char OldChar[255] = {};
		OldName.copy(OldChar,255);
		OldName.clear();
		//    cout << " Removing: " << OldChar << endl;
		remove(OldChar);
		OldName = filename;
	}
	if(gdata.rank == 0 && Problem->checkout(2)){
		cout << "......................................................" << endl;
	}

	if(isfloat) {
		numflt_done++;
	} else {
		numdbl_done++;
	}


}


//! @brief Routine to write a hdf5 output file
void Environment::Output(Data &gdata, bool isfloat, bool terminate, bool progress_counters)
{
	/*
	  Routine to write a hdf5 output file for double data -- can be used
	  for restart
	*/
	if(gdata.rank == 0 && Problem->checkout(2)) {
		cout << "......................................................" << endl;
		if(isfloat) {
			cout << "   Output (flt) ";
		} else {
			cout << "   Output (dbl) ";
		}
		cout << " triggered at time: " << gdata.time << endl; 
	}

	string filename;
	char dirname[255];

	if(isfloat) {
		sprintf(dirname, "%s/%s_float", getenv("poub"),getenv("pname"));
	} else {
		sprintf(dirname, "%s/%s_double", getenv("poub"),getenv("pname"));
	}
	for (char *tmp=dirname+strlen(dirname)-1; *tmp=='0'; tmp--) *tmp=0;

	//mkdir(dirname, 511);
	filesystem::create_directory(dirname);
	filesystem::permissions(dirname, filesystem::perms::owner_read | filesystem::perms::owner_exec | filesystem::perms::owner_write | filesystem::perms::group_exec | filesystem::perms::others_exec);
	//utime(dirname, NULL);
	auto ftime = filesystem::last_write_time(dirname);
	filesystem::file_time_type currentTime(decltype(ftime)::clock::now());
	filesystem::last_write_time(dirname, currentTime);

	filename = dirname;
	filename += "/";
	filename += MakeFilename(gdata.coords, gdata.tstep, terminate, isfloat);

#if(FLUID_TYPE == CRONOS_MULTIFLUID)
	int n_Omega = gdata.fluids->get_N_OMEGA();
	int n_omIntUser = gdata.fluids->get_N_OMINT_USER();
	int n_omIntAll = gdata.fluids->get_N_OMINT_ALL();
#else
	int n_Omega = N_OMEGA;
	int n_omIntUser = N_OMINT_USER;
	int n_omIntAll = N_OMINT_ALL;
#endif
  
	int numout(n_omIntAll);
	if(terminate) {
		numout = n_Omega+ n_omIntAll -N_ADD;
	}

	if(!isfloat) {
		numout = n_omIntAll + N_SUBS;
	}
//		numout = n_omIntAll + N_SUBS;

	//	Hdf5Stream h5out(filename, numout, gdata.rank);

	Hdf5Stream *h5out;

#if(HDF_PARALLEL_IO == CRONOS_ON)
	h5out = new Hdf5Stream(filename, numout, gdata.rank, true);
#else // New Output version
	h5out = new Hdf5Stream(filename, numout, gdata.rank, false);
#endif

	// // Apply possible user-change to fields
	// if(isfloat) {
	// 	Problem->TransConsCharUser(gdata);
	// }

	//	cout << " dataout at " << gdata.tstep << endl;
#if(HDF_PARALLEL_IO == CRONOS_ON)
#if(CRONOS_OUTPUT_COMPATIBILITY == CRONOS_ON)
	gfunc->dataout(gdata, *h5out, *Problem, numout, isfloat, isfloat);
#else
	gfunc->dataout(gdata, *h5out, *Problem, numout, isfloat, true);
#endif
#else
	gfunc->dataout(gdata, *h5out, *Problem, numout, isfloat, false);
#endif

	// // Apply possible user-change to fields
	// if(isfloat) {
	// 	Problem->TransCharConsUser(gdata);
	// }

	// Now add the data from the rksolver:
	rksolver->addConstants_toH5(gdata, *h5out);

	// Now add the data from the eulersolver:

	/*for (auto &EulerSolver : EulerSolvers) {
		EulerSolver->addConstants_toH5(gdata, *h5out);
	}*/

	delete h5out;

	// if(!isfloat) {
	// 	char OldChar[255] = {};
	// 	OldName.copy(OldChar,255);
	// 	OldName.clear();
	// 	//    cout << " Removing: " << OldChar << endl;
	// 	remove(OldChar);
	// 	OldName = filename;
	// }
	if(gdata.rank == 0 && Problem->checkout(2)){
		cout << "......................................................" << endl;
	}

	if (progress_counters) {
		if(isfloat) {
			numflt_done++;
		} else {
			numdbl_done++;
		}
	}

}

//! @brief Routine to load data from a previous timestep
void Environment::LoadData(Data &gdata)
{

	/*
	  Routine to load data from a previous timestep.

	*/
	this->restart_step = static_cast<int>(value((char*)"restart_step"));
	gdata.tstep = this->restart_step;

	char dirname[255];
	sprintf(dirname, "%s/%s_double", getenv("poub"),getenv("pname"));
	for (char *tmp=dirname+strlen(dirname)-1; *tmp=='0'; tmp--) *tmp=0;

	double factor = 1.;

#if((CRONOS_OUTPUT_COMPATIBILITY == CRONOS_ON) || (HDF_PARALLEL_IO == CRONOS_OFF) )
	int coordsOld[DIM];
	int numFiles[DIM] = {1,1,1};
	// Checking if number of processes is changed:
	try {
		getCoords(gdata, static_cast<string>(dirname), coordsOld, numFiles);
	} catch (CException exep) {
		Abort(gdata, exep);
	}

	// Iterating over necessary number of files


	int coordsLoad[DIM];
	for(int i=0; i<numFiles[0]; ++i) {
		coordsLoad[0] = coordsOld[0] + i;
		for(int j=0; j<numFiles[1]; ++j) {
			coordsLoad[1] = coordsOld[1] + j;
			for(int k=0; k<numFiles[2]; ++k) {
				coordsLoad[2] = coordsOld[2] + k;
	

				string filename = dirname;
				filename += "/";
				filename += MakeFilename(coordsLoad, restart_step, false, false);
	
				if(gdata.rank == 0) {
					if (!file_exists(filename)) {
						cerr << "Could not find restart file" << endl;
						exit(1);
					}
				}
				try {
					Hdf5iStream h5in(filename);
					gfunc->datain(gdata, h5in, filename, *Problem);

					// Get additional constants for Runge Kutta solver
					rksolver->getConstants_fromH5(gdata, h5in);

				} catch (CException exep) {
					Abort(gdata, exep);
				}
			}
		}
	}
#else
	string filename = dirname;
	filename += "/";
	filename += MakeFilename(gdata.coords, restart_step, false, false);
	cout << " Name: " << filename << endl;

	if(gdata.rank == 0) {
		if (!file_exists(filename)) {
			cerr << "Could not find restart file" << endl;
			exit(1);
		}
	}

	Hdf5iStream h5in(filename);
	factor = gfunc->datain_collective(gdata, h5in, filename, *Problem);

	// Get additional constants for Runge Kutta solver
	rksolver->getConstants_fromH5(gdata, h5in);
#endif

	// Recompute times for data output:
	double dt_eps = gdata.dt*1.e-2;
	numdbl_done = static_cast<int>((gdata.time+dt_eps)/dt_dbl) + 1;
	numflt_done = static_cast<int>((gdata.time+dt_eps)/dt_flt) + 1;
	numascii_done = static_cast<int>((gdata.time+dt_eps)/dt_ascii) + 1;
	numinfo_done = static_cast<int>((gdata.time+dt_eps)/dt_info) +1;
	numinfo_pass = 2;
#if (CRONOS_MOVIE == CRONOS_ON)
	nummov_done = static_cast<int>((gdata.time+dt_eps)/dt_mov) + 1;
#endif

	compute_dt(gdata, factor);
}

#if(FLUID_TYPE != CRONOS_MULTIFLUID)
//! @brief Routine to load floating point data from a previous timestep
void Environment::LoadData_flt(Data &gdata, int load_step)
{

	/*
	  Routine to load data from a previous timestep.
	*/

	gdata.tstep = load_step;

	char dirname[255];
	sprintf(dirname, "%s/%s_float", getenv("poub"),getenv("pname"));
	for (char *tmp=dirname+strlen(dirname)-1; *tmp=='0'; tmp--) *tmp=0;


	string filename = dirname;
	filename += "/";
	filename += MakeFilename(NULL, load_step, false, true);

	if(gdata.rank == 0) {
		if (!file_exists(filename)) {
			cerr << "Could not find restart file" << endl;
			exit(1);
		}
	}
	try {
		Hdf5iStream h5in(filename, gdata.rank);
		gfunc->load_flt(gdata, h5in, filename, *Problem);

	} catch (CException exep) {
		Abort(gdata, exep);
	}


	// Recompute times for data output:
	double dt_eps = gdata.dt*1.e-2;
	numdbl_done = static_cast<int>((gdata.time+dt_eps)/dt_dbl) + 1;
	numflt_done = static_cast<int>((gdata.time+dt_eps)/dt_flt) + 1;
	numascii_done = static_cast<int>((gdata.time+dt_eps)/dt_ascii) + 1;
	numinfo_done = static_cast<int>((gdata.time+dt_eps)/dt_info) +1;
	numinfo_pass = 2;
#if (CRONOS_MOVIE == CRONOS_ON)
	nummov_done = static_cast<int>((gdata.time+dt_eps)/dt_mov) + 1;
#endif

	compute_dt(gdata);
}
#endif


string Environment::MakeFilename(const int coords[DIM], int tstep,
                                 bool terminate, bool isfloat)
{


	std::stringstream sstr_fname;
	sstr_fname << getenv("pname");

	if(isfloat) {
		sstr_fname << "_flt_step";
	} else {
		sstr_fname << "_dbl_step";
	}

	if(terminate) {
		sstr_fname << "error";
	} else {
#ifdef FNAME_DIGITS
		sstr_fname << std::setfill('0') << std::setw(FNAME_DIGITS) << tstep;
#else
		sstr_fname << tstep;
#endif
	}

	string filename = sstr_fname.str();
  
	filename += ".h5";

	return filename;
}

void Environment::getCoords(const Data &gdata, string dirname,
                            int coords[DIM], int numFiles[DIM])
{

	int facloc[DIM];
	for(int i=0; i<DIM; ++i) {
		facloc[i] = 1;
		numFiles[i] = 1;
	}

	// Giving the equivalent coordinates of the old jobs.
	for(int i=0; i<DIM; ++i) {
		if(facloc[i] > 1) {
			coords[i] = gdata.coords[i]/facloc[i];
		} else if (facloc[i] < 1) {
			coords[i] = -gdata.coords[i]*facloc[i];
		} else {
			coords[i] = gdata.coords[i];
		}
	}

}
