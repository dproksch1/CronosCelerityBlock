#include <stdlib.h>
#include <fstream>
#include <matrix.H>

#define ISODD(i) ((i)&1)     

#include <sys/types.h>
#include <sys/stat.h>   
#ifdef _MSC_VER
#include "timewrapper.H"
#else
#include <utime.h>
#endif

#include <fenv.h>
#include <signal.h>

#include "gridgen.H"
#ifdef parallel
#include "mpi.h"
#endif

#include "queue.H"

cronos_ostream ccout(std::cout, 0);
cronos_ostream ccerr(std::cerr, 0);

int main(int argc, char* argv[])
{

#ifdef parallel
	int ntasks;
	MPI_Init(&argc, &argv);
	MPI_Comm_size(MPI_COMM_WORLD, &ntasks);
	cout << " Number of tasks: " << ntasks << endl;

	ccout.checkRank();
	ccerr.checkRank();
#endif
	/*
	  Additional arguments:
	  (a) Programm name
	  (b) data directory, Program name
	*/

	auto setenvWrapper = [](const char* varName, char* value) {
#ifdef _MSC_VER
		_putenv_s(varName, value);
#else
		setenv(varName, value, 1);
#endif
	};

	if (argc == 2) {
		//setenv("pname", argv[1], 1);
		setenvWrapper("pname", argv[1]);
	} else if (argc == 3) {
		//setenv("poub",argv[1],1);
		setenvWrapper("poub", argv[1]);
		//setenv("pname",argv[2],1);
		setenvWrapper("pname", argv[2]);
	}

	feclearexcept(FE_ALL_EXCEPT);  		// Clear all exceptions
	signal(SIGFPE, SIG_IGN);       		// Ignore Floating point errors (zero division)
	signal(SIGSEGV, stacktraceHandler); // install the stacktracer for SEGFAULTS
	signal(SIGABRT, stacktraceHandler); // install the stacktracer for aborts (e.g. by failed asserts)

#ifdef RK_STEPS
	// N_P is given in constants.H - number of poisson-steps
	if (N_P) {
		cerr << "N_P should be 0 for Runge Kutta" << endl;
		exit(2);
	}
#endif

	int EndProgram(0);

	{

		DeviceSelector device_selector;
		Queue queue(device_selector);

	Data gdata;

	// Environment *solver;
	// solver = new Environment(gdata);
	Environment solver(gdata);
	//  solver->setType(gdata);


	bool RESTART = bool(value((char*)"restart"));

	if (RESTART) {
		solver.LoadData(gdata);
		try {
			solver.RKSolver->restart(gdata, *solver.gfunc, *solver.Problem);
		} catch (CException exep) {
			solver.Abort(gdata, exep);
		}
	} else {
		gdata.time = 0.;
		try {
			solver.RKSolver->init(gdata, *solver.gfunc, *solver.Problem);
			// solver->EulerSolver->init(gdata, *solver->gfunc, *solver->Problem);
		} catch (CException exep) {
			solver.Abort(gdata, exep);
		}
		solver.InitOutput(gdata);
	}

		std::cout << "Running on " << queue.get_device().get_info<cl::sycl::info::device::name>() << std::endl;

		// for (;;) {
		while (EndProgram == 0) {
			EndProgram = solver.integrate(gdata, queue);
		//	// EndProgram = solver.Finalize(gdata, static_cast<string>("Stopping for gdm"));

		}


#ifdef parallel
	MPI_Barrier(gdata.comm3d);
	MPI_Finalize();
#endif
// #ifdef parallel
// 	MPI_Finalize();
// #endif
	//queue.wait();
	std::cout << "finished finalizing, still inside SYCL scope" << std::endl << std::flush;
	}
	std::cout << "finished finalizing, outside SYCL scope, exit code " << EndProgram << std::endl << std::flush;

	if (EndProgram == 1) {
		// Everything is fine -> exit normally
		return 0;
	} else {
		return EndProgram;
	}
}










