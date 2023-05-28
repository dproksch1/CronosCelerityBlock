#include <stdio.h>
#include <stdlib.h>
#include <fstream>
#include <matrix.H>

#define ISODD(i) ((i)&1)
////#define REAL double
//#define NONE -1

#include <sys/types.h>
#include <sys/stat.h>   
#ifdef _MSC_VER
#include "timewrapper.H"
#else
#include <utime.h>
#endif
//#include "celerity_test.H"

#include <fenv.h>
#include <signal.h>
#include "CException.H"
#include "data.H"
#include "gridgen.H"
#include "queue.H"

cronos_ostream ccout(std::cout, 0);
cronos_ostream ccerr(std::cerr, 0);

int main(int argc, char* argv[])
{
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

	//DeviceSelector device_selector;
	//Queue queue(device_selector);
	Queue queue;

	Data gdata;
	Environment solver(gdata);

	//cout << gdata.fluid.get_fieldName(qReconst) << "\n";

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
			solver.RKSolver->init(gdata, *solver.gfunc, *solver.Problem, queue);
			// solver->EulerSolver->init(gdata, *solver->gfunc, *solver->Problem);
		} catch (CException exep) {
			solver.Abort(gdata, exep);
		}
		solver.InitOutput(gdata);
	}

	// for (;;) {
	while (EndProgram == 0) {
		EndProgram = solver.integrate(gdata, queue);
	//	// EndProgram = solver.Finalize(gdata, queue, static_cast<string>("Stopping for gdm"));

	}

// queue.slow_full_sync();
	if (EndProgram == 1) {
		// Everything is fine -> exit normally
		return 0;
	} else {
		return EndProgram;
	}
}










