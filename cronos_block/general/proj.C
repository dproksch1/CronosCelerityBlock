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
		setenvWrapper("pname", argv[1]);
	} else if (argc == 3) {
		setenvWrapper("poub", argv[1]);
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

	Queue queue;

	Data gdata;
	Environment solver(gdata);

	bool RESTART = bool(value((char*)"restart"));

	if (RESTART) {
		solver.LoadData(gdata);
		try {
			solver.rksolver->restart(queue, gdata, *solver.gfunc, *solver.Problem);
		} catch (CException exep) {
			solver.Abort(gdata, exep);
		}
	} else {
		gdata.time = 0.;
		try {
			solver.rksolver->init(queue, gdata, *solver.gfunc, *solver.Problem);
			// solver->EulerSolver->init(gdata, *solver->gfunc, *solver->Problem);
		} catch (CException exep) {
			solver.Abort(gdata, exep);
		}
		solver.InitOutput(queue, gdata);
	}

	while (EndProgram == 0) {
		EndProgram = solver.integrate(gdata, queue);
	}

	// manual execution of the queue's destructor to prevent desynchronized process termination
	queue.~distr_queue();

	if (EndProgram == 1) {
		// Everything is fine -> exit normally
		return 0;
	} else {
		return EndProgram;
	}
}
