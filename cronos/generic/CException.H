#ifndef CRONOS_CEXCEPTION_H
#define CRONOS_CEXCEPTION_H = 1

#include <string>
using namespace std;

class CException
{
public:
	CException(string);
	CException(string, int);
	CException(string, int, int, int, int);
	CException(string, int, int, int, int, int);
	string returnError();
	string returnLocation();
	string returnThrower();

	string returnReport();

	const char* what() const noexcept;
	int type;
private:
	bool stacktrace = true;

	void finish();
	string ErrorMessage;
	string Thrower;
	string whatMessage;
	string report;
	int loc[3];
	int q;
};

#ifdef __linux__
#include <execinfo.h>
#include <sstream>
#include <cxxabi.h>
#endif

/** Returns the demangled stack backtrace of the caller function */
std::string get_stacktrace(unsigned int max_frames = 64);


#include <iostream>
#include <string.h>

/** Signal handler that prints a stacktrace if called */
void stacktraceHandler(int sig);

#endif
