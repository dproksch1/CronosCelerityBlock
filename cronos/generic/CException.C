#include "CException.H"
#include <sstream>
#include <stdio.h>

CException::CException(string ErrorMessage)
{
	this->ErrorMessage = ErrorMessage;
	this->type = 0;
	finish();
};

CException::CException(string ErrorMessage, int pos)
{
	std::ostringstream dummy;
	char cpos[255];
	sprintf(cpos,"%2.2d",pos);
	dummy << ErrorMessage << " at checkpoint: " << cpos;
	this->ErrorMessage = dummy.str();
	this->type = 0;
	finish();
};

CException::CException(string ErrorMessage, int q, int i, int j, int k)
{
	this->ErrorMessage = ErrorMessage;
	this->q      = q;
	this->loc[0] = i;
	this->loc[1] = j;
	this->loc[2] = k;
	this->type = 1;
	finish();
};

CException::CException(string ErrorMessage, int pos, int q, int i, int j, int k)
{
	std::ostringstream dummy;
	char cpos[255];
	sprintf(cpos,"%2.2d",pos);
	dummy << ErrorMessage << " at checkpoint: " << cpos;
	this->ErrorMessage = dummy.str();
	this->q      = q;
	this->loc[0] = i;
	this->loc[1] = j;
	this->loc[2] = k;
	this->type = 1;
	finish();
};

void CException::finish()
{
	whatMessage = ErrorMessage;
	if (Thrower != "") {
		whatMessage += " -> " + Thrower;
	}

	if (stacktrace) {
		whatMessage += get_stacktrace();
	}

	stringstream buffer;
	buffer << " ------ " << returnError() << " ------ " << endl;
	if (type == 1){
		buffer << "  which occurred for om" << returnLocation() << endl;
	}
	if (Thrower != "") {
		buffer << "  -> " << Thrower << endl;
	}
	if (stacktrace) {
		buffer << get_stacktrace();
	}
	report = buffer.str();
}


string CException::returnError()
{
	return this->ErrorMessage;
};

string CException::returnLocation()
{
	string Location;
	char charloc[255];
	Location+="[";
	sprintf(charloc,"%i",q);
	Location+=charloc;
	Location += "](";
	sprintf(charloc,"%i",loc[0]);
	Location += charloc;
	Location += ",";
	sprintf(charloc,"%i",loc[1]);
	Location += charloc;
	Location += ",";
	sprintf(charloc,"%i",loc[2]);
	Location += charloc;
	Location += ")";
	return Location;
};

string CException::returnReport() {
//	finish();
//	cout << " My report " << report << endl;
	return report;
}


string CException::returnThrower()
{
	return Thrower;
};

const char* CException::what() const noexcept
{
	return whatMessage.c_str();
}

#ifdef __linux__
std::string get_stacktrace(unsigned int max_frames)
{
	std::stringstream buffer;
	buffer << "stack trace:" << std::endl;

	// storage array for stack trace address data
	void* addrlist[max_frames];

	// retrieve current stack addresses
	int addrlen = backtrace(addrlist, max_frames);

	if (addrlen == 0) {
		buffer << "  <empty, possibly corrupt>" << std::endl;
		return buffer.str();
	}

	// resolve addresses into strings containing "filename(function+address)",
	// this array must be free()-ed
	char** symbollist = backtrace_symbols(addrlist, addrlen);

	// allocate string which will be filled with the demangled function name
	size_t funcnamesize = 256;
	char * funcname = (char*)malloc(funcnamesize);

	// iterate over the returned symbol lines. skip the first, it is the
	// address of this function.
	for (int i = 1; i < addrlen; i++)
	{
		char *begin_name = 0, *begin_offset = 0, *end_offset = 0;

		// find parentheses and +address offset surrounding the mangled name:
		// ./module(function+0x15c) [0x8048a6d]
		for (char *p = symbollist[i]; *p; ++p)
		{
			if (*p == '(')
				begin_name = p;
			else if (*p == '+')
				begin_offset = p;
			else if (*p == ')' && begin_offset) {
				end_offset = p;
				break;
			}
		}

		if (begin_name && begin_offset && end_offset
			&& begin_name < begin_offset)
		{
			*begin_name++ = '\0';
			*begin_offset++ = '\0';
			*end_offset = '\0';

			// mangled name is now in [begin_name, begin_offset) and caller
			// offset in [begin_offset, end_offset). now apply
			// __cxa_demangle():

			int status;
			char* ret = abi::__cxa_demangle(begin_name, funcname, &funcnamesize, &status);

			if (status == 0) {
				funcname = ret; // use possibly realloc()-ed string
				buffer << "  " << symbollist[i] << " : " << funcname << "+" << begin_offset << std::endl;
			}
			else {
				// demangling failed. Output function name as a C function with
				// no arguments.
				buffer << "  " << symbollist[i] << " : " << begin_name << "()+" << begin_offset << std::endl;
			}
		}
		else
		{
			// couldn't parse the line? print the whole line.
			buffer << "  " << symbollist[i] << std::endl;
		}
	}

	free(funcname);
	free(symbollist);
	return buffer.str();
}
#else
std::string get_stacktrace(unsigned int max_frames) {
	return "Stack trace not available on non-linux platforms.";
}
#endif

#ifdef __linux__
void stacktraceHandler(int sig) {
	cout.flush(); // Ensure that the error is printed at the end
	char * signal = strsignal(sig);
	cerr << "Error: signal "<< sig << " - " << signal << endl;
	cerr << get_stacktrace() << endl;
	exit(1);
}
#else
void stacktraceHandler(int sig) { }
#endif