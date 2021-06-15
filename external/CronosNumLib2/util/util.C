#include "util.H"
#include <iostream>
#include <stdio.h>

#include <string.h>
#include <stdlib.h>
#include <errno.h>


using namespace std;


bool cronos_isnan(double val)
{
	return val != val;
}


ParameterFileReader ProjParameters;

bool value_exists(const char *str, const char* ext)
{
	return ProjParameters.value_exists(string(str));
}

double value(const char *str, const char* ext)
{
	return ProjParameters.value(string(str));
}

string svalue(const char *str, const char* ext)
{
	return ProjParameters.svalue(string(str));
}

bool value_exists(std::string str, const char* ext) {
	return ProjParameters.value_exists(str);
}

double value(std::string str, const char* ext) {
	return ProjParameters.value(str);
}

std::string svalue(std::string str, const char* ext) {
	return ProjParameters.svalue(str);
}


std::string ltrim(const std::string& s)
{
	size_t start = s.find_first_not_of(WHITESPACE);
	return (start == std::string::npos) ? "" : s.substr(start);
}

std::string rtrim(const std::string& s)
{
	size_t end = s.find_last_not_of(WHITESPACE);
	return (end == std::string::npos) ? "" : s.substr(0, end + 1);
}

std::string trim(const std::string& s)
{
	return rtrim(ltrim(s));
}
