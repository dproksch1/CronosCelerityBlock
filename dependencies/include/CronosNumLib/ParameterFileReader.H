#ifndef PARAMETER_FILE_READER_H
#define PARAMETER_FILE_READER_H

#include <cstring>
#include <cstdlib>
#include <cerrno>

#include <sstream>
#include <fstream>
#include <string>
#include <unordered_map>
#include <iostream>
#include <vector>
#include <limits>
#include <algorithm>

#include "util.H"


class ParameterFileReader {
public:
	typedef std::unordered_map<std::string, std::string>  paramMapT;

	ParameterFileReader();
	ParameterFileReader(const std::string fname);
	~ParameterFileReader();

	// parameter access
	bool value_exists(const std::string key);
	double value(const std::string key);
	std::string svalue(const std::string key);

	const paramMapT & get_all();


	// dynamic versions (triggering reload)
	bool value_exists_dynamic(const std::string key);
	double value_dynamic(const std::string key);
	std::string svalue_dynamic(const std::string key);

	// Trigger file readin
	bool read_file();

	void set(const std::string key, const std::string value);

protected:
	// internal parameter access
	bool value_exists(const std::string key, paramMapT & paramMap);
	double value(const std::string key, paramMapT & paramMap);
	std::string svalue(const std::string key, paramMapT & paramMap);

	// read in
	void get_envFName();
	bool read_file(paramMapT & _paramMap);

	// control flags
	bool use_envFName;
	bool paramsLoaded_static;

	// paramter storage
	std::string fname;
	paramMapT paramMap_static;
	paramMapT paramMap_dynamic;
};

extern ParameterFileReader ProjParameters;
#endif
