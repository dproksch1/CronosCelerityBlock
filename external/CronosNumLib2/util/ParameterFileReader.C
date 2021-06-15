#include "ParameterFileReader.H"

using namespace std;


ParameterFileReader::ParameterFileReader()
{
	// do not read file at construction
	// env could not be set yet

	use_envFName = true;
	paramsLoaded_static = false;
}


ParameterFileReader::ParameterFileReader(const string fname)
{
	this->fname = fname;
	use_envFName = false;
	paramsLoaded_static = false;
}


ParameterFileReader::~ParameterFileReader()
{
}


void ParameterFileReader::get_envFName()
{
	char _fname[256];
	if (getenv("poub") == NULL) {
		fprintf(stderr,"Environment variable %s not defined.\n","poub");
		exit(1);
	}
	strcpy(_fname, getenv("poub"));
	if (getenv("pname") == NULL) {
		fprintf(stderr,"Environment variable %s not defined.\n","pname");
		exit(1);
	}
	sprintf(_fname,"%s/%s.%s",getenv("poub"),getenv("pname"), "cat");
	fname = _fname;
}

bool ParameterFileReader::read_file()
{
	paramsLoaded_static = read_file(paramMap_static);
	return paramsLoaded_static;
}


bool BothAreSpaces(char lhs, char rhs)
{
	return (lhs == rhs) && (lhs == ' ');
}

bool ParameterFileReader::read_file(paramMapT & _paramMap)
{
	if (use_envFName) {
		get_envFName();
	}

	ifstream file(fname.c_str());

	if (file.fail()) {
		cerr << "value: cannot open file " << fname << endl;
		cerr << "   - error: " << strerror(errno) << endl;
		exit(1);
	} else {

		// clear all existing parameters (in case of a reload)
		_paramMap.clear();

		string line;
		istringstream iss_line;
		string key, value;

		vector<string> parts;
		string part;
		string cleanedLine;

		while( getline(file, line) )
		{
			// Trim left whitespaces
			line = ltrim(line);

			// jump empty lines
			// not doing this can lead to problems
			if (line.size() == 0) {
				continue;
			}

			// Jump comment lines
			if (line[0] == '#') {
				continue;
			}

			// Remove double whitespaces (needed for the parameter separation)
			std::string::iterator new_end = std::unique(line.begin(), line.end(), BothAreSpaces);
			line.erase(new_end, line.end());

			// Jump lines that were containing only whitespaces
			if (line.length() == 0) {
				continue;
			}

			iss_line.str(line);
			iss_line.clear();	


			// parameter separation
			// this is to account for whitespaces within the parameter line
			// only the first parameter (separated by whitespaces) after the '=' is used
			//
			parts.clear();
			while( getline(iss_line, part, ' ') ) {
				parts.push_back(part);
			}
			cleanedLine = parts[0];
			for (int i=1; i < parts.size(); ++i) {
				cleanedLine += parts[i];
				if (parts[i-1].find('=') != std::string::npos) { // contains '='
					break;
				}
			}
			iss_line.str(cleanedLine);
			iss_line.clear();	
			// parameter separation end
		

			getline(iss_line, key, '=');
			getline(iss_line, value, ' ');
			
			if ( _paramMap.count(key) > 0 ) {
				cerr << "variable " << key << " already in the dictionary" << endl;
				exit(1);
			} else {
				_paramMap[key] = value;
			}
		}
	}

	file.close();

	return true;
}


bool ParameterFileReader::value_exists(const std::string key, paramMapT & paramMap)
{
	return (paramMap.count(key) > 0);
}


double ParameterFileReader::value(const std::string key, paramMapT & paramMap)
{
	if (value_exists(key, paramMap)) {
		try{
			return stod(paramMap[key]);
		} catch(std::invalid_argument & e) {
			stringstream err;
			err << "ParameterFileReader: Error converting parameter " << key << "=" << paramMap[key] << " to double";
			cerr << err.str() << endl;
			throw invalid_argument(err.str());
		}
	} else {
		stringstream err;
		err << "ParameterFileReader: key not found: " << key;
		cerr << err.str() << endl;
		throw invalid_argument(err.str());
	}
}

string ParameterFileReader::svalue(const std::string key, paramMapT & paramMap)
{

	if (value_exists(key, paramMap)) {
		return paramMap[key];
	} else {
		stringstream err;
		err << "ParameterFileReader: key not found: " << key;
		cerr << err.str() << endl;
		throw invalid_argument(err.str());
	}
}

void ParameterFileReader::set(const std::string key, const std::string value) {
#ifndef NDEBUG
	if (!value_exists(key, paramMap_static)) {
		cerr << "variable " << key << " not existent - it was be created" << endl;
	}
#endif
	paramMap_static[key] = value;
}


bool ParameterFileReader::value_exists(const std::string key)
{
	// load if not loaded already
	if (!paramsLoaded_static) {
		paramsLoaded_static = read_file(paramMap_static);
	}

	return value_exists(key, paramMap_static);
}


double ParameterFileReader::value(const std::string key)
{
	// load if not loaded already
	if (!paramsLoaded_static) {
		paramsLoaded_static = read_file(paramMap_static);
	}

	return value(key, paramMap_static);
}

string ParameterFileReader::svalue(const std::string key)
{
	// load if not loaded already
	if (!paramsLoaded_static) {
		paramsLoaded_static = read_file(paramMap_static);
	}

	return svalue(key, paramMap_static);
}

bool ParameterFileReader::value_exists_dynamic(const std::string key)
{
	read_file(paramMap_dynamic);
	return value_exists(key, paramMap_dynamic);
}

double ParameterFileReader::value_dynamic(const std::string key)
{
	read_file(paramMap_dynamic);
	return value(key, paramMap_dynamic);
}

string ParameterFileReader::svalue_dynamic(const std::string key)
{
	read_file(paramMap_dynamic);
	return svalue(key, paramMap_dynamic);
}



const ParameterFileReader::paramMapT & ParameterFileReader::get_all()
{
	return paramMap_static;
}

