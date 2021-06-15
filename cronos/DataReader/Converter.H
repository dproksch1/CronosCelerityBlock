#include <sys/stat.h>
#include <iostream>
#include <dirent.h>
//#include "EndianType.H" // needed ???
#include <matrix.H>
#include "Hdf5File_cbase.H"
#include <vector>
#include <cstdlib>    // for atoi
#include <algorithm>

using std::string;
using std::cout;
using std::cerr;

class Attribute {
public:
	Attribute (string attr_name);
	template<class T>
	void fill_n (int size, T value[], int type);
	void fill_s (int size, string value);
	short type;
	string name;
	int size;
	string * svalue;
	int    * ivalue;
	float  * fvalue;
	double * dvalue;
};

class Converter{
public:
	Converter (string poub, string pname);
	void convert (int itime);
	void get_itimeList (std::vector<int> & it_list);
	bool findfname (string &, int, int, int, int);
	bool verbose;
private:
	void ReadHDF (string);
	void WriteHDF (int);
	int GetNumFields (string);
	NumMatrix<float,3>* om;
	std::vector<Attribute> attr_list;
	int mx[3], nproc[3];
	int numFields;
	int rim;
	double xb[3], dx[3];
	string poub, pname;
	string *omNames;
	bool EdgeGridding;
};
