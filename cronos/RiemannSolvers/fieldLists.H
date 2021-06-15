#ifndef CRONOS_FIELDLISTS_H
#define CRONOS_FIELDLISTS_H = 1

#include <string>
#include <map>
#include "constants.H"
#include <stdio.h>
#include <assert.h>
#include "fluid.H"
#include "multifluid.H"


//#if (FLUID_TYPE == CRONOS_MHD)
//const int N_OM_PAR = 1;
//const int N_OM_PERP = 2;
//#else
//const int N_OM_PAR = 0;
//const int N_OM_PERP = 0;
//#endif
//const int N_OM_NORM = N_OMINT - N_OM_PERP - N_OM_PAR;

class fieldLists {
public:
//#if(FLUID_TYPE == CRONOS_MULTIFLUID)
//	fieldLists(const CronosMultifluid &fluids, int);
//#else
	fieldLists(const CronosFluid &fluid, int);
//#endif
//	fieldLists(int);
	std::string getGenType(int q);
	bool isGenType(int q, const std::string &);
	bool isGenTypeCentered(int q);
	bool isGenTypeParallel(int q);
	bool isGenTypePerpendicular(int q);
private:
	std::map<int, std::string> GenType;
	std::map<int, std::string> Name;
	std::map<int, std::string> SpecificType;
//	int N_OM_PAR;
//	int N_OM_
};


#endif
