#ifndef CRONOS_GRIDPARAM_H
#define CRONOS_GRIDPARAM_H = 1

#ifdef parallel
#include "mpi.h"
#endif
#include "constants.H"

class Interface_GridFunction {
public:
	Interface_GridFunction() {}
	virtual REAL get_gridFunc(REAL) = 0;
	virtual ~Interface_GridFunction() {}
};


class GridFunction_Exp: public Interface_GridFunction {
public:
	GridFunction_Exp(REAL, REAL, REAL);
	virtual REAL get_gridFunc(REAL);
	REAL domain_begin, domain_end, domain_len;
};


class GridFunction_Ratio: public Interface_GridFunction {
public:
	GridFunction_Ratio(REAL, REAL, REAL);
	virtual REAL get_gridFunc(REAL);
	REAL _ref, _min_ratio, _max_ratio;
};

class GridFunction_Percent: public Interface_GridFunction {
public:
	GridFunction_Percent(REAL, REAL, REAL, int);
	virtual REAL get_gridFunc(REAL);
private:
	REAL _domain_beg;
	REAL _delta, _increase;
	int _num;
};

class GridFunction_PercentVar: public Interface_GridFunction {
public:
	GridFunction_PercentVar(REAL, REAL, REAL, int);
	virtual REAL get_gridFunc(REAL);
private:
	REAL _domain_beg;
	REAL _delta, _increase;
	int _num;
};



class GridFunction_Log: public Interface_GridFunction {
public:
	GridFunction_Log(REAL, REAL);
	virtual REAL get_gridFunc(REAL);
	REAL domain_begin, domain_end;
};


class GridFunction_Lin: public Interface_GridFunction {
public:
	GridFunction_Lin() {}
	virtual REAL get_gridFunc(REAL ratio) {return ratio;}
};


#endif
