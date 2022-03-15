#ifndef CRONOS_GRIDPARAM_H
#define CRONOS_GRIDPARAM_H = 1

#include "constants.H"

class Interface_GridFunction {
public:
	Interface_GridFunction() {}
	virtual double get_gridFunc(double) = 0;
	virtual ~Interface_GridFunction() {}
};


class GridFunction_Exp: public Interface_GridFunction {
public:
	GridFunction_Exp(double, double, double);
	virtual double get_gridFunc(double);
	double domain_begin, domain_end, domain_len;
};


class GridFunction_Ratio: public Interface_GridFunction {
public:
	GridFunction_Ratio(double, double, double);
	virtual double get_gridFunc(double);
	double _ref, _min_ratio, _max_ratio;
};

class GridFunction_Percent: public Interface_GridFunction {
public:
	GridFunction_Percent(double, double, double, int);
	virtual double get_gridFunc(double);
private:
	double _domain_beg;
	double _delta, _increase;
	int _num;
};

class GridFunction_PercentVar: public Interface_GridFunction {
public:
	GridFunction_PercentVar(double, double, double, int);
	virtual double get_gridFunc(double);
private:
	double _domain_beg;
	double _delta, _increase;
	int _num;
};



class GridFunction_Log: public Interface_GridFunction {
public:
	GridFunction_Log(double, double);
	virtual double get_gridFunc(double);
	double domain_begin, domain_end;
};


class GridFunction_Lin: public Interface_GridFunction {
public:
	GridFunction_Lin() {}
	virtual double get_gridFunc(double ratio) {return ratio;}
};


#endif
