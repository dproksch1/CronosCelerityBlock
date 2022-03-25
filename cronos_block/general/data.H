#ifndef CRONOS_DATA_H
#define CRONOS_DATA_H = 1

#include <map>
#include "matrix.H"
#ifdef _MSC_VER
#include "timewrapper.H"
#else
#include <time.h>
#endif

#include "Hdf5File_cbase.H"
#include "grid.H"
#include "fluid.H"
#include "constants.H"
#include "queue.H"

#ifndef DATA_TYPE
#define DATA_TYPE double
#define REAL_FLAG
#endif

//#define REAL double

class Pot : public NumMatrix<DATA_TYPE,DIM>
{
public:
	Pot();
	Pot(const int *);
	Pot(const NumMatrix<DATA_TYPE,DIM>&);
	void resize(const int *);
	void resize(const int *, int);
	void resize_ll(const int *, const int *);
	void set_max(const double &);
	void set_min(const double &);
	Pot& operator=(const double&);  
	double get_max();
	double get_min();
// 	std::string get_FieldType();
// 	void set_FieldType(std::string);
// private:
// 	std::string FieldType; // Indicates the name used in the program
};

class id_handler {
public:
	id_handler();
	unsigned int get_fieldId(NumMatrix<double,3> &omField);
	unsigned int get_userFieldId(NumMatrix<double,3> &omField);
	bool is_generic(NumMatrix<double,3> &omField);
	void enter_fieldId(std::string);
	void print_map();
	std::map<std::string, unsigned int> fieldIds;
	unsigned int max_generic;
	unsigned int max_generic_int;
	unsigned int max_user;
};

class Data: public Grid {
public:
	Data();
	~Data();
	void floatom_out(int R, Hdf5Stream & file, int om_lo, int om_hi);
	void float_out(int, bool);
	void CheckNan(int);
	double computeInt(int);
	double computeRMS(int);
	NumMatrix<double,3>  computeAbs(int, int) const;
	double getMin(int);
	double getMin(int [], int);
	double getMax(int);
	double getMax(int [], int);

	void set_fieldIds();
	bool is_userField(NumMatrix<double,3> &omField);

	NumMatrix<float,3> float_data(int, int, bool);
	CronosFluid fluid;

	Pot *om;
	std::vector<CelerityBuffer<double, 3>> omSYCL;
	NumMatrix<double,3> *nom;
//	Pot om[N_OM+N_P];
//	Pot om[31];
//	std::vector<Pot> om;
//	std::vector<Pot*> om;
#if (OMS_USER == TRUE)
	Pot *om_user;
#endif

	// This is set to NULL for OMS_USER == false (has to exist due to design decisions)
	NumMatrix<double,3> *nom_user;

#if (GEOM == SPHERICAL)
	NumMatrix<double,1> massFlux;
	double radialArea;
#endif

	NumMatrix<double,3> pTherm; // can store thermal pressure
	NumMatrix<int,3> carbuncleFlag; // flag, where carbuncle problem might occur

	double time, dt, t_end;
	double cfl;
	int tstep;

	bool mag, storePressure, use_carbuncleFlag;

private:
	/**
	 * Display time in user friendly fashion
	 * @PARAM t_sec full time in seconds
	 * */
	std::string git_humanReadable(double t_sec);

	/** Start and stop time of program for performance analysis
	 * */
	timeval tick_start,tock_end;
	id_handler om_ids;
};

#endif