#ifndef CRONOS_TIME_STEPPING_H
#define CRONOS_TIME_STEPPING_H = 1


#include "matrix.H"
#include "problem.H"
#include "data.H"
#include "constants.H"
#include <vector>


class TimeIntegrator {
public:
	TimeIntegrator(const int);
	~TimeIntegrator();
	void set_IntRange(int [], int []);
	void set_corrField(int);
	void save_data(const Pot om[], const int substep);
	virtual REAL get_dt(Data & gdata, const int substep);
	virtual void Substep(const Data &gdata, ProblemType &, NumMatrix<REAL, 3> &nom,
	             Pot om[], const int substep) = 0;

protected:
	int ibeg[DIM], iend[DIM];
	int qch;
	const int n_saves;
	NumMatrix<REAL, DIM> *om_save;
};


class RKSteps: public TimeIntegrator {
public:
	RKSteps();
	virtual REAL get_dt(Data & gdata, int substep);
	virtual void Substep(const Data &gdata, ProblemType &, NumMatrix<REAL, 3> &nom,
	             Pot om[], const int substep);

private:
	REAL oneThird, twoThirds;
};


class VanLeerIntegrator: public TimeIntegrator {
public:
	VanLeerIntegrator();
	virtual REAL get_dt(Data & gdata, int substep);
	virtual void Substep(const Data &gdata, ProblemType &, NumMatrix<REAL, 3> &nom,
	             Pot om[], const int substep);

};



#endif
