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
	virtual ~TimeIntegrator() = 0;
	void set_IntRange(int [], int []);
	void set_corrField(int);
	void save_data(const Pot om[], const int substep);
	virtual double get_dt(Data & gdata, const int substep);
	virtual void Substep(const Data &gdata, ProblemType &, NumMatrix<double, 3> &nom,
	             Pot om[], const int substep) = 0;

protected:
	int ibeg[DIM], iend[DIM];
	int qch;
	const int n_saves;
	NumMatrix<double, DIM> *om_save;
};


class RKSteps: public TimeIntegrator {
public:
	RKSteps();
	virtual double get_dt(Data & gdata, int substep);
	virtual void Substep(const Data &gdata, ProblemType &, NumMatrix<double, 3> &nom,
	             Pot om[], const int substep);

private:
	double oneThird, twoThirds;
};


class VanLeerIntegrator: public TimeIntegrator {
public:
	VanLeerIntegrator();
	virtual double get_dt(Data & gdata, int substep);
	virtual void Substep(const Data &gdata, ProblemType &, NumMatrix<double, 3> &nom,
	             Pot om[], const int substep);

};



#endif
