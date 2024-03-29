/**
 * @file
 * @brief Time integration and cell update applications
 *
 * This file provides time integration functionality that allows for the application of cell updates
 * determined by the main kernel.
 */

#ifndef CRONOS_TIME_STEPPING_H
#define CRONOS_TIME_STEPPING_H = 1


#include "matrix.H"
#include "problem.H"
#include "data.H"
#include "constants.H"
#include "utils.H"
#include <vector>

class TimeIntegrator {
public:
	TimeIntegrator(const int);
	virtual ~TimeIntegrator() = 0;
	void set_IntRange(int [], int []);
	void init_omBuffer(Queue &, const int []);
	void set_corrField(int);
	virtual void Substep(Queue &queue, const Data &gdata, CelerityRange<3> omRange,
                      CelerityBuffer<nom_t, 3> nomSYCL,
                      double dt, const int substep) = 0;
	

protected:
	int ibeg[DIM], iend[DIM];
	const int n_saves;
	std::vector<CelerityBuffer<double, 3>> omSYCL_save;
};


class RKSteps: public TimeIntegrator {
public:
	RKSteps();
	virtual void Substep(Queue &queue, const Data &gdata, CelerityRange<3> omRange,
                      CelerityBuffer<nom_t, 3> nomSYCL,
                      double dt, const int substep);

private:
	double oneThird, twoThirds;
};


class VanLeerIntegrator: public TimeIntegrator {
public:
	VanLeerIntegrator();
	virtual void Substep(Queue &queue, const Data &gdata, CelerityRange<3> omRange,
                      CelerityBuffer<nom_t, 3> nomSYCL,
                      double dt, const int substep);

};



#endif
