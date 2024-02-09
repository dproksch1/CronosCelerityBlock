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
                      double dt, const int substep, size_t nom_max[]) = 0;
	

protected:
	int ibeg[DIM], iend[DIM];
	int qch;
	const int n_saves;
	std::vector<CelerityBuffer<double, 3>> omSYCL_save;
};


class RKSteps: public TimeIntegrator {
public:
	RKSteps();
	virtual void Substep(Queue &queue, const Data &gdata, CelerityRange<3> omRange,
                      CelerityBuffer<nom_t, 3> nomSYCL,
                      double dt, const int substep, size_t nom_max[]);

private:
	double oneThird, twoThirds;
};


class VanLeerIntegrator: public TimeIntegrator {
public:
	VanLeerIntegrator();
	virtual void Substep(Queue &queue, const Data &gdata, CelerityRange<3> omRange,
                      CelerityBuffer<nom_t, 3> nomSYCL,
                      double dt, const int substep, size_t nom_max[]);

};



#endif
