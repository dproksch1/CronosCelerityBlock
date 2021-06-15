#ifndef CRONOS_GRAVITY_H
#define CRONOS_GRAVITY_H

#include "normalisation.H"
#include "units.H"
#include <matrix.H>
#include <vector>


class Gravity {
public:
	Gravity(normalisation &norm, std::vector<Quantity> posObject_SI, Quantity MassObject_SI);

	Gravity(normalisation &norm, Quantity MassObject_SI);
//	Gravity(normalisation &norm, NumArray<double> posObject_num, double MassObject_num);

	/** Compute absolute value of gravitational force (normalised units)
	 * */
	double get_ForceAbs(const NumArray<double> &pos_num);

	/** Compute absolute value of gravitational force (normalised units)
	 * Here, position of object is provided dynamically
	 * */
	double get_ForceAbs(const NumArray<double> &pos_num, const NumArray<double> &posObj_num);

	/** Compute gravitational force vector (normalised units)
	 * */
	void get_Force(const NumArray<double> &pos_num, NumArray<double> &force_num);

	/** Compute gravitational force vector (normalised units)
	 * Here, position of object is provided dynamically
	 * */
	void get_Force(const NumArray<double> &pos_num, const NumArray<double> &posObj_num, NumArray<double> &force_num);
private:
	/** Position of massive object in normalised units */
	NumArray<double> posObject_num;
	/** Mass times G in normalised units */
	double GTimesM_num;
};


#endif
