#include "Gravity.H"
#include "physical_constants.H"

Gravity::Gravity(normalisation &norm, std::vector<Quantity> posObject_SI, Quantity MassObject_SI) {
	Quantity GPhys = CRONOS_CONSTANTS::GravitationalConstant;

	// Compute normalised values from physical quantities
	GTimesM_num = norm.get_num(norm.GTimesM, GPhys*MassObject_SI);
	for (unsigned i_dir=0; i_dir < posObject_SI.size(); i_dir++) {
		posObject_num[i_dir] = norm.get_num(posObject_SI[i_dir]);
	}
}

Gravity::Gravity(normalisation &norm, Quantity MassObject_SI) {
	Quantity GPhys = CRONOS_CONSTANTS::GravitationalConstant;

	// Compute normalised values from physical quantities
	GTimesM_num = norm.get_num(norm.GTimesM, GPhys*MassObject_SI);
	posObject_num.clear();
}


double Gravity::get_ForceAbs(const NumArray<double> &pos_num) {
	return get_ForceAbs(pos_num, posObject_num);
}

double Gravity::get_ForceAbs(const NumArray<double> &pos_num,
		const NumArray<double> &posObj_num) {
	// Get distance to object
	double sqr_dist = sqr(pos_num[0] - posObj_num[0]) +
			sqr(pos_num[1] - posObj_num[1]) +
			sqr(pos_num[2] - posObj_num[2]);
//	cout << " o0 " << posObj_num[0] << endl;
//	cout << " o1 " << posObj_num[1] << endl;
//	cout << " o2 " << posObj_num[2] << endl;
//
//	cout << " obj0 " << posObject_num[0] << endl;
//	cout << " obj1 " << posObject_num[1] << endl;
//	cout << " obj2 " << posObject_num[2] << endl;

	// Compute gravitational force in numerical units
	double FGravAbs = GTimesM_num/sqr_dist;

	return FGravAbs;
}

void Gravity::get_Force(const NumArray<double> &pos_num, NumArray<double> &force_num) {
	get_Force(pos_num, posObject_num, force_num);
}

void Gravity::get_Force(const NumArray<double> &pos_num,
		const NumArray<double> &posObj_num, NumArray<double> &force_num) {
	// Get distance to object
	double dist = sqrt(sqr(pos_num[0] - posObj_num[0]) +
			sqr(pos_num[1] - posObj_num[1]) +
			sqr(pos_num[2] - posObj_num[2]));
	double factor = GTimesM_num/cube(dist);

	force_num[0] = (posObj_num[0] - pos_num[0])*factor;
	force_num[1] = (posObj_num[1] - pos_num[1])*factor;
	force_num[2] = (posObj_num[2] - pos_num[2])*factor;
}
