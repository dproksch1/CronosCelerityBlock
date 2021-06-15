#include "KeplerOrbit.H"
#include "Orbit.H"
#include "physical_constants.H"

//KeplerOrbit::KeplerOrbit(double massA, double massB,
//		double semiMajor_Axis, double orbitalEccentricity, double orbitalPhaseOffset) {
//	this->massA = massA;
//	this->massB = massB;
//	this->semiMajor_Axis = semiMajor_Axis;
//	this->orbitalEccentricity = orbitalEccentricity;
//	this->orbitalPhaseOffset = orbitalPhaseOffset;
//
//	orbitalPeriod = get_orbitalPeriod(massA, massB, semiMajor_Axis);
//}


KeplerOrbit::KeplerOrbit(const Quantity &massA, const Quantity &massB, const Quantity &orbitalPeriod,
		const Quantity &semiMajor_Axis, double orbitalEccentricity, double orbitalPhaseOffset) {

	// Assume all parameters to be set in SI units
	this->massA = massA;
	this->massB = massB;
	this->orbitalPeriod = orbitalPeriod;
	this->semiMajor_Axis = semiMajor_Axis;
	this->orbitalEccentricity = orbitalEccentricity;
	this->orbitalPhaseOffset = orbitalPhaseOffset;

	// Check if parameters are set correctly:
	if(! check_orbitalParameters() ) {
		cerr << " Error: orbital parameters are not set consistently " << endl;
		exit(-3);
	}

}

KeplerOrbit KeplerOrbit::newFromPeriod(const Quantity &massA, const Quantity &massB,
		const Quantity &orbitalPeriod, double orbitalEccentricity, double orbitalPhaseOffset) {

	// Compute semi-major axis
	Quantity semiMajor_Axis = get_semiMajorAxis(massA, massB, orbitalPeriod);

	return KeplerOrbit(massA, massB, orbitalPeriod, semiMajor_Axis, orbitalEccentricity, orbitalPhaseOffset);
}


KeplerOrbit KeplerOrbit::newFromSemiMajorAxis(const Quantity &massA, const Quantity &massB, const Quantity &semiMajor_Axis,
		double orbitalEccentricity, double orbitalPhaseOffset) {

	// Compute orbital period
	Quantity orbitalPeriod = get_orbitalPeriod(massA, massB, semiMajor_Axis);

	return KeplerOrbit(massA, massB, orbitalPeriod, semiMajor_Axis, orbitalEccentricity, orbitalPhaseOffset);
}


void KeplerOrbit::get_stellarPositionsCart(normalisation &norm, NumArray<double> &posA_num,
			NumArray<double> &posB_num, double time_num ) {
	// Transform time to SI units
	Quantity time_SI = norm.get_phys(norm.TIME, time_num );
//	cout << " time_SI " << time_SI << endl;

	std::vector<Quantity> posA_SI(3), posB_SI(3);
	get_stellarPositionsCart(posA_SI, posB_SI, time_SI);

	// Cartesian case
	// Transform positions into numerical units
	for(int i_dir=0; i_dir<2; ++i_dir) {
		posA_num[i_dir] = norm.get_num(posA_SI[i_dir]);
		posB_num[i_dir] = norm.get_num(posB_SI[i_dir]);
	}
}

void KeplerOrbit::get_stellarPositionsCart(std::vector<Quantity> &posA,
			std::vector<Quantity> &posB, const Quantity &time) {

	// Quantity zurückgeben

	Quantity rho_A, rho_B, phi;
	get_stellarPositionsBase(rho_A, rho_B, phi, time );

	Quantity xA = cos(phi.get_value()) * rho_A;
	Quantity yA = sin(phi.get_value()) * rho_A;

	Quantity xB = -cos(phi.get_value()) * rho_B;
	Quantity yB = -sin(phi.get_value()) * rho_B;

	posA[0] = xA;
	posA[1] = yA;

	posB[0] = xB;
	posB[1] = yB;
}

void KeplerOrbit::get_stellarPositionsCyl(normalisation &norm, NumArray<double> &posA,
			NumArray<double> &posB, double time_num) {
	// Transform time to SI units
	Quantity time_SI = norm.get_phys(norm.TIME, time_num);

	std::vector<Quantity> posA_SI(3), posB_SI(3);
	get_stellarPositionsCyl(posA_SI, posB_SI, time_SI);

	// Cartesian case
	// Transform positions into numerical units - transform only necessary for radial component
	posA(0) = norm.get_num(posA_SI[0]);
	posA(1) = norm.get_num(posA_SI[1]);

	posB(0) = norm.get_num(posB_SI[0]);
	posB(1) = norm.get_num(posB_SI[1]);
}

void KeplerOrbit::get_stellarPositionsCyl(std::vector<Quantity> &posA,
			std::vector<Quantity> &posB, const Quantity &time_SI) {

	Quantity rho_A_SI, rho_B_SI, phi;
//	double phi;
	get_stellarPositionsBase(rho_A_SI, rho_B_SI, phi, time_SI);

	Quantity phase_shift;
	phase_shift.set_value(M_PI);

	posA[0] = rho_A_SI;
	posA[1] = phi;

	posB[0] = rho_B_SI;
	if(phi<M_PI) {
		posB[1] = phi+phase_shift;
	} else {
		posB[1] = phi-phase_shift;
	}
}


void KeplerOrbit::get_stellarPositionsBase(Quantity &rho_A,
		Quantity &rho_B, Quantity &phi, const Quantity &time) {

	const double orbitalPhase = ORBIT::PHASE::time2Phase(time, orbitalPeriod) + orbitalPhaseOffset;

	const double val_meanAnomaly = ORBIT::PHASE::meanAnomaly(orbitalPhase);
	const double val_trueAnomaly = ORBIT::PHASE::trueAnomaly(orbitalPhase, orbitalEccentricity);
	const double orbitalSeparation =
		ORBIT::PHASE::binarySeparation(orbitalPhase, semiMajor_Axis, orbitalEccentricity);

	const double separationFraction = massA / (massA + massB);
//	cout << " phase " << orbitalPhase << " " << val_meanAnomaly << " " << val_trueAnomaly << endl;

	rho_A = orbitalSeparation * (1. - separationFraction) * CRONOS_CONSTANTS::Meter;
	phi.set_value(val_trueAnomaly);

	rho_B = orbitalSeparation * separationFraction * CRONOS_CONSTANTS::Meter;

}


void KeplerOrbit::get_stellarVelocitiesCart(normalisation &norm, NumArray<double> &veloA,
			NumArray<double> &veloB, double time) {
	// Transform time to SI units
	Quantity time_SI = norm.get_phys(norm.TIME, time);

	std::vector<Quantity> v_A_SI(3), v_B_SI(3);
	get_stellarVelocitiesCart(v_A_SI, v_B_SI, time_SI);

	// Cartesian case
	// Transform positions into numerical units
	veloA(0) = norm.get_num(v_A_SI[0]);
	veloA(1) = norm.get_num(v_A_SI[1]);
	veloB(0) = norm.get_num(v_B_SI[0]);
	veloB(1) = norm.get_num(v_B_SI[1]);
}


void KeplerOrbit::get_stellarVelocitiesCart(std::vector<Quantity> &veloA,
			std::vector<Quantity> &veloB, const Quantity &time_SI) {
	Quantity v_rho_A_SI, v_rho_B_SI, dot_phi_SI;
	get_stellarVelocitiesBase(v_rho_A_SI, v_rho_B_SI, dot_phi_SI, time_SI);

	// Also need positions to calculate cartesian coordinates
	Quantity rho_A_SI, rho_B_SI, phi;
	get_stellarPositionsBase(rho_A_SI, rho_B_SI, phi, time_SI );

	// Cartesian velocity of first star
	Quantity vA_x = cos(phi.get_value())*v_rho_A_SI -
			sin(phi.get_value())*rho_A_SI*dot_phi_SI;
	Quantity vA_y = sin(phi.get_value())*v_rho_A_SI +
			cos(phi.get_value())*rho_A_SI*dot_phi_SI;

	// Cartesian velocity of second star
	Quantity vB_x = -cos(phi.get_value())*v_rho_B_SI +
			sin(phi.get_value())*rho_B_SI*dot_phi_SI;
	Quantity vB_y = -sin(phi.get_value())*v_rho_B_SI -
			cos(phi.get_value())*rho_B_SI*dot_phi_SI;


	veloA[0] = vA_x;
	veloA[1] = vA_y;

	veloB[0] = vB_x;
	veloB[1] = vB_y;

}



void KeplerOrbit::get_stellarVelocitiesCyl(normalisation &norm, NumArray<double> &v_A,
			NumArray<double> &v_B, double time) {
	// Transform time to SI units
	Quantity time_SI = norm.get_phys(norm.TIME, time);

	std::vector<Quantity> v_A_SI(3), v_B_SI(3);
	get_stellarVelocitiesCyl(v_A_SI, v_B_SI, time_SI);

	// Cartesian case
	// Transform positions into numerical units
	v_A[0] = norm.get_num(v_A_SI[0]);
	v_A[1] = norm.get_num(v_A_SI[1]);
	v_B[0] = norm.get_num(v_B_SI[0]);
	v_B[1] = norm.get_num(v_B_SI[1]);
}


void KeplerOrbit::get_stellarVelocitiesCyl(std::vector<Quantity> &veloA,
			std::vector<Quantity> &veloB, const Quantity &time) {
	Quantity v_rho_A, v_rho_B, dot_phi;
	get_stellarVelocitiesBase(v_rho_A, v_rho_B, dot_phi, time);

	veloA[0] = v_rho_A;
	veloA[1] = dot_phi;

	veloB[0] = v_rho_B;
	veloB[1] = dot_phi;
}


void KeplerOrbit::get_stellarVelocitiesBase(Quantity &vRho_A, Quantity &vRho_B,
		Quantity &omega, const Quantity &time) {
	const double orbitalPhase = ORBIT::PHASE::time2Phase(time, orbitalPeriod) + orbitalPhaseOffset;

	const double val_meanAnomaly = ORBIT::PHASE::meanAnomaly(orbitalPhase);
	const double val_trueAnomaly = ORBIT::PHASE::trueAnomaly(orbitalPhase, orbitalEccentricity);

	double v_rel = ORBIT::PHASE::radialVelocity(orbitalPhase, orbitalPeriod,
			semiMajor_Axis, orbitalEccentricity);
	omega = ORBIT::PHASE::trueAngularVelocity(orbitalPhase, orbitalPeriod,
			orbitalEccentricity) / CRONOS_CONSTANTS::Second;

	const double separationFraction = massA / (massA + massB);
	vRho_A = v_rel*(1. - separationFraction) * CRONOS_CONSTANTS::Meter/CRONOS_CONSTANTS::Second;
	vRho_B = v_rel* separationFraction * CRONOS_CONSTANTS::Meter/CRONOS_CONSTANTS::Second;

}



Quantity KeplerOrbit::get_orbitalPeriod(const double massA, const double massB, const double semiMajorAxis) {
	double GGrav = CRONOS_CONSTANTS::GravitationalConstant;

	double semiMajorAxis_SI = semiMajorAxis;
	double massA_SI = massA;
	double massB_SI = massB;

	double orbitalPeriod = 2.*M_PI*sqrt(ORBIT::cube(semiMajorAxis_SI)/(GGrav*(massA_SI + massB_SI)));

	// return orbital period in years
	return orbitalPeriod * CRONOS_CONSTANTS::Second;
}


Quantity KeplerOrbit::get_semiMajorAxis(const double massA, const double massB, const double orbitalPeriod) {
	double GGrav = CRONOS_CONSTANTS::GravitationalConstant;

	double orbitalPeriod_SI = orbitalPeriod;
	double massA_SI = massA;
	double massB_SI = massB;

	double semiMajorAxis = cbrt(sqr(orbitalPeriod_SI/(2.*M_PI))*GGrav*(massA_SI + massB_SI));

	return semiMajorAxis * CRONOS_CONSTANTS::Meter;


}


bool KeplerOrbit::check_orbitalParameters() {

	double eps = 1.e-8;

	double period = get_orbitalPeriod(massA, massB, semiMajor_Axis);

	double err_period = std::abs(period - orbitalPeriod);

	if(err_period < eps*orbitalPeriod) {
		return true;
	} else {
		return false;
	}
}
