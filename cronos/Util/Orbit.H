#ifndef ORBIT_H
#define ORBIT_H

#include <gsl/gsl_errno.h>
#include <gsl/gsl_roots.h>
#include <gsl/gsl_const_mksa.h>

/**
 * Example Usage:
 * Binary system in the x - y plane, periastron along x axis
 *
 * Given: massA, massB, orbitalPeriod, orbitalEccentricity, semiMajorAxis
 *
		const double orbitalPhase = ORBIT::PHASE::time2Phase(gdata.time, orbitalPeriod) + orbitalPhaseOffset;

		const double meanAnomaly = ORBIT::PHASE::meanAnomaly(orbitalPhase);
		const double trueAnomaly = ORBIT::PHASE::trueAnomaly(orbitalPhase, orbitalEccentricity); const double orbitalSeparation =
		ORBIT::PHASE::binarySeparation(orbitalPhase, semiMajorAxis_num, orbitalEccentricity);

		const double separationFraction = massA / (massA + massB);

		xA = cos(trueAnomaly) * orbitalSeparation * (1. - separationFraction);
		yA = sin(trueAnomaly) * orbitalSeparation * (1. - separationFraction);

		xB = -cos(trueAnomaly) * orbitalSeparation * separationFraction;
		yB = -sin(trueAnomaly) * orbitalSeparation * separationFraction;

 * */

namespace ORBIT {

double sqr(const double x) { return x * x; }

double cube(const double x) { return x * x * x; }


/** compute the mean anomaly from a given eccentric anomaly */
double meanAnomaly(const double eccAnomaly, const double eccentricity) {
	return (eccAnomaly - eccentricity * sin(eccAnomaly));
}

/** compute the true anomaly from a given eccentric anomaly */
double trueAnomaly(const double eccAnomaly, const double eccentricity) {
	const double eccAnomaly_cos = cos(eccAnomaly);
	double trueAnomaly = acos((eccAnomaly_cos - eccentricity) /
			(1. - eccentricity * eccAnomaly_cos));
	const int orbitSegment = floor(eccAnomaly / M_PI);
	const int orbitNumber = floor(0.5 * (orbitSegment + 1));
	if (orbitSegment % 2 == 0) {
		trueAnomaly = orbitNumber * 2 * M_PI + trueAnomaly;
	} else {
		trueAnomaly = orbitNumber * 2 * M_PI - trueAnomaly;
	}
	return trueAnomaly;
}

/** compute the binary separation from a given eccentric anomaly */
double binarySeparation(const double eccAnomaly, const double semiMajorAxis,
		const double eccentricity) {
	return semiMajorAxis * (1. - eccentricity * cos(eccAnomaly));
}

namespace DETAILS {
double meanAnomaly_df(const double eccAnomaly, const double eccentricity) {
	return (1 - eccentricity * cos(eccAnomaly));
}

double trueAnomaly_df(const double eccAnomaly, const double eccentricity) {
	return sqrt((1. - sqr(eccentricity))) / (1. - eccentricity * cos(eccAnomaly));
}

double trueAnomaly_dMeanAnomaly(const double eccAnomaly,
		const double eccentricity) {
	return DETAILS::trueAnomaly_df(eccAnomaly, eccentricity) /
			DETAILS::meanAnomaly_df(eccAnomaly, eccentricity);
}
}  // namespace DETAILS

/** compute the eccentric anomaly from a given mean anomaly */
double eccentricAnomaly_mean(const double meanAnomaly,
		const double eccentricity) {
	double eccAnomaly = meanAnomaly;

	int MAX_ITERATIONS = 100;
	double tol = 1e-10;
	int iter = 0;

	for (iter = 0; iter < MAX_ITERATIONS; ++iter) {
		const double f = ORBIT::meanAnomaly(eccAnomaly, eccentricity) - meanAnomaly;
		const double df = ORBIT::DETAILS::meanAnomaly_df(eccAnomaly, eccentricity);
		const double dx = -f / df;
		eccAnomaly += dx;
		if (abs(dx) < abs(tol * meanAnomaly) || abs(f) < tol) {
//			cout << " break after " << iter << " " << dx << " " << tol*eccAnomaly << " " << f << endl;
			break;
		}
	}

	return eccAnomaly;
}

/** compute the eccentric anomaly from a given true anomaly */
double eccentricAnomaly_true(const double trueAnomaly,
		const double eccentricity) {
	double eccAnomaly = trueAnomaly;

	int MAX_ITERATIONS = 100;
	double tol = 1e-10;
	int iter = 0;

	for (int iter = 0; iter < MAX_ITERATIONS; ++iter) {
		const double f = ORBIT::trueAnomaly(eccAnomaly, eccentricity) - trueAnomaly;
		const double df = ORBIT::DETAILS::trueAnomaly_df(eccAnomaly, eccentricity);
		const double dx = -f / df;
		eccAnomaly += dx;
		if (abs(dx) < abs(tol * eccAnomaly) || abs(f) < tol) {
			break;
		}
	}

	return eccAnomaly;
}

/** Convenient user functions in terms of phase (0..1)*/
namespace PHASE {

double eccAnomaly2Phase(const double eccAnomaly, const double eccentricity) {
	return ORBIT::meanAnomaly(eccAnomaly, eccentricity) / (2. * M_PI);
}

double trueAnomaly2Phase(const double trueAnomaly, const double eccentricity) {
	const double eccAnomaly =
			ORBIT::eccentricAnomaly_true(trueAnomaly, eccentricity);
	return eccAnomaly2Phase(eccAnomaly, eccentricity);
}

double time2Phase(const double time, const double period) {
	return time / period;
}

double meanAnomaly(const double phase) { return phase * 2. * M_PI; }

double eccentricAnomaly(const double phase, const double eccentricity) {
	return ORBIT::eccentricAnomaly_mean(PHASE::meanAnomaly(phase), eccentricity);
}

double trueAnomaly(const double phase, const double eccentricity) {
	const double eccAnomaly = PHASE::eccentricAnomaly(phase, eccentricity);
	return ORBIT::trueAnomaly(eccAnomaly, eccentricity);
}

double trueAngularVelocity(const double phase, const double period,
		const double eccentricity) {
	const double eccAnomaly = PHASE::eccentricAnomaly(phase, eccentricity);
	return DETAILS::trueAnomaly_dMeanAnomaly(eccAnomaly, eccentricity) * 2. *
			M_PI / period;
}

/** Compute the radial velocity between both components for given orbital phase*/
double radialVelocity(const double phase, const double period,
		const double semiMajorAxis, const double eccentricity) {
	const double eccAnomaly = PHASE::eccentricAnomaly(phase, eccentricity);
	return 2.*M_PI*semiMajorAxis*eccentricity/(period)*sin(eccAnomaly)/(1-eccentricity*cos(eccAnomaly));
}


/** compute the binary separation for a given orbtial phase*/
double binarySeparation(const double phase, const double semiMajorAxis,
		const double eccentricity) {
	const double eccAnomaly = PHASE::eccentricAnomaly(phase, eccentricity);
	return ORBIT::binarySeparation(eccAnomaly, semiMajorAxis, eccentricity);
}
}  // namespace PHASE




}  // namespace ORBIT





#endif
