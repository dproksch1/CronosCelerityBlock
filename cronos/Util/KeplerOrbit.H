#ifndef CRONOS_KEPLER_ORBIT_H
#define CRONOS_KEPLER_ORBIT_H

#include "normalisation.H"
#include "units.H"
#include <matrix.H>

class KeplerOrbit {
public :
	/** Constructor for Kepler orbit class
	 * assuming all variables to be given in SI units
	 * */
//	KeplerOrbit(double massA, double massB, double semiMajor_Axis,
//			double orbitalEccentricity, double orbitalPhaseOffset=0.);
	KeplerOrbit(const Quantity &massA, const Quantity &massB, const Quantity &orbitalPeriod,
			const Quantity &semiMajor_Axis, double orbitalEccentricity, double orbitalPhaseOffset=0.);

	/** Constructor for Kepler orbit class
	 * assuming all variables to be given in SI units
	 * Here, masses, eccentricity and semi-major axis are provided, while orbital period is computed
	 * */
	static KeplerOrbit newFromPeriod(const Quantity &massA, const Quantity &massB, const Quantity &semiMajor_Axis,
			double orbitalEccentricity, double orbitalPhaseOffset=0.);

	/** Constructor for Kepler orbit class
	 * assuming all variables to be given in SI units
	 * Here, masses, eccentricity and orbital period are provided, while semi-major axis is computed
	 * */
	static KeplerOrbit newFromSemiMajorAxis(const Quantity &massA, const Quantity &massB, const Quantity &orbitalPeriod,
			double orbitalEccentricity, double orbitalPhaseOffset=0.);

	/** compute orbital period from masses and semi-major axis
	 * @PARAM massA mass of first star (units of Solar masses)
	 * @PARAM massB mass of second star (units of Solar masses)
	 * @PARAM semiMajorAxis semi-major axis of the system (units of AU)
	 * */
	static Quantity get_orbitalPeriod(const double massA,
			const double massB, const double semiMajorAxis);

	/** compute semi-major axis from masses and orbital period
	 * @PARAM massA mass of first star (units of Solar masses)
	 * @PARAM massB mass of second star (units of Solar masses)
	 * @PARAM orbitalPeriod orbital period of the system (units of s)
	 * */
	static Quantity get_semiMajorAxis(const double massA,
			const double massB, const double orbitalPeriod);

	/** compute Cartesian position of stars in orbital plane with
	 * periostron along x-axis for given time (SI units)
	 * @return Cartesian position in SI units*/
	void get_stellarPositionsCart(std::vector<Quantity> &posA,
			std::vector<Quantity> &posB, const Quantity &time);
	/** compute Cartesian position of stars in orbital plane with
	 * periostron along x-axis for given time (units determined by normalisation)
	 * @return Cartesian position in units determined by normalisation */
	void get_stellarPositionsCart(normalisation &norm, NumArray<double> &posA,
			NumArray<double> &posB, double time);

	/** compute cylindrical position of stars in orbital plane with
	 * periostron along x-axis for given time (SI units)
	 * @return Cylindrical position in SI units*/
	void get_stellarPositionsCyl(std::vector<Quantity> &posA,
			std::vector<Quantity> &posB, const Quantity &time);
	/** compute cylindrical position of stars in orbital plane with
	 * periostron along x-axis for given time (units determined by normalisation)
	 * @return Cylindrical position in units determined by normalisation */
	void get_stellarPositionsCyl(normalisation &norm, NumArray<double> &posA,
			NumArray<double> &posB, double time);

	/** compute Cartesian velocities of stars in orbital plane with
	 * periostron along x-axis for given time (SI units)
	 * @return Cartesian position in SI units*/
	void get_stellarVelocitiesCart(std::vector<Quantity> &posA,
			std::vector<Quantity> &posB, const Quantity &time);
	/** compute Cartesian position of stars in orbital plane with
	 * periostron along x-axis for given time (units determined by normalisation)
	 * @return Cartesian position in units determined by normalisation */
	void get_stellarVelocitiesCart(normalisation &norm, NumArray<double> &posA,
			NumArray<double> &posB, double time);

	/** compute cylindrical velocities of stars in orbital plane with
	 * periostron along x-axis for given time (SI units)
	 * @return Cylindrical velocities in SI units*/
	void get_stellarVelocitiesCyl(std::vector<Quantity> &posA,
			std::vector<Quantity> &posB, const Quantity &time);
	/** compute cylindrical velocities of stars in orbital plane with
	 * periostron along x-axis for given time (units determined by normalisation)
	 * @return Cylindrical velocities in units determined by normalisation */
	void get_stellarVelocitiesCyl(normalisation &norm, NumArray<double> &posA,
			NumArray<double> &posB, double time);

private :
	/** compute position of stars in orbital plane with
	 * periostron along x-axis for given time (SI units)
	 * @return position in SI units*/
	void get_stellarPositionsBase(Quantity &rho_A, Quantity &rho_B,
			Quantity &phi, const Quantity &time);
	/** Compute orbital velocity of stars in cylindrical coordinates
	 * */
	void get_stellarVelocitiesBase(Quantity &vRho_A, Quantity &vRho_B,
			Quantity &omega, const Quantity &time);

	/** Check if orbital parameters are set correctly */
	bool check_orbitalParameters();

	/** masses of both stars in units of kg */
	double massA, massB;
	/** semi-major axis of stellar system in units of m */
	double semiMajor_Axis;
	double orbitalEccentricity;
	double orbitalPhaseOffset;
	/** orbital period of system given / computed in s */
	double orbitalPeriod;
};



#endif
