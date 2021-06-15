#ifndef CRONOS_PHYS_CONST_H
#define CRONOS_PHYS_CONST_H

#include "units.H"


//
// SI BASE UNITS
//

namespace CRONOS_BASE_UNITS {

const Quantity Meter 			= UnitQuantity(Quantity::LEN);
const Quantity KiloGram 		= UnitQuantity(Quantity::MASS);
const Quantity Second 			= UnitQuantity(Quantity::TIME);
const Quantity Kelvin 			= UnitQuantity(Quantity::TEMP);
const Quantity Ampere 			= UnitQuantity(Quantity::CURR);

const Quantity Velocity 		= Meter / Second;
const Quantity Acceleration 	= Velocity / Second;
const Quantity Newton 			= Acceleration * KiloGram;
const Quantity Joule 			= Newton * Meter;
const Quantity Watt 			= Joule / Second;
const Quantity Coulomb 			= Ampere * Second;

const Quantity Volt 			= Joule / Coulomb;
const Quantity Tesla 			= KiloGram / ( Coulomb * Second);
}



//
// SOME DERIVED UNITS AND CONSTANTS
//

#ifdef _MSC_VER
// defined here rather than using _USE_MATH_DEFINES before including cmath.h because of header include order
#define M_PI 3.14159265358979323846
#endif

const double   pi						= M_PI;

namespace CRONOS_CONSTANTS {

using namespace CRONOS_BASE_UNITS;

const Quantity SpeedOfLight 			= 299792458.0 * Velocity;
const Quantity VacuumPermeability		= 4.e-7 * pi * KiloGram * Meter / sqr(Coulomb);
const Quantity VacuumPermittivity		= 1. / (VacuumPermeability * sqr(SpeedOfLight));
const Quantity ThomsonCrossSection		= 6.6524587158e-29 * sqr(Meter);
const Quantity GravitationalConstant	= 6.6741e-11 * cube(Meter) / (KiloGram * sqr(Second));
const Quantity BoltzmannConstant		= 1.38064852e-23 * Joule / Kelvin;
const Quantity PlanckConstant			= 6.626070040e-34 * Joule * Second;
const Quantity StefanBoltzmannConstant 	= 5.670367e-8 * Watt / sqr(Meter) / sqr(sqr(Kelvin));
const double   FineStructureConstant	= 7.2973525664e-3;

const Quantity AtomicMassUnit			= 1.660539040e-27 * KiloGram;
const Quantity ElementaryCharge			= 1.602176621e-19 * Coulomb;
const Quantity AstronomicalUnit			= 1.49598e11 * Meter;
const Quantity SolarRadius              = 6.955e8 * Meter;
const Quantity SolarMass				= 1.98855e30 * KiloGram;
const Quantity SolarLuminosity			= 3.83e26 * Watt;

const Quantity ElectronVolt				= ElementaryCharge * Volt;
const Quantity Parsec					= 3.0856775814913672789e16 * Meter;
const Quantity Year						= 31536000.0 * Second;
const Quantity Gauss					= 1e-4 * Tesla;
const Quantity Erg						= 1e-7 * Joule;

const Quantity ElectronMass				= 9.1093836e-31 * KiloGram;
const Quantity ProtonMass				= 1.6726219e-27 * KiloGram;
const Quantity HydrogenMass				= ElectronMass + ProtonMass;
const Quantity CrossSection_pp			= 3.2e-30 * sqr(Meter);
}



//
// SOME CONSTANTS IN SI
// This is needed for some legacy code to work.
// These constants should not be used for up-to-date projects
// Work with the above quantities instead!
//

namespace CRONOS_LEGACY_CONSTANTS {

const double kB = 1.3807e-23;
const double mu0 = 4.*pi*1.e-7;
const double parsec = 3.08568025e16;
const double AstronomicalUnit = 149.6e9; // 1 Astronomical Unit
const double G_Grav = 6.673e-11; // Gravitational constant
const double AtomicMassUnit = 1.66053878e-27;
const double m_H = 1.00794*AtomicMassUnit; // Hydrogen mass
const double m_e = 5.489e-4*AtomicMassUnit; // Electron mass
const double m_p = m_H - m_e; // Proton mass (binding energy ignored)
const double year = 31556952; // one year in units of seconds (365,2425 days)
const double stefan =  5.670373e-8;  // [J m⁻² s⁻¹  K⁻⁴ ] Stefan-Boltzmann constant
const double echarge = 1.60217657e-19; // [C] elementary charge
const double c_light = 2.99792458e8; // [m s⁻¹] speed of light
const double sigma_Th = 6.6524586e-29; // [m²]
const double sigma_pp = 3.2e-30; // [m²]
const double MeV_to_Joule = 1.60217733e-13;
const double fineStructure = 7.297e-3;
const double r_el = 2.817960e-15; // [m] classical electron radius  
const double Planck = 6.626070e-34; // [Js]
const double eps_0 = 8.8541878176e-12; // [F m^-1]

// Solar system units
const double MSun = 1.9891e30; // Solar mass
const double LSun = 3.846e26; // [J/s] =  [kg m² s⁻³]
const double RSun = 6.955e8; // [m]
}







#endif
