#ifndef CRONOS_NORMALISATION_H
#define CRONOS_NORMALISATION_H



#ifndef MODULE_TEST
#include "constants.H"
#include <util.H>
#endif

#include "units.H"
#include "physical_constants.H"
#include "CException.H"

#include <iostream>
#include <string>
#include <cassert>
#include <cmath>
#include <gsl/gsl_linalg.h>

/*********************************************************\
 *
 * CLASS: normalisation:
 *
 * Supplies routines to convert between
 * physical (in SI units) and numerical quantities
 *
\*********************************************************/

class normalisation
{
public:
	static const int N_BASE = 6;
	typedef enum {
		// primary types
		// these have to remain always in the same order and first in the list
		LEN,	// m
		MASS,	// kg
		TIME,	// s
		TEMP,	// K
		CURR,	// A
		NUM_DENS,	// part m^(-3)

		// derived
		CHAR,	// C
		DENS,	// kg m^(-3)
		VEL,	// m s^(-1)
		PRESS,	// kg m^(-1) s^(-2)
		FORCE,  // kg m s^(-2)
		ACCEL,  // m s^(-2)

		MOMENTUM, // kg m s^(-1)
		MOM_DENS, // kg m^(-2) s^(-1)

		ENERGY, // J
		E_DENS,  // J m^(-3)

		MAG_IND, // T
		MAG_VEC_POT, // V s/m <-> m T

		SRC_E,   // kg m^(-1)s^(-3)
		SRC_MOM, // kg m^(-2)s^(-2)

		SIGMA,	// kg m^(-2)
		CROSS_SEC,	// m^(2)
		CROSS_SEC_S, // m^(2) kg^{-1}
		GTimesM, // m^3 s^(-2)

		RATE, 		// s^(-1)
		DIFF_EMISS, // W m^(-3) sr^(-1) J^(-1)
		DIFF_E_FLUX, // W m^(-2) J^(-1)

		MASS_LOSS_RATE, // kg s^{-1}

		COOLING_RATE, // kg m^5 s^(-3)

		N_TYPES
	} Type;

	normalisation();
	virtual ~normalisation();

	void set_parent(
			normalisation const * parent_norm
			);

	virtual void setup() = 0;

	inline const Quantity & unitQuantity(
			const Type type
			) const;

	inline const Quantity & normQuantity(
			const Type type
			) const;

	inline Quantity get_phys(
			const Type & i,
			const double & num
			) const ;

//	inline double get_num(
//			const Type & i,
//			const double & phys
//			) const ;

	inline double get_num(
			const Quantity & quantity
			) const ;

	inline double get_num(
			const Type & i,
			const Quantity & quantity
			) const ;

	inline double renormalise(
			const Type & type,
			const normalisation * norm,
			const double & num
			) const ;

	inline double renormalise(
			const Type & old_type,
			const normalisation * norm,
			const Type & new_type,
			const double & num
			) const ;

	inline double renormalise(
			const Type & type,
			const normalisation & norm,
			const double & num
			) const ;

	inline double renormalise(
			const Type & old_type,
			const normalisation & norm,
			const Type & new_type,
			const double & num
			) const ;


	// For Info output
	void print_norms(
			std::ostream& os
			);

	// For HDF5
	double get_normConst(			// return normalisation constant to transform back to phys units
			std::string type
			);

	std::string get_unitNormConst(
			std::string type
			);


protected:
	std::string units[N_TYPES];

	void set_definition(
			const Type type
			);

	void set_definition(
			const Type type,
			const Quantity & quantity
			);

	void set_definitionFromParent(
			const Type type
			);

	const Quantity & get_parentDefinition(
			const Type type
			) const;

	void finish_setup();


private:
	void _set_unitDims(
			const Type i,
			const std::initializer_list<double> unitDimension
			);

	void _get_derivedUnits(
			Quantity quantities[N_TYPES]
			);

	void _get_baseUnitsFromDefinition(
			);

	Type _physicalTypesMap[Quantity::N_PHYSICAL_DIMENSIONS];
	normalisation const * _parent_norm;


	Quantity _unitQuantity[N_TYPES];
	double _unitDims[N_TYPES][N_BASE];

	bool _norm_is_defined(
			const Type type
			) const;

	Quantity _normQuantity[N_TYPES];
	Type _definedType[N_BASE];
	int _currentDefinition;
	bool _setup_finished;

};

class NoNormalisation : public normalisation {
public:
	virtual void setup();
};

#ifdef CRONOS_CONSTANTS_H_INCLUDED
class NormSound : public normalisation {
public:
	virtual void setup();
};

class NormMagInd : public normalisation {
public:
	virtual void setup();
};
#endif





//
//	Inline Method definitions
//

inline const Quantity & normalisation::unitQuantity(
		const Type type
		) const
{
	return _unitQuantity[type];
}

inline const Quantity & normalisation::normQuantity(
		const Type type
		) const
{
	assert(_norm_is_defined(type));
	return _normQuantity[type];
}

inline Quantity normalisation::get_phys(
		const Type & i,
		const double & num
		) const
{
	return num * _normQuantity[i];
}

//inline double normalisation::get_num(
//		const Type & i,
//		const double & phys
//		) const
//{
//	return phys / _normQuantity[i].get_value();
//}

inline double normalisation::get_num(
		const Quantity & quantity
		) const
{
	double value = quantity.get_value();
	for( int iDim=0; iDim < Quantity::N_PHYSICAL_DIMENSIONS; ++iDim) {
		value *= pow(
					_normQuantity[_physicalTypesMap[iDim]].get_value(),
					- quantity.get_unitDimension(iDim)
				);
	}
	return value;
}

inline double normalisation::get_num(
		const Type & i,
		const Quantity & quantity
		) const
{
	return quantity.normalizeTo(_normQuantity[i]);
}

inline double normalisation::renormalise(
		const Type & type,
		const normalisation * norm,
		const double & num
		) const
{
	return renormalise(type, *norm, type, num);
}

inline double normalisation::renormalise(
		const Type & old_type,
		const normalisation * norm,
		const Type & new_type,
		const double & num
		) const
{
	return renormalise(old_type, *norm, new_type, num);
}

inline double normalisation::renormalise(
		const Type & type,
		const normalisation & norm,
		const double & num
		) const
{
	return renormalise(type, norm, type, num);
}

inline double normalisation::renormalise(
		const Type & old_type,
		const normalisation & norm,
		const Type & new_type,
		const double & num
		) const
{
	return get_num(new_type, norm.get_phys(old_type, num));
}

#endif
