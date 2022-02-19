#ifndef CRONOS_UNITS_H
#define CRONOS_UNITS_H

#include <iomanip>
#include <sstream>
#include <iostream>
#include <string>
#include <cmath>
#include <cassert>
#include <initializer_list>

#ifndef double
//#define REAL double
#endif


class Quantity
{
public:
	typedef enum {
		// primary (in parameter file)
		LEN,	// m
		MASS,	// kg
		TIME,	// s
		TEMP,	// K
		CURR,	// A

		N_PHYSICAL_DIMENSIONS
	} physDim;

	Quantity(
			);

	Quantity(
			const double value,
			const double unitDimension[N_PHYSICAL_DIMENSIONS]
			);

	Quantity(
			const double value,
			const std::initializer_list<double> unitDimension
			);

	inline double normalizeTo(
			const Quantity & quantity
			) const;

	inline const Quantity & set_value(
			const double value
			);

	inline const double & get_value(
			) const;

	const double * get_unitDimension() const;
	const double & get_unitDimension(const int & idx) const;


	std::string get_si_units() const;

	bool sameUnitDimension(const Quantity & otherQuant) const;
	bool isNumber() const;

	// Typecasting
	template<typename T> operator T() const {return T(value);}

protected:
	double value;
	double unitDimension[N_PHYSICAL_DIMENSIONS];


private:
	friend std::ostream& operator<<(std::ostream & os, const Quantity & quantity);

	friend Quantity & operator+=(Quantity & qL, const Quantity & qR);
	friend Quantity & operator-=(Quantity & qL, const Quantity & qR);
	friend Quantity & operator*=(Quantity & qL, const Quantity & qR);
	friend Quantity & operator/=(Quantity & qL, const Quantity & qR);

	friend Quantity operator+(const Quantity & qL, const Quantity & qR);
	friend Quantity operator-(const Quantity & qL, const Quantity & qR);
	friend Quantity operator*(const Quantity & qL, const Quantity & qR);
	friend Quantity operator/(const Quantity & qL, const Quantity & qR);

	friend Quantity & operator*=(Quantity & q, const double & x);
	friend Quantity & operator/=(Quantity & q, const double & x);

	friend Quantity operator*(const Quantity & q, const double & x);
	friend Quantity operator/(const Quantity & q, const double & x);
	friend Quantity operator*(const double & x, const Quantity & q);
	friend Quantity operator/(const double & x, const Quantity & q);

	friend bool operator>(const double & x, const Quantity & q);
	friend bool operator>(const Quantity & q, const double & x);
	friend bool operator>(const Quantity & qL, const Quantity & qR);
	friend bool operator<(const double & x, const Quantity & q);
	friend bool operator<(const Quantity & q, const double & x);
	friend bool operator<(const Quantity & qL, const Quantity & qR);

	friend Quantity sqr(const Quantity & q);
	friend Quantity cube(const Quantity & q);
	friend Quantity sqrt(const Quantity & q);
	friend Quantity pow(const Quantity & q, const double & exp);

};

class UnitQuantity: public Quantity {
public:
	UnitQuantity(
			);
	UnitQuantity(
			const Quantity & quantity
			);
	UnitQuantity(
			const int id
			);
	UnitQuantity(
			const double unitDimension[N_PHYSICAL_DIMENSIONS]
			);
	UnitQuantity(
			const std::initializer_list<double> unitDimension
			);
};



//
//	Inline Method definitions
//

inline double Quantity::normalizeTo(
		const Quantity & quantity
		) const
{
	assert(this->sameUnitDimension(quantity));
	return value / quantity.value;
}

inline const Quantity & Quantity::set_value(
		const double value
		)
{
	this->value = value;
	return * this;
}

inline const double & Quantity::get_value(
		) const
{
	return value;
}


#endif
