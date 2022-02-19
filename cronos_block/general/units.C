#include "units.H"
using namespace std;



Quantity::Quantity(
		)
{
	this->value = 0;

	for ( int i=0; i < N_PHYSICAL_DIMENSIONS; ++i) {
		this->unitDimension[i] = 0;
	}
}

Quantity::Quantity(
		const double value,
		const double unitDimension[N_PHYSICAL_DIMENSIONS]
		)
{
	this->value = value;

	for ( int i=0; i < N_PHYSICAL_DIMENSIONS; ++i) {
		this->unitDimension[i] = unitDimension[i];
	}
}

Quantity::Quantity(
		const double value,
		const std::initializer_list<double> unitDimension
		)
{
	assert(unitDimension.size() == N_PHYSICAL_DIMENSIONS);

	this->value = value;
	std::copy(unitDimension.begin(), unitDimension.end(), this->unitDimension);
}

const double * Quantity::get_unitDimension(
		) const
{
	return unitDimension;
}

const double & Quantity::get_unitDimension(
		const int & idx
		) const
{
	assert(idx >= 0 && idx < N_PHYSICAL_DIMENSIONS);
	return unitDimension[idx];
}


bool Quantity::sameUnitDimension(
		const Quantity & otherQuant
		) const
{
	bool sameDimension = true;
	for ( int i=0; i < N_PHYSICAL_DIMENSIONS; ++i) {
		sameDimension &= (unitDimension[i] == otherQuant.unitDimension[i]);
	}
	return sameDimension;
}

bool Quantity::isNumber(
		) const
{
	for ( int i=0; i < N_PHYSICAL_DIMENSIONS; ++i) {
		if (unitDimension[i] == 0)
			continue;
		else
			return false;
	}
	return true;
}

std::string Quantity::get_si_units(
		) const
{
	stringstream stream;
	stream << std::fixed << std::setprecision(1);

	std::string symb_loc;

	for ( int i=0; i < N_PHYSICAL_DIMENSIONS; ++i) {
		if (unitDimension[i] == 0)
			continue;

		switch(i) {
			case LEN:  symb_loc = "m" ; break;
			case MASS: symb_loc = "kg"; break;
			case TIME: symb_loc = "s" ; break;
			case TEMP: symb_loc = "K" ; break;
			case CURR: symb_loc = "A" ; break;
		}

		if (unitDimension[i] == 1) {
			stream << symb_loc;
		} else {
			int exp = unitDimension[i];

			if (double(exp) == unitDimension[i]) {
				stream << symb_loc << "^" << exp;
			}
			else {
				stream << symb_loc << "^" << unitDimension[i];
			}
		}
		stream << " ";
	}

	return stream.str();
}

Quantity & operator +=(
		Quantity & qL,
		const Quantity & qR
		)
{
	assert(qL.sameUnitDimension(qR));

	qL.value += qR.value;
	return qL;
}

Quantity & operator -=(
		Quantity & qL,
		const Quantity & qR
		)
{
	assert(qL.sameUnitDimension(qR));

	qL.value -= qR.value;
	return qL;
}

Quantity & operator *=(
		Quantity & qL,
		const Quantity & qR
		)
{
	for ( int i=0; i < Quantity::N_PHYSICAL_DIMENSIONS; ++i) {
		qL.unitDimension[i] += qR.unitDimension[i];
	}
	qL.value *= qR.value;
	return qL;
}

Quantity & operator /=(
		Quantity & qL,
		const Quantity & qR
		)
{
	for ( int i=0; i < Quantity::N_PHYSICAL_DIMENSIONS; ++i) {
		qL.unitDimension[i] -= qR.unitDimension[i];
	}
	qL.value /= qR.value;
	return qL;
}

Quantity & operator *=(
		Quantity & q,
		const double & x
		)
{
	q.value *= x;
	return q;
}

Quantity & operator /=(
		Quantity & q,
		const double & x
		)
{
	q.value /= x;
	return q;
}


Quantity operator+(
		const Quantity & qL,
		const Quantity & qR
		)
{
	Quantity newQuant(qL);
	newQuant += qR;
	return newQuant;
}

Quantity operator-(
		const Quantity & qL,
		const Quantity & qR
		)
{
	Quantity newQuant(qL);
	newQuant -= qR;
	return newQuant;
}

Quantity operator*(
		const Quantity & qL,
		const Quantity & qR
		)
{
	Quantity newQuant(qL);
	newQuant *= qR;
	return newQuant;
}

Quantity operator/(
		const Quantity & qL,
		const Quantity & qR
		)
{
	Quantity newQuant(qL);
	newQuant /= qR;
	return newQuant;
}

Quantity operator*(
		const Quantity & q,
		const double & x
		)
{
	Quantity newQuant(q);
	newQuant *= x;
	return newQuant;
}

Quantity operator/(
		const Quantity & q,
		const double & x
		)
{
	Quantity newQuant(q);
	newQuant /= x;
	return newQuant;
}

Quantity operator*(
		const double & x,
		const Quantity & q
		)
{
	return q * x;
}

Quantity operator/(
		const double & x,
		const Quantity & q
		)
{
	Quantity newQuant;
	for ( int i=0; i < Quantity::N_PHYSICAL_DIMENSIONS; ++i) {
		newQuant.unitDimension[i] = -q.unitDimension[i];
	}
	newQuant.value = x / q.value;
	return newQuant;
}


std::ostream& operator<<(
		std::ostream& os,
		const Quantity& quantity
		)
{
    if (!os.good())
        return os;

    os << quantity.value;
//	for ( int i=0; i < Quantity::N_PHYSICAL_DIMENSIONS; ++i) {
//		os << quantity.unitDimension[i] << "\t";
//	}
//    os << "|";
    os << " " << quantity.get_si_units();

    return os;
}


bool operator>(
		const double & x,
		const Quantity & q
		)
{
	return x > q.value;
}

bool operator>(
		const Quantity & q,
		const double & x
		)
{
	return q.value > x;
}


bool operator>(
		const Quantity & qL,
		const Quantity & qR
		)
{
	return qL.value > qR.value;
}


bool operator<(
		const double & x,
		const Quantity & q
		)
{
	return x < q.value;
}

bool operator<(
		const Quantity & q,
		const double & x
		)
{
	return q.value < x;
}

bool operator<(
		const Quantity & qL,
		const Quantity & qR
		)
{
	return qL.value < qR.value;
}



Quantity sqr(
		const Quantity & q
		)
{
	return q * q;
}

Quantity cube(
		const Quantity & q
		)
{
	return q * q * q;
}

Quantity pow(
		const Quantity & q,
		const double & exp
		)
{
	Quantity newQuant(q);
	for ( int i=0; i < Quantity::N_PHYSICAL_DIMENSIONS; ++i) {
		newQuant.unitDimension[i] *= exp;
	}
	newQuant.value = std::pow(q.value, exp);
	return newQuant;
}

Quantity sqrt(
		const Quantity & q
		)
{
	return pow(q, 0.5);
}


UnitQuantity::UnitQuantity(
		):
				Quantity()
{
	value = 1;
}

UnitQuantity::UnitQuantity(
		const int id
		):
				Quantity()
{
	value = 1;
	unitDimension[id] = 1;
}

UnitQuantity::UnitQuantity(
		const double unitDimension[Quantity::N_PHYSICAL_DIMENSIONS]
		):
				Quantity(1., unitDimension)
{
}

UnitQuantity::UnitQuantity(
		const std::initializer_list<double> unitDimension
		):
				Quantity(1., unitDimension)
{
}

UnitQuantity::UnitQuantity(
		const Quantity & quantity
		):
				Quantity(quantity)
{
	this->value = 1.;
}

