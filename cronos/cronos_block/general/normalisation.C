#include "normalisation.H"


using namespace std;
using namespace CRONOS_CONSTANTS;

normalisation::normalisation()
{
	for (int iQuant=0; iQuant < N_TYPES; ++iQuant) {
		for (int iDim=0; iDim < N_BASE; ++iDim) {
			_unitDims[iQuant][iDim] = 0;
		}
	}

	_parent_norm = 0;
	_currentDefinition = 0;
	_setup_finished = false;

	_physicalTypesMap[Quantity::LEN]  = LEN;
	_physicalTypesMap[Quantity::MASS] = MASS;
	_physicalTypesMap[Quantity::TIME] = TIME;
	_physicalTypesMap[Quantity::TEMP] = TEMP;
	_physicalTypesMap[Quantity::CURR] = CURR;


	// store units of normalisation constants
//	units[LEN] = "m";
//	units[MASS] = "kg";
//	units[NUM_DENS] = "m^{-3}";
//	units[TEMP] = "K";
//	units[CURR] = "A";
//	units[CHAR] = "C";
//	units[TIME] = "s";
//	units[ENERGY] = "J";
//	units[DENS] = "kg m^{-3}";
//	units[VEL] = "m/s";
//	units[MAG_IND] = "T";
//	units[MAG_VEC_POT] = "V s m^{-1}";
//	units[E_DENS] = "J m^{-3}";
//	units[SRC_E] = "kg m^(-1)s^(-3)";
//	units[SRC_MOM] = "kg m^(-2)s^(-2)";
//	units[PRESS] = "kg m^(-1) s^(-2)";
//	units[FORCE] = "kg m s^(-2)";
//	units[ACCEL] = "m s^(-2)";
//	units[SIGMA] = "kg m^(-2)";
//	units[GTimesM] = "m^3 s^(-2)";
//	units[MOMENT] = "kg m s^(-1)";
//	units[MOM_DENS] = "kg m^(-2) s^(-1)";


	_unitQuantity[LEN]  = Meter;
	_unitQuantity[MASS] = KiloGram;
	_unitQuantity[TIME] = Second;
	_unitQuantity[TEMP] = Kelvin;
	_unitQuantity[CURR] = Ampere;
	_unitQuantity[NUM_DENS] = 1. / cube(Meter);


	// 							LEN		MASS	TIME	TEMP	CURR	NUM_DENS
	_set_unitDims(LEN,			{1, 	0, 		0,		0,		0,		0		});
	_set_unitDims(MASS,			{0, 	1, 		0,		0,		0,		0		});
	_set_unitDims(TIME,			{0, 	0, 		1,		0,		0,		0		});
	_set_unitDims(TEMP,			{0, 	0, 		0,		1,		0,		0		});
	_set_unitDims(CURR,			{0, 	0, 		0,		0,		1,		0		});
	_set_unitDims(NUM_DENS,		{0, 	0, 		0,		0,		0,		1		});

	_set_unitDims(VEL,			{1, 	0, 		-1,		0,		0,		0		});
	_set_unitDims(ACCEL,		{1, 	0, 		-2,		0,		0,		0		});
	_set_unitDims(FORCE,		{1, 	1, 		-2,		0,		0,		0		});
	_set_unitDims(ENERGY,		{2, 	1, 		-2,		0,		0,		0		});
	_set_unitDims(E_DENS,		{2, 	1, 		-2,		0,		0,		1		});
	_set_unitDims(SRC_E,		{2, 	1, 		-3,		0,		0,		1		});
	_set_unitDims(MOMENTUM,		{1, 	1, 		-1,		0,		0,		0		});
	_set_unitDims(MOM_DENS,		{1, 	1, 		-1,		0,		0,		1		});
	_set_unitDims(SRC_MOM,		{1, 	1, 		-2,		0,		0,		1		});
	_set_unitDims(DENS,			{0, 	1, 		0,		0,		0,		1		});
	_set_unitDims(SIGMA,		{1, 	1, 		0,		0,		0,		1		});
	_set_unitDims(CROSS_SEC,		{-1, 	0, 		0,		0,		0,		-1		});
	_set_unitDims(CROSS_SEC_S,		{-1, 	-1, 		0,		0,		0,		-1		});
	_set_unitDims(PRESS,		{2, 	1, 		-2,		0,		0,		1		});
	_set_unitDims(GTimesM,		{3, 	0, 		-2,		0,		0,		0		});
	_set_unitDims(CHAR,			{0, 	0, 		1,		0,		1,		0		});
	_set_unitDims(MAG_IND,		{0, 	1, 		-2,		0,		-1,		0		});
	_set_unitDims(MAG_VEC_POT,	{1, 	1, 		-2,		0,		-1,		0		});
	_set_unitDims(RATE,			{0, 	0, 		-1,		0,		0,		0		});
	_set_unitDims(DIFF_EMISS,	{0, 	0, 		-1,		0,		0,		1		});
	_set_unitDims(DIFF_E_FLUX,	{-2, 	0, 		-1,		0,		0,		0		});
	_set_unitDims(MASS_LOSS_RATE,	{3, 	1, 		-1,		0,		0,		1		});
	_set_unitDims(COOLING_RATE,	{2, 	1, 		-3,		0,		0,		-1		});


	_get_derivedUnits(_unitQuantity);

	for (int iQuant=0; iQuant < N_TYPES; ++iQuant) {
		units[iQuant] = _unitQuantity[iQuant].get_si_units();
	}

}

void normalisation::_set_unitDims(
	const Type i,
	const std::initializer_list<REAL> unitDimension
	)
{
	// Defines a derived unit (with type i) with respect to the base units
	assert(unitDimension.size() == N_BASE);
	std::copy(unitDimension.begin(), unitDimension.end(), _unitDims[i]);
}

normalisation::~normalisation(
		)
{
}

void normalisation::set_parent(
		normalisation const * parent_norm
		)
{
	_parent_norm = parent_norm;
}

void normalisation::print_norms(std::ostream& os) {
	os<<" --- normalization constants ---"		<< endl;
	os<<" Len:  \t"	<<_normQuantity[LEN]		<< endl;
	os<<"       \t"	<<_normQuantity[LEN] / AstronomicalUnit << "au" << endl;
	os<<" Mass: \t"	<<_normQuantity[MASS]		<< endl;
	os<<"       \t"	<<_normQuantity[MASS] * sqr(SpeedOfLight)/ElectronVolt << " eV / c^2" << endl;
	os<<" Time: \t"	<<_normQuantity[TIME]		<< endl;
	os<<"       \t"	<<_normQuantity[TIME] / Year <<"yr" << endl;
	os<<" Temp: \t"	<<_normQuantity[TEMP]		<< endl;
	os<<" Curr: \t"	<<_normQuantity[CURR]		<< endl;
	os<<" NDens:\t"	<<_normQuantity[NUM_DENS]	<< endl;
	os<<" ..............................."<<endl;
	os<<" Mag:  \t"	<<_normQuantity[MAG_IND] / Tesla << "T" << endl;
	os<<" EDens:\t"	<<_normQuantity[E_DENS]	/ (Joule / cube(Meter)) << " J / m^3" << endl;
	os<<"       \t"	<<_normQuantity[E_DENS]	/ (ElectronVolt / cube(Meter)) << " eV / m^3" << endl;
	os<<" -------------------------------"<<endl;
}

REAL normalisation::get_normConst(string type) {
	// return normalisation constant according to field type
	if(type=="rho") {
		return _normQuantity[DENS];
	} else if(type=="v_x" || type=="v_y" || type=="v_z") {
		return _normQuantity[VEL];
	} else if(type=="B_x" || type=="B_y" || type=="B_z") {
		return _normQuantity[MAG_IND];
	} else if(type=="A_x" || type=="A_y" || type=="A_z") {
		return _normQuantity[MAG_VEC_POT];
	} else if(type=="Etherm") {
		return _normQuantity[E_DENS];
	} else if(type=="Temp") {
		return _normQuantity[TEMP];
	} else {
		return UnitQuantity();
	}
}


string normalisation::get_unitNormConst(std::string type) {
	if(type=="rho") {
		return units[DENS];
	} else if(type=="v_x" || type=="v_y" || type=="v_z") {
		return units[VEL];
	} else if(type=="B_x" || type=="B_y" || type=="B_z") {
		return units[MAG_IND];
	} else if(type=="A_x" || type=="A_y" || type=="A_z") {
		return units[MAG_VEC_POT];
	} else if(type=="Etherm") {
		return units[E_DENS];
	} else if(type=="Temp") {
		return units[TEMP];
	} else {
		return "";
	}
}


void normalisation::set_definition(
		const Type type
		)
{
	set_definition(type, _unitQuantity[type]);
}

void normalisation::set_definition(
		const Type type,
		const Quantity & quantity
		)
{
	// Defines the normalisation for a desired type

	assert(quantity.sameUnitDimension(_unitQuantity[type]) && " Definition type / unit mismatch ");
	assert(quantity > 0. && " Only positive normalisations allowed ");

	_definedType[_currentDefinition] = type;
	_normQuantity[type] = quantity;

	_currentDefinition++;
}

void normalisation::set_definitionFromParent(
		const Type type
		)
{
	set_definition(type, get_parentDefinition(type));
}

const Quantity & normalisation::get_parentDefinition(
		const Type type
		) const
{
	assert( _parent_norm != 0 );
	assert( _parent_norm->_setup_finished);

	return _parent_norm->_normQuantity[type];
}


void normalisation::finish_setup() {
	if (_currentDefinition != N_BASE)
		throw CException(" Invalid amount of definitions ");

	_get_baseUnitsFromDefinition();
	_get_derivedUnits(_normQuantity);
	_setup_finished = true;
}




bool normalisation::_norm_is_defined(
		const Type type
		) const
{
	if (_setup_finished) {
		return true;
	} else {
		bool is_defined = false;
		for (int iQuant=0; iQuant < N_BASE; ++iQuant) {
			is_defined |= type == _definedType[iQuant];
		}
		return is_defined;
	}
}


void normalisation::_get_derivedUnits(
		Quantity quantities[N_TYPES]
		)
{
	// Computes the derived quantities from base quantities

	for (int iQuant=N_BASE; iQuant < N_TYPES; ++iQuant) {
		quantities[iQuant] = UnitQuantity(); // Set it to 1
		for (int iBase=0; iBase < N_BASE; ++iBase) {
			quantities[iQuant] *= pow(quantities[iBase], _unitDims[iQuant][iBase]);
		}
	}
}


void normalisation::_get_baseUnitsFromDefinition(
		)
{
	// Computes the normalisation factors for the base units
	// from a set of given normalisations

	REAL _unitDimsDef[N_BASE*N_BASE];
	REAL b_data[N_BASE];
	REAL x_data[N_BASE];

	for (int iQuant=0; iQuant < N_BASE; ++iQuant) {
		Type definedtype = _definedType[iQuant];
		b_data[iQuant] = std::log(_normQuantity[definedtype].get_value());

		for (int iBase=0; iBase < N_BASE; ++iBase) {
			_unitDimsDef[iQuant * N_BASE + iBase] = _unitDims[definedtype][iBase];
		}
	}

	gsl_matrix_view A = gsl_matrix_view_array(_unitDimsDef, N_BASE, N_BASE);
	gsl_vector_view b = gsl_vector_view_array(b_data, N_BASE);
	gsl_vector_view x = gsl_vector_view_array(x_data, N_BASE);

	int s;

	gsl_permutation * p = gsl_permutation_alloc (N_BASE);

	gsl_linalg_LU_decomp (&A.matrix, p, &s);
	gsl_linalg_LU_solve (&A.matrix, p, &b.vector, &x.vector);
	gsl_permutation_free (p);


	for (int iQuant=0; iQuant < N_BASE; ++iQuant) {
		_normQuantity[iQuant] = std::exp(x_data[iQuant]) * _unitQuantity[iQuant];
	}
}



void NoNormalisation::setup() {
	set_definition(LEN);
	set_definition(MASS);
	set_definition(TIME);
	set_definition(TEMP);
	set_definition(CURR);
	set_definition(NUM_DENS);

	finish_setup();
}


#ifdef CRONOS_CONSTANTS_H_INCLUDED
void NormSound::setup(
		)
{
	// Standard normalisation scheme.
	// Set 6 linearly independent normalisations to define a global normalisation.
	// Inside this function only normQuantity(type) of an already specified type is defined.

	REAL L_phys = value("LenPhys");
	REAL L_norm(1.);
	if(value_exists("LenScale")) {
		L_norm = value("LenScale");
	}

	REAL m_phys = value("MolecularMass");
	REAL numDens_phys = value("NumDensPhys");
	REAL Temp_phys = value("TemperaturePhys");
	REAL gamma = value("Adiabatic_exponent");

	set_definition(LEN, L_phys/L_norm * Meter);
	set_definition(MASS, m_phys * HydrogenMass);
	set_definition(NUM_DENS, numDens_phys / cube(Meter));

	set_definition(TEMP, Temp_phys * Kelvin);

	// Velocity - here: speed of sound
	set_definition(VEL, sqrt(gamma * BoltzmannConstant * normQuantity(TEMP) / normQuantity(MASS)));
	set_definition(MAG_IND, sqrt(VacuumPermeability * normQuantity(NUM_DENS) * normQuantity(MASS)) * normQuantity(VEL));

	finish_setup();
}


void NormMagInd::setup(
		)
{
	// Standard normalisation scheme.
	// Set 6 linearly independent normalisations to define a global normalisation.
	// Inside this function only normQuantity(type) of an already specified type is defined.

	REAL L_phys = value("LenPhys");
	REAL L_norm(1.);
	if(value_exists("LenScale")) {
		L_norm = value("LenScale");
	}

	REAL m_phys = value("MolecularMass");
	REAL numDens_phys = value("NumDensPhys");
	REAL gamma = value("Adiabatic_exponent");
	REAL Mag_phys = value("MagIndPhys");

	set_definition(LEN, L_phys/L_norm * Meter);
	set_definition(MASS, m_phys * HydrogenMass);
	set_definition(NUM_DENS, numDens_phys / cube(Meter));

	set_definition(MAG_IND, Mag_phys * Tesla);

	// Velocity - here: speed of sound
	set_definition(VEL, normQuantity(MAG_IND) / sqrt(VacuumPermeability * normQuantity(NUM_DENS) * normQuantity(MASS)));
	set_definition(TEMP, normQuantity(MASS) * sqr(normQuantity(VEL)) / (gamma * BoltzmannConstant));

	finish_setup();
}
#endif
