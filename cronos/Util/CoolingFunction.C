#include "CoolingFunction.H"
#include "buildinfo.H"
#include "physical_constants.H"

#include <fstream>

using namespace CRONOS_BASE_UNITS;

CoolingFunction::CoolingFunction(normalisation &norm, ProblemType &problem, int choice_curve) : ref_norm(norm) {

	Quantity cooling_rate_phys = 1.*KiloGram*pow(Meter, 5)*pow(Second, -3);
	cooling_rate_factor_phys2norm = norm.get_num(norm.COOLING_RATE, cooling_rate_phys);

//	cout << " Teste norm " << norm.get_unitNormConst("rho") << endl;
//	cout << " Teste ref " << ref_norm.get_unitNormConst("rho") << endl;

	Quantity energy_dens = 1.*KiloGram*pow(Meter, -1)*pow(Second, -2);
	double energy_num = norm.get_num(norm.E_DENS, energy_dens);

	// no minimum temperature as default value
	this->minTemp_phys = 0.;

	if(choice_curve==1) {
		load_CoolingTable_SchureEtAl2009(norm);
		load_CoolingTable_DalgarnoMcCray1972(norm);
	} else {
		cerr << " currently only cooling curve by Schure et al. (2009) available " << endl;
		exit(3);
	}

	this->gamma = problem.gamma;

//	double XVal = 1./1.4;
//	double YVal = 1.-XVal;
//	double ZVal = 1.-XVal - YVal;

	double XVal = problem.massFraction_X;
	double YVal = problem.massFraction_Y;
	double ZVal = 1.-XVal - YVal;

	//	Quantity mu = CRONOS_CONSTANTS::AtomicMassUnit/(2.*XVal + 3.*(1.-XVal -ZVal)/4. + ZVal/2.);
	Quantity mu = problem.meanParticleMass;
	Quantity mu_e = CRONOS_CONSTANTS::AtomicMassUnit/(1.+XVal);
	Quantity mu_H = CRONOS_CONSTANTS::AtomicMassUnit/XVal;

	mu_num = norm.get_num(mu);
	mu_H_num = norm.get_num(mu_H);
//	cout << " normalised molecular weight (total and hydrogen) " << mu_num << " " << mu_H_num << endl;

	// Compute normalised Boltzmann constant
	kB_num = norm.get_num(CRONOS_CONSTANTS::BoltzmannConstant);
//	cout << " Boltzmann constant " << CRONOS_CONSTANTS::BoltzmannConstant << " " << kB_num << endl;

	// Compute factor preceding cooling time
	Quantity coolingTime_factor_phys = CRONOS_CONSTANTS::BoltzmannConstant*mu_H*mu_H/mu / (gamma-1.);
	coolingTime_factor_num = norm.get_num(coolingTime_factor_phys);

//	cout << endl << " cooling-time factor ";
//	cout << coolingTime_factor_phys << " " << coolingTime_factor_num << endl << endl;
//	rho_0 = norm.get_phys(norm.DENS,1.);
//	m_0 = norm.get_phys(norm.MASS,1.);
//	v_0 = norm.get_phys(norm.VEL, 1.);
//	L_0 = norm.get_phys(norm.LEN, 1.);

}

CoolingFunction::~CoolingFunction() {

}




int CoolingFunction::load_CoolingTable_SchureEtAl2009(normalisation &norm) {

	coolTabSchure_logTemp.resize(111);
	coolTabSchure_Temp_num.resize(111);
	coolTabSchure_logLambda.resize(111);

	ifstream coolingFile;
//	coolingFile.open("./cooling_curve_SchureEtAl2009.txt");
//	coolingFile.open("/home/kissmrbu/src/CronosGit/CronosCode.master-beta/cronos/Util/cooling_curve_SchureEtAl2009.txt");

	string path = CronosRootDir() + "/cronos/Util/cooling_curve_SchureEtAl2009.txt";

	coolingFile.open(path);

	if(coolingFile.fail()) {
		cerr << " cooling table (Schure et al 2009) file does not exist " << endl;
		return -1;
	}  else {
		for (int iEntry=0; iEntry <111; iEntry++) {
			coolingFile >> coolTabSchure_logTemp(iEntry);
			coolingFile >> coolTabSchure_logLambda(iEntry);
			Quantity Tloc_phys = pow(10.,coolTabSchure_logTemp(iEntry))*Kelvin;
			double TLoc_num = norm.get_num(Tloc_phys);
			coolTabSchure_Temp_num(iEntry) = TLoc_num;

		}
	}

//		cout << " Teste 2 norm " << norm.get_unitNormConst("rho") << endl;
//		cout << " Teste ref " << ref_norm.get_unitNormConst("rho") << endl << endl << endl;
//		exit(3);

	return 0;

}


int CoolingFunction::load_CoolingTable_DalgarnoMcCray1972(normalisation &norm) {

	coolTabDGMcC_logTemp.resize(76);
	coolTabDGMcC_Temp_num.resize(76);
	coolTabDGMcC_logLambda.resize(Index::set(0,0),Index::set(75,3));

	ifstream coolingFile;
//	coolingFile.open("./cooling_curve_DalgarnoMcCray1972.txt");
//	coolingFile.open("/home/kissmrbu/src/CronosGit/CronosCode.master-beta/cronos/Util/cooling_curve_DalgarnoMcCray1972.txt");

	string path = CronosRootDir() + "/cronos/Util/cooling_curve_DalgarnoMcCray1972.txt";
	coolingFile.open(path);


	if(coolingFile.fail()) {
		cerr << " cooling table (Dalgarno & McCray 1972) file does not exist " << endl;
		return -1;
	}  else {
		for (int iEntry=0; iEntry <76; iEntry++) {
			coolingFile >> coolTabDGMcC_logTemp(iEntry);
			coolingFile >> coolTabDGMcC_logLambda(iEntry,0);
			coolingFile >> coolTabDGMcC_logLambda(iEntry,1);
			coolingFile >> coolTabDGMcC_logLambda(iEntry,2);
			coolingFile >> coolTabDGMcC_logLambda(iEntry,3);
			Quantity Tloc_phys = pow(10.,coolTabDGMcC_logTemp(iEntry))*Kelvin;
			double TLoc_num = norm.get_num(Tloc_phys);
			coolTabDGMcC_Temp_num(iEntry) = TLoc_num;
		}
	}

	return 0;

}

void CoolingFunction::get_discreteCoolingTable(NumArray<double> &Temp_num, NumArray<double> &logTemp_arr, NumArray<double> &logLam_arr, int choice_curve) {
	// Currently only implemented for cooling curve by Schure et al (2009)
	if(choice_curve==1) {
		logTemp_arr = coolTabSchure_logTemp;
		logLam_arr = coolTabSchure_logLambda;
		Temp_num = coolTabSchure_Temp_num;

	} else {
		cerr << " No further cooling curves implemented, yet " << endl;
		exit(-33);
	}
}


void CoolingFunction::set_minTemp(Quantity minTemp_phys) {
	this->minTemp_phys = minTemp_phys.normalizeTo(CRONOS_BASE_UNITS::Kelvin);
}


double CoolingFunction::get_coolingTime(double Dens_num, double Temp_num, int ionisation_fraction) {

//	cout << " temp num " << Temp_num << endl;
	double coolingRate_num = get_CoolingRate( Temp_num, ionisation_fraction );

	double t_cool_num = coolingTime_factor_num*Temp_num/( Dens_num*coolingRate_num );

//	double t_cool_numVar = kB_num/(gamma-1.) * sqr(mu_H_num)/mu_num / Dens_num * Temp_num/coolingRate_num;
//	cout << " timescale: " << t_cool_num << " " << t_cool_numVar << endl;


//	cout << " factor et al. " << coolingTime_factor_num << " " << coolingRate_num << endl;
	return t_cool_num;
}


double CoolingFunction::get_EnergyLossRate(double Dens_num, double Temp_num, int ionisation_fraction) {

	double coolingRate_num = get_CoolingRate(Temp_num, ionisation_fraction);
//	double energyLossRate_num = coolingRate_num * sqr(numDens_num);

	// Compute normalised hydrogen number density
	double nH_num = Dens_num/mu_H_num;
	double energyLossRate_num = coolingRate_num * sqr(Dens_num/nH_num);

	return energyLossRate_num;

}

double CoolingFunction::get_TempLossRate(double Dens_num, double Temp_num, int ionisation_fraction) {

	// get numerical cooling rate
	double coolingRate_num = get_CoolingRate(Temp_num, ionisation_fraction);
//	double energyLossRate_num = coolingRate_num * sqr(numDens_num);

	double tempLossRate_num = (gamma-1.)*Dens_num * mu_num/sqr(mu_H_num) / kB_num * coolingRate_num;

	return tempLossRate_num;

}


double CoolingFunction::get_CoolingRate(double Temp_num, int ionisation_fraction) {

	// Get physical temperature
//	cout << " phys cooling rate " << endl;
	double Temp_phys = ref_norm.get_phys(ref_norm.TEMP, Temp_num).normalizeTo(CRONOS_BASE_UNITS::Kelvin);
//	cout << " Phys temp: " << Temp_phys << " " << cooling_rate_num << endl;
//	cout << " Temp " << Temp_phys << endl;

	// if temperature is below minimum temperature set by user, no cooling is applied
	if(Temp_phys < minTemp_phys) {
		return 0.;
	}

	// Compute physical cooling rate in SI units [J m^3 s^-1]
	REAL Cooling_rate_phys = 0.;
	// get radiative Cooling by log-linear interpolation of table
//	if ( (Temp_phys > 1.01e4) && (Temp_phys < 1.44e8)) {
	if ( (Temp_phys > 6.309e3) && (Temp_phys < 1.44e8)) {
		REAL logTtest = log10(Temp_phys);
		int iEntry = (logTtest-3.80)*25;
		REAL grad = (coolTabSchure_logLambda(iEntry+1)-coolTabSchure_logLambda(iEntry))/(coolTabSchure_logTemp(iEntry+1)-coolTabSchure_logTemp(iEntry));
		Cooling_rate_phys = coolTabSchure_logLambda(iEntry) + grad*(logTtest-coolTabSchure_logTemp(iEntry));
		Cooling_rate_phys = pow(10.,Cooling_rate_phys); //SI
		//or by extrapolation at flat end
	} else if (Temp_phys >= 1.44e8) {
		Cooling_rate_phys = pow(10,-35.3928);
	} else {
		if(ionisation_fraction==0) {
			Cooling_rate_phys = 0.;
		} else {
			int ion_index(0);
			if(ionisation_fraction==1) {
				ion_index = 3;
			} else if (ionisation_fraction==2) {
				ion_index = 2;
			} else if (ionisation_fraction==3) {
				ion_index = 1;
			} else if (ionisation_fraction==4) {
				ion_index = 0;
			} else {
				cerr << " No such ionisation fraction available " << endl;
				exit(3);
			}

			REAL logTtest = log10(Temp_phys);
			int iEntry = (logTtest-3.80)*25;

			// FindIntervall?
			iEntry = max(iEntry, coolTabDGMcC_logTemp.getLength() - 2);
			REAL grad = (coolTabDGMcC_logLambda(iEntry+1,ion_index)-coolTabDGMcC_logLambda(iEntry,ion_index))/(coolTabDGMcC_logTemp(iEntry+1)-coolTabDGMcC_logTemp(iEntry));
			Cooling_rate_phys = coolTabDGMcC_logLambda(iEntry,ion_index) + grad*(logTtest-coolTabDGMcC_logTemp(iEntry));
			Cooling_rate_phys = pow(10.,Cooling_rate_phys); //SI

		}
	}

//	cout << " cooling " << Cooling_rate_phys << endl;
	double Cooling_rate_num = Cooling_rate_phys * cooling_rate_factor_phys2norm;
//	cout << " cooling (norm) " << Cooling_rate_num << endl;

//	double CoolingKlaus = Cooling_rate_phys* sqr(rho_0/m_0)/rho_0
//					/cube(v_0)*L_0;
//	cout << " nach Klaus: " << CoolingKlaus << endl;



//	Cooling *= sqr(gdata.om[q_rho](ix,jy,kz)*rho_norm/m_H)/rho_norm
//			/cube(c_sound)*RSun;



	//}

	return Cooling_rate_num;


	return 1.;
}
