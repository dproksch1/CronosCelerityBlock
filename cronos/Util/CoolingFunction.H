#ifndef COOLING_FUNCTION_H
#define COOLING_FUNCTION_H

#include "normalisation.H"
//#include "units.H" // contained in normalisation.H
#include "problem.H"
#include "matrix.H"

class CoolingFunction {
public:
	CoolingFunction(normalisation &norm, ProblemType &problem, int choice_curve=1);
	~CoolingFunction();
	double get_EnergyLossRate(double numDens_num, double Temp_num, int ionisation_fraction=0);

	double get_TempLossRate(double numDens_num, double Temp_num, int ionisation_fraction=0);

	/**!
	 * Compute cooling rate (currently from Schure et al (2009), only)
	 * -> values multiplied by electron to proton number density => no
	 * compensation for temperature-dependent electron number density is
	 * necessary, here
	 * \return cooling rate Lambda_hd in normalised units
	 * */
	double get_CoolingRate(double Temp_num, int ionisation_fraction=0);

	/** Set a minimum temperature (in physical units) below which no cooling is applied
	 * */
	void set_minTemp(Quantity minTemp_phys);

	double get_coolingTime(double numDens_num, double Temp_num, int ionisation_fraction=0);

	void get_discreteCoolingTable(NumArray<double> &Temp_num, NumArray<double> &logTemp_arr,
			NumArray<double> &logLam_arr, int choice_curve);
private:
	/** Load tabulate cooling-curve data
	 * -> extracted from Schure et al. (2009) A&A 508, 751
	 * Careful: currently only temperatures above 10^4 K are included
	 * */
	int load_CoolingTable_SchureEtAl2009(normalisation &norm);

	/** Load tabulate cooling-curve data relating to Dalgarno & McCray (1972) ARA&A 10, 375
	 * -> extracted from Schure et al. (2009) A&A 508, 751
	 * */
	int load_CoolingTable_DalgarnoMcCray1972(normalisation &norm);

	/** Reference to normalisation class
	 * */
	normalisation &ref_norm;
	NumArray<double> coolTabSchure_logTemp, coolTabSchure_logLambda, coolTabSchure_Temp_num;
	NumArray<double> coolTabDGMcC_logTemp, coolTabDGMcC_Temp_num;
	NumMatrix<double,2> coolTabDGMcC_logLambda;
	double minTemp_phys;
	double cooling_rate_factor_phys2norm, coolingTime_factor_num;
	/** Normalised value of effective molecular weight (total and hydrogen)
	 * (see after Eq. (5) in Townsend (2009))
	 * */
	double mu_num, mu_H_num;
	double kB_num;
	/**Adiabatic exponent
	 * */
	double gamma;
	double rho_0, m_0, v_0, L_0;
};


#endif
