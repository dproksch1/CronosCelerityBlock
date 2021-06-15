#include "CoolingSolver.H"

using namespace CRONOS_BASE_UNITS;


CoolingSolver::CoolingSolver(Data &gdata) {
	q_rho = gdata.fluid.get_q_rho();
	q_Eges = gdata.fluid.get_q_Eges();
	q_wind1 = N_OMINT_USER-2;
	q_wind2 = N_OMINT_USER-1;
	solver_type = 1;
	ionisation_index = 2;
	fCool = NULL;
	T_ref = 0.;
	T_min_num = 0.;

	// Use the exact solver according to Townsend (2009)
	solver_type = 2;
}

CoolingSolver::~CoolingSolver() {
	delete fCool;
}


void CoolingSolver::init(Data &gdata, gridFunc &gfunc, ProblemType &problem) {

	int choice_coolingCurve = 1;
	fCool = new CoolingFunction(*problem.TrafoNorm, problem, choice_coolingCurve);
	prepare_exactSolver();

	// Check for lower temperature boundary
	if(problem.force_min(q_Eges)) {
		T_min_num = problem.min_Val(q_Eges);
	} else {
		T_min_num = 0.;
	}

	// Set minimum also for cooling rate
	Quantity T_min_phys = problem.TrafoNorm->get_phys(problem.TrafoNorm->TEMP, T_min_num);
	if(gdata.rank==0) {
		cout << " Setting minimum Temperature " << T_min_phys << " --> " << T_min_num << endl;
	}
	fCool->set_minTemp(T_min_phys);

	bool do_test = false;

	if(do_test) {
		test_me();
	}

}


REAL CoolingSolver::singlestep(Data &gdata, gridFunc &gfunc, ProblemType &problem) {

//	// Test for cooling solver
//	gdata.om[q_Eges](12,12,12) = 1.e3;
//	gdata.om[0](12,12,12) = 1.;
	if(gdata.rank==0) {
		cout << " Cooling step with dt " << gdata.dt << endl;
	}


	if(solver_type==1) {
		return singlestep_adaptiveEuler(gdata, gfunc, problem);
	} else {
		return singlestep_exactIntScheme(gdata, gfunc, problem);
	}

}

int CoolingSolver::test_me() {

	REAL rho_num = 1.;
	REAL Temp_num = 1.43e4; //  [K]
	Temp_num = 1.e3;

	REAL time = 0.;
	REAL dt = 1.e-5;


	double Temp_num_new = get_Tnext_exactIntScheme(rho_num, Temp_num, dt);

	double Temp_num_Euler = get_Tnext_adaptiveEuler(rho_num, Temp_num, dt);


	cout << " Test " << Temp_num << " " << Temp_num_new << " " << Temp_num_Euler << endl << endl << endl;

	double YVal = get_Y(Temp_num);
	double Temp_var = get_Yinv(YVal);

	cout << " Hin " << YVal << " und zurueck " << Temp_var << " vs " << Temp_num << endl;

	exit(3);

	return 0;
}


REAL CoolingSolver::singlestep_adaptiveEuler(Data &gdata, gridFunc &gfunc, ProblemType &problem) {

	if(gdata.om[q_Eges].getName() != "Temp") {
		cerr << " Currently cooling is implemented for the energy, only " << endl;
		exit(-54);
	}

	for(int k=0; k<=gdata.mx[2]; ++k) {
		for(int j=0; j<=gdata.mx[1]; ++j) {
			for(int i=0; i<=gdata.mx[0]; ++i) {

				double rho_num = gdata.om[q_rho](i,j,k);

				REAL time = gdata.time;
				REAL dt = gdata.dt;

				while (time < gdata.time+gdata.dt){
					REAL Temp_num = gdata.om[q_Eges](i,j,k); //  [K]

					if(Temp_num < T_min_num) {
						break;
						Temp_num = T_min_num;
					}

					// get initial temperature-loss rate
					double TLossRate = fCool->get_TempLossRate(rho_num, Temp_num, ionisation_index);
					double CoolingTime = fCool->get_coolingTime(rho_num, Temp_num, ionisation_index);
//					cout << " numbers " << dt << " " << CoolingTime << " " << TLossRate << endl;
//					cout << " Vals " << Temp_num << " " << rho_num << endl;

					// if time step is too large, need to reduce time step size
					if(dt/CoolingTime > 0.00001) {
						dt *= 0.5;
					} else {
						gdata.om[q_Eges](i,j,k) -= TLossRate*dt;
						time += dt;
					}
//					cout << " Cooling down " << time << " " << gdata.om[q_Eges](i,j,k) << endl;

					if(Temp_num < T_min_num) {
						Temp_num = T_min_num;
						break;
					}

				}



			}
		}
	}

	return 0.;
}


REAL CoolingSolver::get_Tnext_adaptiveEuler(double rho_num, double T_old, double del_t) {
	double dt_full = del_t;
	double dt_done;
	REAL Temp_num = T_old; //  [K]

	while(dt_done < dt_full) {
		// get initial temperature-loss rate
		double TLossRate = fCool->get_TempLossRate(rho_num, Temp_num, ionisation_index);
		double CoolingTime = fCool->get_coolingTime(rho_num, Temp_num, ionisation_index);
//		cout << " numbers " << del_t << " " << CoolingTime << " " << TLossRate << endl;
//		cout << " Vals " << Temp_num << " " << rho_num << endl;

		// if time step is too large, need to reduce time step size
		if(del_t/CoolingTime > 0.00001) {
			del_t *= 0.5;
		} else {
			Temp_num -= TLossRate*del_t;
			dt_done += del_t;
		}
//		cout << " Cooling down " << time << " " << Temp_num << endl;
	}
	return Temp_num;
}


REAL CoolingSolver::singlestep_exactIntScheme(Data &gdata, gridFunc &gfunc, ProblemType &problem) {

	if(gdata.om[q_Eges].getName() != "Temp") {
		cerr << " Currently cooling is implemented for the energy, only " << endl;
		exit(-54);
	}

	for(int k=0; k<=gdata.mx[2]; ++k) {
		for(int j=0; j<=gdata.mx[1]; ++j) {
			for(int i=0; i<=gdata.mx[0]; ++i) {

				REAL rho_num = gdata.om[q_rho](i,j,k);
				REAL Temp_num = gdata.om[q_Eges](i,j,k); //  [K]

				REAL time = gdata.time;
				REAL dt = gdata.dt;


				gdata.om[q_Eges](i,j,k) = get_Tnext_exactIntScheme(rho_num, Temp_num, dt);

			}
		}
	}

	return 0.;
}


REAL CoolingSolver::get_Tnext_exactIntScheme(REAL rho_num, REAL T_old_num, REAL del_t) {

	if(T_old_num < T_min_num) {
		return T_min_num;
	}

	// Implementation of Eq. (26) in Townsend (2009)
	double yVal = get_Y(T_old_num);
	double yMod = T_old_num/T_ref*fCool->get_CoolingRate(T_ref)/fCool->get_CoolingRate(T_old_num)*del_t/fCool->get_coolingTime(rho_num,T_old_num, ionisation_index);
	yVal += yMod;
//	cout << " YVal problem ";
//	cout << get_Y(T_old_num) << " " << T_ref << " " << T_old_num << " " << endl;
	double yOld = get_Y(T_old_num);
	double lamOld = fCool->get_CoolingRate(T_old_num);
	double lamRef = fCool->get_CoolingRate(T_ref);
	double tCool = fCool->get_coolingTime(rho_num,T_old_num, ionisation_index);
//	cout << " YVal ";
//	cout << yOld << " ";
//	cout << T_old_num << " " << T_ref << " " << lamRef << " " << lamOld << " ";
//	cout << tCool << endl;
//	cout << " Cooling time " << tCool << " " << del_t << endl;
//	cout << " Rate: " << lamOld << " " << lamRef << endl;
//	cout << " Y " << yVal-1. << " " << yOld-1. << " " << yMod << endl;
//	cout << " Mod: ";
//	cout << T_old_num/T_ref << " ";
//	cout << lamRef/lamOld << " ";
//	cout << del_t/tCool << " ";
//	cout << T_old_num << " " << T_ref << " ";
//	cout << endl;
//	exit(3);

	REAL T_new = get_Yinv(yVal);

	if(T_new < T_min_num) {
		T_new = T_min_num;
	}

	return T_new;
}

int CoolingSolver::get_tempInterval(double Temp_num) {
	// Find relevant temperature bin:
	int numPoints = discTemp_num.getLength();
	int i_point;
	for(i_point=0; i_point<numPoints; ++i_point) {
		if(Temp_num<discTemp_num(i_point)) {
			break;
		}
	}

	i_point -= 1;

	// Check if okay
//	cout << " Temps " << discTemp_num(i_point) << " < " << Temp_num << " < " << discTemp_num(i_point+1) << endl;

	return i_point;
}

int CoolingSolver::find_Interval(NumArray<double> &val_arr, double value) {
	// Find relevant temperature bin:
	int numPoints = val_arr.getLength();
//	cout << " Trying to iterate " << numPoints << endl;
	int i_point;
	for(i_point=0; i_point<=numPoints-2; ++i_point) {
//		cout << " iter " << i_point << " " << value << " " << val_arr(i_point) << " " << val_arr(i_point+1) << " " << (value-val_arr(i_point))*(value-val_arr(i_point+1)) << endl;
//		if(value<val_arr(i_point)) {
		if(((value-val_arr(i_point))*(value-val_arr(i_point+1))<0.) || std::abs(value-val_arr(i_point))<1.e-15) {
			return i_point;
		}
	}

//	i_point -= 1;

	// Check if okay
//	cout << " Interval " << i_point << " " <<  val_arr(i_point) << " < " << value << " < " << val_arr(i_point+1) << endl;

	return numPoints-2;
}


REAL CoolingSolver::get_Y(REAL Temp_num) {

//	double i_point = get_tempInterval(Temp_num);
	// Find temperature interval, for which T_k <= Temp_num <= T_k+1
//	cout << " finding this " << Temp_num << endl;
	int i_point = find_Interval(discTemp_num, Temp_num);

	double YTemp = 0.;
	double alpha_k = powerLawLambda(i_point);

	double lambdaN = discLam_phys(discLam_phys.getLength()-1);
	double lambdaK = discLam_phys(i_point);

	double TempN_num = discTemp_num(discTemp_num.getLength()-1);
	double TempK_num = discTemp_num(i_point);


//	cout << " Both temps " << TempN_num << " " << TempK_num  << endl;
//	cout << " Both cooling rates " << lambdaN << " " << lambdaK << endl;

	if(std::abs(alpha_k-1.)>1.e-9) {
		YTemp = 1./(1.-alpha_k)*lambdaN/lambdaK*TempK_num/TempN_num*(1.-pow(TempK_num/Temp_num, alpha_k-1.));
	} else {
		YTemp = lambdaN/lambdaK*TempK_num/TempN_num*log(TempK_num/Temp_num);
	}
//	cout << " YTemp " << YTemp+ Yk_bin(i_point) << endl;

	return YTemp + Yk_bin(i_point);

}

REAL CoolingSolver::get_Yinv(double YVal) {
//	cout << " Temp num " << Temp_num << endl;

	// Find Y interval, for which Y_k <= Y <= Y_k+1
//	cout << " Trying to find Y interval " << YVal << endl;
//	cout << Yk_bin << endl;
	double i_point = find_Interval(Yk_bin, YVal);
//	double i_point = get_tempInterval(Temp_num);

	double alpha_k = powerLawLambda(i_point);

	double lambdaN = discLam_phys(discLam_phys.getLength()-1);
	double lambdaK = discLam_phys(i_point);

	double TempN_num = discTemp_num(discTemp_phys.getLength()-1);
	double TempK_num = discTemp_num(i_point);

	double YInv = 0.;

	if(std::abs(alpha_k-1.)>1.e-9) {
		YInv = TempK_num*pow(1.-(1.-alpha_k)*lambdaK/lambdaN * TempN_num/TempK_num*(YVal - Yk_bin(i_point)), 1./(1-alpha_k));
	} else {
		YInv = TempK_num*exp(-lambdaK/lambdaN*TempN_num/TempK_num*(YVal - Yk_bin(i_point)));
	}
	return YInv;

}



void CoolingSolver::prepare_exactSolver() {
	NumArray<double> logTemp, logLambda;
	// get discrete cooling curve by Schure et al (2009)
	fCool->get_discreteCoolingTable(discTemp_num, logTemp, logLambda,1);

	// get number of discrete points of cooling curve
	int numPoints = logTemp.getLength();
	int numBins = numPoints-1;

	// --> point grid from 0 to numPoints-1
	// --> bins from 0 to numPoints-2
	powerLawLambda.resize(numBins);

	// Here, we use SI units, since normalisation is not relevant -> we only consider ratios
	NumArray<double> lambda_phys(numPoints), Temp_phys(numPoints);
	for(int i_point=0; i_point<numPoints; i_point++) {
		lambda_phys(i_point) = pow(10.,logLambda(i_point));
		Temp_phys(i_point) = pow(10.,logTemp(i_point));
	}
	discTemp_phys = Temp_phys;
	discLam_phys = lambda_phys;

	// compute power-law index alpha in each bin
	cout << " Power-law index " << endl;
	for(int i_bin=0; i_bin<numPoints-1; i_bin++) {
		double alpha = (logLambda(i_bin+1) - logLambda(i_bin))/(logTemp(i_bin+1) - logTemp(i_bin));
		powerLawLambda(i_bin) = alpha;
//		cout << i_bin << " " << alpha << endl;
	}

	// compute Yk for each bin according to Eq. (A6) in Townsend (2009) -> given at discrete points of cooling curve
	Yk_bin.resize(numPoints);
	Yk_bin(numPoints-1) = 0.; //  YN = 0 according to Townsend (2009)

	double TempN_phys = Temp_phys(numPoints-1);
	double lamN_phys = lambda_phys(numPoints-1);
	T_ref = discTemp_num(numPoints-1);

	for(int i_point=numPoints-2; i_point>=0; i_point--) {
		double alpha_k = powerLawLambda(i_point);
		Yk_bin(i_point) = 1.;
		double del_Yk = 0.;
		if(std::abs(alpha_k-1.) > 1.e-9) {
			del_Yk = 1./(1.-alpha_k)*lamN_phys/lambda_phys(i_point)*Temp_phys(i_point)/TempN_phys*(1.-pow(Temp_phys(i_point)/Temp_phys(i_point+1), alpha_k-1.));
		} else {
			del_Yk = lamN_phys/lambda_phys(i_point)*Temp_phys(i_point)/TempN_phys*log(Temp_phys(i_point)/Temp_phys(i_point+1));
		}
		Yk_bin(i_point) = Yk_bin(i_point+1) - del_Yk;
		//		cout << " Yk(" << i_point << ")="<< Yk_bin(i_point)-1. << endl;
	}
//	exit(3);

}


