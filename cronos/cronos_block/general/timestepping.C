#include "timestepping.H"

TimeIntegrator::TimeIntegrator(const int n_saves):
	n_saves(n_saves)
{
	om_save = new NumMatrix<REAL, DIM>[n_saves];
}


TimeIntegrator::~TimeIntegrator()
{
	delete [] om_save;
}


void TimeIntegrator::set_corrField(int qch)
{
	this->qch = qch;

}

void TimeIntegrator::set_IntRange(int ibeg[], int iend[])
{
	for(int d=0; d<DIM; ++d) {
		this->ibeg[d] = ibeg[d];
		this->iend[d] = iend[d];
	}

	for (int i=0; i<n_saves; i++) {
		om_save[i].resize(ibeg, iend);
	}
}


void TimeIntegrator::save_data(const Pot om[], const int substep)
{
	for (int k = ibeg[2]; k <= iend[2]; ++k){
		for (int j = ibeg[1]; j <= iend[1]; ++j){
			for (int i = ibeg[0]; i <= iend[0]; ++i){
				om_save[substep](i,j,k) = om[qch](i,j,k);
			}
		}
	}
}




REAL TimeIntegrator::get_dt(Data & gdata, const int substep)
{
	if (substep == TIME_SUBSTEPS - 1)	return gdata.dt;
	else return 0.;
}





RKSteps::RKSteps():
		TimeIntegrator(1)
{
	this->oneThird  = 1./3.;
	this->twoThirds = 2./3.;
}


REAL RKSteps::get_dt(Data & gdata, const int substep)
{
#if (RK_STEPS == 3)
	if 		(substep == 0) 	return gdata.dt;
	else if (substep == 1) 	return -0.5*gdata.dt;
	else if (substep == 2) 	return 0.5*gdata.dt;
#else
	if 		(substep == 0)	return gdata.dt;
#endif
	else return 0.;
}


void RKSteps::Substep(const Data &gdata, ProblemType &Problem,
                      NumMatrix<REAL, 3> &nom,
                      Pot om[], const int substep) {

#if  (RK_STEPS == 3)
	if(substep == 0) {
		save_data(om, 0);

		for (int k = ibeg[2]; k <= iend[2]; ++k){
			for (int j = ibeg[1]; j <= iend[1]; ++j){
				for (int i = ibeg[0]; i <= iend[0]; ++i){
					om[qch](i,j,k) -= gdata.dt*nom(i,j,k);
				}
			}
		}
	} else if (substep == 1) {
		for (int k = ibeg[2]; k <= iend[2]; ++k) {
			for (int j = ibeg[1]; j <= iend[1]; ++j) {
				for (int i = ibeg[0]; i <= iend[0]; ++i) {
//					om[qch](i,j,k) = 0.75*om[qold](i,j,k)
					om[qch](i,j,k) = 0.75*om_save[0](i,j,k)
						+0.25*om[qch](i,j,k)
						-0.25*gdata.dt*nom(i,j,k);
				}
			}
		}
	} else {
		REAL twth   = 2./3.;
		REAL twthdt = twth*gdata.dt;
		
		for (int k = ibeg[2]; k <= iend[2]; ++k){
			for (int j = ibeg[1]; j <= iend[1]; ++j){
				for (int i = ibeg[0]; i <= iend[0]; ++i){
//					om[qch](i,j,k) = 1./3.*om[qold](i,j,k)
					om[qch](i,j,k) = 1./3.*om_save[0](i,j,k)
						+twth*om[qch](i,j,k)
						-twthdt*nom(i,j,k);
				}
			}
		}
	}
#else

//#if (FLUID_TYPE != CRONOS_MHD)
#if (FLUID_TYPE == CRONOS_HYDRO)
	//	if(Problem.mag == 0) {
	if(om[qch].getName() == "B_x" ||
	   om[qch].getName() == "B_y" ||
	   om[qch].getName() == "B_z" ||
	   om[qch].getName() == "A_x" ||
	   om[qch].getName() == "A_y" ||
	   om[qch].getName() == "A_z") {
		return;
	}
#endif
	//	}

	if (substep == 0) { // First Runge Kutta step
		save_data(om, 0);

		for (int k = ibeg[2]; k <= iend[2]; ++k){
			for (int j = ibeg[1]; j <= iend[1]; ++j){
				for (int i = ibeg[0]; i <= iend[0]; ++i){
					om[qch](i,j,k) -= gdata.dt*nom(i,j,k);
				}
			}
		}
	} else if (substep == 1) { // Second Runge Kutta step
		for (int k = ibeg[2]; k <= iend[2]; ++k) {
			for (int j = ibeg[1]; j <= iend[1]; ++j) {
				for (int i = ibeg[0]; i <= iend[0]; ++i) {
					om[qch](i,j,k) =
						 0.5*om_save[0](i,j,k)
//						 0.5*om[qold](i,j,k)
						+0.5*om[qch](i,j,k)
						-0.5*gdata.dt*nom(i,j,k);
				}
			}
		}
	}
#endif
}



VanLeerIntegrator::VanLeerIntegrator():
		TimeIntegrator(1)
{
}


REAL VanLeerIntegrator::get_dt(Data & gdata, const int substep)
{
	if 		(substep == 0)	return 0.5 * gdata.dt;
	else if	(substep == 1)	return 0.5 * gdata.dt;
	else return 0.;

}

void VanLeerIntegrator::Substep(const Data &gdata, ProblemType &Problem,
                      NumMatrix<REAL, 3> &nom,
                      Pot om[], const int substep) {
	if (substep == 0) {
		save_data(om, 0);
		for (int k = ibeg[2]; k <= iend[2]; ++k){
			for (int j = ibeg[1]; j <= iend[1]; ++j){
				for (int i = ibeg[0]; i <= iend[0]; ++i){
					om[qch](i,j,k) -= 0.5*gdata.dt*nom(i,j,k);
				}
			}
		}

	} else if (substep == 1) {
		for (int k = ibeg[2]; k <= iend[2]; ++k) {
			for (int j = ibeg[1]; j <= iend[1]; ++j) {
				for (int i = ibeg[0]; i <= iend[0]; ++i) {
					om[qch](i,j,k) =
//							om[qold](i,j,k)
							om_save[0](i,j,k)
							- gdata.dt*nom(i,j,k);
				}
			}
		}
	}
}
