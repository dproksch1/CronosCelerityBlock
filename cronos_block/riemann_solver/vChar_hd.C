#include "RiemannSolverHD.H"

RiemannSolverHD::RiemannSolverHD(const Data &gdata, int dir, int Fluid_Type) : RiemannSolver(gdata, dir, Fluid_Type) {
#if(FLUID_TYPE==CRONOS_HYDRO)
	this->q_rho = gdata.fluid.get_q_rho();
	this->q_sx  = gdata.fluid.get_q_sx();
	this->q_sy  = gdata.fluid.get_q_sy();
	this->q_sz  = gdata.fluid.get_q_sz();
	this->q_Eges = gdata.fluid.get_q_Eges();
	this->q_Eadd = gdata.fluid.get_q_Eadd();
#endif
}

void RiemannSolverHD::get_vChar(const Data &gdata, const ProblemType &Problem,
		const phys_fields_0D &pfL, const phys_fields_0D &pfR, num_fields_0D &f_num,
		int /*ix*/, int /*iy*/, int /*iz*/, int dir, double &cfl_lin) const {
	//! Compute characteristic velocities

	int shift_vec[3] = {0,0,0};
	shift_vec[dir] = -1;

	//int iPos[3] = {ix, iy, iz};

	double rhoinv_p = 1./pfL.uCon[q_rho];
	double rhoinv_m = 1./pfR.uCon[q_rho];

	// Flow velocity
	double u_p = pfL.uPri[dir+q_sx];
	double u_m = pfR.uPri[dir+q_sx];

	double pres_p = pfL.ptherm;
	double pres_m = pfR.ptherm;

	// Sound speed
	double cs_p   = sqrt(Problem.gamma*pres_p*rhoinv_p);
	double cs_m   = sqrt(Problem.gamma*pres_m*rhoinv_m);

	double v_ch_p = std::max(std::max(cs_p+u_p,cs_m+u_m),0.);
	double v_ch_m = std::max(std::max(cs_p-u_p,cs_m-u_m),0.);

	f_num.v_ch_p = v_ch_p;
	f_num.v_ch_m = v_ch_m;

	double vmax = std::max(v_ch_p, v_ch_m);

	//		fields.v_ch_p(i) = std::max(std::max(cs_p+u_p,cs_m+u_m),0.);
	//		fields.v_ch_m(i) = std::max(std::max(cs_p-u_p,cs_m-u_m),0.);
	//
	//		double vmax = std::max(fields.v_ch_p(i),fields.v_ch_m(i));

	// Local computation of cfl number
#if  (NON_LINEAR_GRID == CRONOS_OFF)
	double cfl_loc = vmax*gdata.idx[dir];
#else
	double cfl_loc = vmax*gdata.getCen_idx(dir, iPos[dir]);
#endif

#ifdef GEOM
#if GEOM != CARTESIAN
	if(dir==0) {
		cfl_loc /= gdata.h0(ix, iy, iz, shift_vec[0],shift_vec[1],shift_vec[2]);
	} else if (dir==1) {
		cfl_loc /= gdata.h1(ix, iy, iz, shift_vec[0],shift_vec[1],shift_vec[2]);
	} else {
		cfl_loc /= gdata.h2(ix, iy, iz, shift_vec[0],shift_vec[1],shift_vec[2]);
	}
#endif
#endif
	cfl_lin = std::max(cfl_lin, cfl_loc);

}

void get_vChar2(const Data &gdata, const ProblemType &Problem,
		const phys_fields_0D &pfL, const phys_fields_0D &pfR, num_fields_0D &f_num, int dir, double &cfl_lin) {
	//! Compute characteristic velocities

	int shift_vec[3] = {0,0,0};
	shift_vec[dir] = -1;

	//int iPos[3] = {ix, iy, iz};
	int q_rho = gdata.fluid.get_q_rho();
	int q_sx = gdata.fluid.get_q_sx();

	double rhoinv_p = 1./pfL.uCon[q_rho];
	double rhoinv_m = 1./pfR.uCon[q_rho];

	// Flow velocity
	double u_p = pfL.uPri[dir+q_sx];
	double u_m = pfR.uPri[dir+q_sx];

	double pres_p = pfL.ptherm;
	double pres_m = pfR.ptherm;

	// Sound speed
	double cs_p   = sqrt(Problem.gamma*pres_p*rhoinv_p);
	double cs_m   = sqrt(Problem.gamma*pres_m*rhoinv_m);

	double v_ch_p = std::max(std::max(cs_p+u_p,cs_m+u_m),0.);
	double v_ch_m = std::max(std::max(cs_p-u_p,cs_m-u_m),0.);

	f_num.v_ch_p = v_ch_p;
	f_num.v_ch_m = v_ch_m;

	double vmax = std::max(v_ch_p, v_ch_m);

	//		fields.v_ch_p(i) = std::max(std::max(cs_p+u_p,cs_m+u_m),0.);
	//		fields.v_ch_m(i) = std::max(std::max(cs_p-u_p,cs_m-u_m),0.);
	//
	//		double vmax = std::max(fields.v_ch_p(i),fields.v_ch_m(i));

	// Local computation of cfl number
#if  (NON_LINEAR_GRID == CRONOS_OFF)
	double cfl_loc = vmax*gdata.idx[dir];
#else
	double cfl_loc = vmax*gdata.getCen_idx(dir, iPos[dir]);
#endif

#ifdef GEOM
#if GEOM != CARTESIAN
	if(dir==0) {
		cfl_loc /= gdata.h0(ix, iy, iz, shift_vec[0],shift_vec[1],shift_vec[2]);
	} else if (dir==1) {
		cfl_loc /= gdata.h1(ix, iy, iz, shift_vec[0],shift_vec[1],shift_vec[2]);
	} else {
		cfl_loc /= gdata.h2(ix, iy, iz, shift_vec[0],shift_vec[1],shift_vec[2]);
	}
#endif
#endif
	cfl_lin = std::max(cfl_lin, cfl_loc);

}

double get_vChar3(const Data &gdata, const ProblemType &Problem,
		const phys_fields_0D &pfL, const phys_fields_0D &pfR, num_fields_0D &f_num, int dir) {
	//! Compute characteristic velocities

	int shift_vec[3] = {0,0,0};
	shift_vec[dir] = -1;

	//int iPos[3] = {ix, iy, iz};
	int q_rho = gdata.fluid.get_q_rho();
	int q_sx = gdata.fluid.get_q_sx();

	double rhoinv_p = 1./pfL.uCon[q_rho];
	double rhoinv_m = 1./pfR.uCon[q_rho];

	// Flow velocity
	double u_p = pfL.uPri[dir+q_sx];
	double u_m = pfR.uPri[dir+q_sx];

	double pres_p = pfL.ptherm;
	double pres_m = pfR.ptherm;

	// Sound speed
	double cs_p   = sqrt(Problem.gamma*pres_p*rhoinv_p);
	double cs_m   = sqrt(Problem.gamma*pres_m*rhoinv_m);

	double v_ch_p = std::max(std::max(cs_p+u_p,cs_m+u_m),0.);
	double v_ch_m = std::max(std::max(cs_p-u_p,cs_m-u_m),0.);

	f_num.v_ch_p = v_ch_p;
	f_num.v_ch_m = v_ch_m;

	double vmax = std::max(v_ch_p, v_ch_m);

	//		fields.v_ch_p(i) = std::max(std::max(cs_p+u_p,cs_m+u_m),0.);
	//		fields.v_ch_m(i) = std::max(std::max(cs_p-u_p,cs_m-u_m),0.);
	//
	//		double vmax = std::max(fields.v_ch_p(i),fields.v_ch_m(i));

	// Local computation of cfl number
#if  (NON_LINEAR_GRID == CRONOS_OFF)
	double cfl_loc = vmax*gdata.idx[dir];
#else
	double cfl_loc = vmax*gdata.getCen_idx(dir, iPos[dir]);
#endif

#ifdef GEOM
#if GEOM != CARTESIAN
	if(dir==0) {
		cfl_loc /= gdata.h0(ix, iy, iz, shift_vec[0],shift_vec[1],shift_vec[2]);
	} else if (dir==1) {
		cfl_loc /= gdata.h1(ix, iy, iz, shift_vec[0],shift_vec[1],shift_vec[2]);
	} else {
		cfl_loc /= gdata.h2(ix, iy, iz, shift_vec[0],shift_vec[1],shift_vec[2]);
	}
#endif
#endif
	
	return cfl_loc;
}

