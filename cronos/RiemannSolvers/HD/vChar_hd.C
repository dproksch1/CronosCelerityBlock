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

void RiemannSolverHD::get_vChar(Queue queue, Data &gdata,
                                ProblemType &Problem,
                                cronos::vector<double> &iPos,
                                phys_fields_1D &pfL,
                                phys_fields_1D &pfR,
                                NumMatrix<REAL,1> &v_ch_mORIG,
                                NumMatrix<REAL,1> &v_ch_pORIG,
								Buffer<REAL, 1> v_ch_mSYCL,
								Buffer<REAL, 1> v_ch_pSYCL,
//                                fields_1D &fields,
                                const int &dir, REAL &cfl_lin) const {

	//assert_sycl_eq(pfL.uConSYCL[q_rho], pfL.uConORIG[q_rho]);
	//assert_sycl_eq(pfR.uConSYCL[q_rho], pfR.uConORIG[q_rho]);
	//assert_sycl_eq(pfL.uPriSYCL[dir + q_sx], pfL.uPriORIG[dir + q_sx]);
	//assert_sycl_eq(pfR.uPriSYCL[dir + q_sx], pfR.uPriORIG[dir + q_sx]);
	//assert_sycl_eq(pfL.pthermSYCL, pfL.pthermORIG);
	//assert_sycl_eq(pfR.pthermSYCL, pfR.pthermORIG);
	//assert_sycl_eq(v_ch_pSYCL, v_ch_pORIG);
	//assert_sycl_eq(v_ch_mSYCL, v_ch_mORIG);
  
	for (int i = -1; i <= gdata.mx[dir]+1; ++i){

		iPos.set(dir, i-0.5);
    
		REAL rhoinv_p = 1./pfL.uConORIG[q_rho](i  );
		REAL rhoinv_m = 1./pfR.uConORIG[q_rho](i-1);
    
		// Flow velocity
		REAL u_p = pfL.uPriORIG[dir+q_sx](i  );
		REAL u_m = pfR.uPriORIG[dir+q_sx](i-1);
    
		REAL pres_p = pfL.pthermORIG(i);
		REAL pres_m = pfR.pthermORIG(i-1);
    
		// Sound speed
		REAL cs_p   = sqrt(Problem.gamma*pres_p*rhoinv_p);
		REAL cs_m   = sqrt(Problem.gamma*pres_m*rhoinv_m);

		v_ch_pORIG(i) = std::max(std::max(cs_p+u_p,cs_m+u_m),0.);
		v_ch_mORIG(i) = std::max(std::max(cs_p-u_p,cs_m-u_m),0.);

//		if(i==399) {
//			cout << " vfp " << cs_p << " " << cs_m;
//			cout << endl;
//			cout << q_sx << endl;
//		}
    
		REAL vmax = std::max(v_ch_pORIG(i), v_ch_mORIG(i));

//		fields.v_ch_p(i) = std::max(std::max(cs_p+u_p,cs_m+u_m),0.);
//		fields.v_ch_m(i) = std::max(std::max(cs_p-u_p,cs_m-u_m),0.);
//
//		REAL vmax = std::max(fields.v_ch_p(i),fields.v_ch_m(i));

		// Local computation of cfl number
#if  (NON_LINEAR_GRID == CRONOS_OFF)
		REAL cfl_loc = vmax*gdata.idx[dir];
#else
		REAL cfl_loc = vmax*gdata.getCen_idx(dir, i);
#endif

#ifdef GEOM
#if GEOM != CARTESIAN
		if(dir==0) {
			cfl_loc /= gdata.getCen_h0(iPos.get(0), iPos.get(1), iPos.get(2));
		} else if (dir==1) {
			cfl_loc /= gdata.getCen_h1(iPos.get(0), iPos.get(1), iPos.get(2));
		} else {
			cfl_loc /= gdata.getCen_h2(iPos.get(0), iPos.get(1), iPos.get(2));
		}
#endif
#endif
		cfl_lin = std::max(cfl_lin, cfl_loc);
    
	}

	// TODO PHILGS: not setting ipos because it's not used in the compilation path GEOM == CARTESIAN

#ifdef USE_SYCL
	size_t size = gdata.mx[dir];
	REAL gamma = Problem.gamma;
	//int idx = gdata.idx[dir];

	Buffer<REAL, 1> cfl_loc_temp = Buffer<REAL, 1>(Range<1>(size + 2));

	// TODO PHILGS: This requires a neighborhood mapper in Celerity
	queue.submit([&](sycl::handler& cgh) {
		auto uConL_acc = pfL.uConSYCL[q_rho].get_access<cl::sycl::access::mode::read>(cgh);
		auto uConR_acc = pfR.uConSYCL[q_rho].get_access<cl::sycl::access::mode::read>(cgh);
		auto uPriL_acc = pfL.uPriSYCL[dir + q_sx].get_access<cl::sycl::access::mode::read>(cgh);
		auto uPriR_acc = pfR.uPriSYCL[dir + q_sx].get_access<cl::sycl::access::mode::read>(cgh);
		auto pthermL_acc = pfL.pthermSYCL.get_access<cl::sycl::access::mode::read>(cgh);
		auto pthermR_acc = pfR.pthermSYCL.get_access<cl::sycl::access::mode::read>(cgh);
		auto v_ch_p_acc = v_ch_pSYCL.get_access<cl::sycl::access::mode::read_write>(cgh);
		auto v_ch_m_acc = v_ch_mSYCL.get_access<cl::sycl::access::mode::read_write>(cgh);
		auto cfl_loc_temp_acc = cfl_loc_temp.get_access<cl::sycl::access::mode::write>(cgh);
		cgh.parallel_for<class CharacteristicVelocities>(Range<1>(size + 3), Id<1>(1), [=](Item<1> item) {
			size_t i = item.get(0);

			REAL rhoinv_p = 1.0 / uConL_acc[i];
			REAL rhoinv_m = 1.0 / uConR_acc[i-1];

			// Flow velocity
			REAL u_p = uPriL_acc[i];
			REAL u_m = uPriR_acc[i - 1];

			REAL pres_p = pthermL_acc[i];
			REAL pres_m = pthermR_acc[i - 1];

			// Sound speed
			REAL cs_p = sqrt(gamma * pres_p * rhoinv_p);
			REAL cs_m = sqrt(gamma * pres_m * rhoinv_m);

			v_ch_p_acc[i] = std::max(std::max(cs_p + u_p, cs_m + u_m), 0.0);
			v_ch_m_acc[i] = std::max(std::max(cs_p - u_p, cs_m - u_m), 0.0);

			REAL vmax = std::max(v_ch_p_acc[i], v_ch_m_acc[i]);

			// Local computation of cfl number
			REAL cfl_loc = vmax * dir;

			// TODO: This is a reduction and cfl_lin is then used to determine dt, need to fix this somehow
			// save cfl_loc to buffer and aggregate host-side somehow?
			//cfl_lin = std::max(cfl_lin, cfl_loc);
			cfl_loc_temp_acc[i] = cfl_loc;
		});
	});

	//assert_sycl_eq(v_ch_pSYCL, v_ch_pORIG);
	//assert_sycl_eq(v_ch_mSYCL, v_ch_mORIG);

	// This requires a host task in Celerity until reductions are available
	auto cfl_loc_temp_acc_host = cfl_loc_temp.get_access<cl::sycl::access::mode::read>();
	for (size_t i = 1; i < size + 3; ++i) {
		cfl_lin = std::max(cfl_lin, cfl_loc_temp_acc_host[i]);
	}

#endif // USE_SYCL

}



void RiemannSolverHD::get_vChar(const Data &gdata, const ProblemType &Problem,
		const phys_fields_0D &pfL, const phys_fields_0D &pfR, num_fields_0D &f_num,
		int /*ix*/, int /*iy*/, int /*iz*/, int dir, REAL &cfl_lin) const {
	//! Compute characteristic velocities

	int shift_vec[3] = {0,0,0};
	shift_vec[dir] = -1;

	//int iPos[3] = {ix, iy, iz};

	REAL rhoinv_p = 1./pfL.uCon[q_rho];
	REAL rhoinv_m = 1./pfR.uCon[q_rho];

	// Flow velocity
	REAL u_p = pfL.uPri[dir+q_sx];
	REAL u_m = pfR.uPri[dir+q_sx];

	REAL pres_p = pfL.ptherm;
	REAL pres_m = pfR.ptherm;

	// Sound speed
	REAL cs_p   = sqrt(Problem.gamma*pres_p*rhoinv_p);
	REAL cs_m   = sqrt(Problem.gamma*pres_m*rhoinv_m);

	REAL v_ch_p = std::max(std::max(cs_p+u_p,cs_m+u_m),0.);
	REAL v_ch_m = std::max(std::max(cs_p-u_p,cs_m-u_m),0.);

	f_num.v_ch_p = v_ch_p;
	f_num.v_ch_m = v_ch_m;

	REAL vmax = std::max(v_ch_p, v_ch_m);

	//		fields.v_ch_p(i) = std::max(std::max(cs_p+u_p,cs_m+u_m),0.);
	//		fields.v_ch_m(i) = std::max(std::max(cs_p-u_p,cs_m-u_m),0.);
	//
	//		REAL vmax = std::max(fields.v_ch_p(i),fields.v_ch_m(i));

	// Local computation of cfl number
#if  (NON_LINEAR_GRID == CRONOS_OFF)
	REAL cfl_loc = vmax*gdata.idx[dir];
#else
	REAL cfl_loc = vmax*gdata.getCen_idx(dir, iPos[dir]);
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
		const phys_fields_0D &pfL, const phys_fields_0D &pfR, num_fields_0D &f_num, int dir, REAL &cfl_lin) {
	//! Compute characteristic velocities

	int shift_vec[3] = {0,0,0};
	shift_vec[dir] = -1;

	//int iPos[3] = {ix, iy, iz};
	int q_rho = gdata.fluid.get_q_rho();
	int q_sx = gdata.fluid.get_q_sx();

	REAL rhoinv_p = 1./pfL.uCon[q_rho];
	REAL rhoinv_m = 1./pfR.uCon[q_rho];

	// Flow velocity
	REAL u_p = pfL.uPri[dir+q_sx];
	REAL u_m = pfR.uPri[dir+q_sx];

	REAL pres_p = pfL.ptherm;
	REAL pres_m = pfR.ptherm;

	// Sound speed
	REAL cs_p   = sqrt(Problem.gamma*pres_p*rhoinv_p);
	REAL cs_m   = sqrt(Problem.gamma*pres_m*rhoinv_m);

	REAL v_ch_p = std::max(std::max(cs_p+u_p,cs_m+u_m),0.);
	REAL v_ch_m = std::max(std::max(cs_p-u_p,cs_m-u_m),0.);

	f_num.v_ch_p = v_ch_p;
	f_num.v_ch_m = v_ch_m;

	REAL vmax = std::max(v_ch_p, v_ch_m);

	//		fields.v_ch_p(i) = std::max(std::max(cs_p+u_p,cs_m+u_m),0.);
	//		fields.v_ch_m(i) = std::max(std::max(cs_p-u_p,cs_m-u_m),0.);
	//
	//		REAL vmax = std::max(fields.v_ch_p(i),fields.v_ch_m(i));

	// Local computation of cfl number
#if  (NON_LINEAR_GRID == CRONOS_OFF)
	REAL cfl_loc = vmax*gdata.idx[dir];
#else
	REAL cfl_loc = vmax*gdata.getCen_idx(dir, iPos[dir]);
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

