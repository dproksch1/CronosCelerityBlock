#include "RiemannSolverMHD.H"

RiemannSolverMHD::RiemannSolverMHD(const Data &gdata, int dir, int Fluid_Type) : RiemannSolver(gdata, dir, Fluid_Type) {
	norm = 1./sqrt(2.);
#if(FLUID_TYPE==CRONOS_MHD)
	this->q_rho = gdata.fluid.get_q_rho();
	this->q_sx  = gdata.fluid.get_q_sx();
	this->q_sy  = gdata.fluid.get_q_sy();
	this->q_sz  = gdata.fluid.get_q_sz();
	this->q_Bx  = gdata.fluid.get_q_Bx();
	this->q_By  = gdata.fluid.get_q_By();
	this->q_Bz  = gdata.fluid.get_q_Bz();
	this->q_Eges = gdata.fluid.get_q_Eges();
	this->q_Eadd = gdata.fluid.get_q_Eadd();
#endif
}

void RiemannSolverMHD::get_vChar(Queue queue, Data &gdata,
                                 ProblemType &Problem,
                                 cronos::vector<double> &iPos,
                                 phys_fields_1D &pfL,
                                 phys_fields_1D &pfR,
                                 NumMatrix<REAL,1> &v_ch_m, NumMatrix<REAL,1> &v_ch_p,
								 Buffer<REAL, 1> v_ch_mSYCL,
								 Buffer<REAL, 1> v_ch_pSYCL,
//                                 fields_1D &fields,
                                 const int &dir, REAL &cfl_lin) const {

  
	for (int i = -1; i <= gdata.mx[dir]+1; ++i){

		iPos.set(dir, i-0.5);
    
		REAL rhoinv_p = 1./pfL.uConORIG[q_rho](i  );
		REAL rhoinv_m = 1./pfR.uConORIG[q_rho](i-1);
    
		// Flow velocity
		REAL u_p = pfL.uPriORIG[dir+q_sx](i  );
		REAL u_m = pfR.uPriORIG[dir+q_sx](i-1);
    
		REAL bsqr_p = (sqr(pfL.uConORIG[q_Bx](i  )) +
		               sqr(pfL.uConORIG[q_By](i  )) +
		               sqr(pfL.uConORIG[q_Bz](i  )));
		REAL bsqr_m = (sqr(pfR.uConORIG[q_Bx](i-1)) +
		               sqr(pfR.uConORIG[q_By](i-1)) +
		               sqr(pfR.uConORIG[q_Bz](i-1)));

		// Alfven velocity (global)
		REAL va2_p  = bsqr_p*rhoinv_p;
		REAL va2_m  = bsqr_m*rhoinv_m;

		// Alfven velocity (directional)
		REAL val2_p = sqr(pfL.uConORIG[dir+q_Bx](i  ))*rhoinv_p;
		REAL val2_m = sqr(pfR.uConORIG[dir+q_Bx](i-1))*rhoinv_m;

		REAL pres_p = pfL.pthermORIG(i);
		REAL pres_m = pfR.pthermORIG(i-1);
    
		// Sound speed
		REAL c2_p   = Problem.gamma*pres_p*rhoinv_p;
		REAL c2_m   = Problem.gamma*pres_m*rhoinv_m;

		// Fast mode speed 
		REAL vf_p = norm*sqrt(va2_p+c2_p + 
		                      sqrt(sqr(va2_p+c2_p)-4.*c2_p*val2_p));
		REAL vf_m = norm*sqrt(va2_m+c2_m + 
		                      sqrt(sqr(va2_m+c2_m)-4.*c2_m*val2_m));

//		if(i==399) {
//			cout << " vfp " << vf_p << " ";
//			cout << va2_p << " " << c2_p << " ";
//			cout << pfL.uCon[dir+q_Bx](i  ) << " ";
//			cout << pfL.uCon[dir+q_By](i  ) << " ";
//			cout << pfL.uCon[dir+q_Bz](i  ) << " ";
//			cout << endl;
//			cout << q_Bx << " " << q_By << " " << q_Bz << endl;
//			cout << q_sx << endl;
//		}

		v_ch_p(i) = std::max(std::max(vf_p+u_p,vf_m+u_m),0.);
		v_ch_m(i) = std::max(std::max(vf_p-u_p,vf_m-u_m),0.);

		REAL vmax = std::max(v_ch_p(i), v_ch_m(i));

//		fields.v_ch_p(i) = std::max(std::max(vf_p+u_p,vf_m+u_m),0.);
//		fields.v_ch_m(i) = std::max(std::max(vf_p-u_p,vf_m-u_m),0.);
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


}



void RiemannSolverMHD::get_vChar(const Data &gdata, const ProblemType &Problem,
	const phys_fields_0D &pfL, const phys_fields_0D &pfR, num_fields_0D &f_num,
		int ix, int iy, int iz, int dir, REAL &cfl_lin) const {
	//! Compute characteristic velocities

	int shift_vec[3] = {0,0,0};
	shift_vec[dir] = -1;

	int iPos[3] = {ix, iy, iz};

	REAL rhoinv_p = 1./pfL.uCon[q_rho];
	REAL rhoinv_m = 1./pfR.uCon[q_rho];

	// Flow velocity
	REAL u_p = pfL.uPri[dir+q_sx];
	REAL u_m = pfR.uPri[dir+q_sx];

	REAL bsqr_p = (sqr(pfL.uCon[q_Bx]) +
			sqr(pfL.uCon[q_By]) +
			sqr(pfL.uCon[q_Bz]));
	REAL bsqr_m = (sqr(pfR.uCon[q_Bx]) +
			sqr(pfR.uCon[q_By]) +
			sqr(pfR.uCon[q_Bz]));

	// Alfven velocity (global)
	REAL va2_p  = bsqr_p*rhoinv_p;
	REAL va2_m  = bsqr_m*rhoinv_m;

	// Alfven velocity (directional)
	REAL val2_p = sqr(pfL.uCon[dir+q_Bx])*rhoinv_p;
	REAL val2_m = sqr(pfR.uCon[dir+q_Bx])*rhoinv_m;

	REAL pres_p = pfL.ptherm;
	REAL pres_m = pfR.ptherm;

	// Sound speed
	REAL c2_p   = Problem.gamma*pres_p*rhoinv_p;
	REAL c2_m   = Problem.gamma*pres_m*rhoinv_m;

	// Fast mode speed
	REAL vf_p = norm*sqrt(va2_p+c2_p +
			sqrt(sqr(va2_p+c2_p)-4.*c2_p*val2_p));
	REAL vf_m = norm*sqrt(va2_m+c2_m +
			sqrt(sqr(va2_m+c2_m)-4.*c2_m*val2_m));

	REAL v_ch_p = std::max(std::max(vf_p+u_p,vf_m+u_m),0.);
	REAL v_ch_m = std::max(std::max(vf_p-u_p,vf_m-u_m),0.);

	f_num.v_ch_p = v_ch_p;
	f_num.v_ch_m = v_ch_m;

	REAL vmax = std::max(v_ch_p, v_ch_m);

	//		fields.v_ch_p(i) = std::max(std::max(vf_p+u_p,vf_m+u_m),0.);
	//		fields.v_ch_m(i) = std::max(std::max(vf_p-u_p,vf_m-u_m),0.);
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



void RiemannSolverMHD::get_vChar2D(Data &gdata,
                                   ProblemType &Problem,
                                   cronos::vector<REAL> &iPos,
                                   phys_fields_2D &pfLL,
                                   phys_fields_2D &pfLR,
                                   phys_fields_2D &pfRL,
                                   phys_fields_2D &pfRR,
                                   fields_2D &fields,
                                   const int &dir, REAL &cfl_lin) {

	int dir0, dir1;
	assert(dir >= 0 && dir < DIM);
	if(dir == 0) {
		dir0 = 1;
		dir1 = 2;
	} else if (dir == 1) {
		dir0 = 2;
		dir1 = 0;
	} else {
		dir0 = 0;
		dir1 = 1;
	}

	for (int i = -1; i <= gdata.mx[dir0]+1; ++i){
		for (int j = -1; j <= gdata.mx[dir1]+1; ++j){
			fields.v_cor_m[0](i,j) = std::max(fields.v_ch_m[0](i,j  ),
			                                  fields.v_ch_m[0](i,j-1));
			fields.v_cor_p[0](i,j) = std::max(fields.v_ch_p[0](i,j  ),
			                                  fields.v_ch_p[0](i,j-1));
			fields.v_cor_m[1](i,j) = std::max(fields.v_ch_m[1](i  ,j),
			                                  fields.v_ch_m[1](i-1,j));
			fields.v_cor_p[1](i,j) = std::max(fields.v_ch_p[1](i  ,j),
			                                  fields.v_ch_p[1](i-1,j));

		}
	}
}

