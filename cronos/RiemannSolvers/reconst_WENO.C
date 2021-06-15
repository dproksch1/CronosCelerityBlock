#include "reconst_WENO.H"
#include <iostream>
#include <stdlib.h>

using namespace std;

template<typename T>
T power(T base, int expo)
{
  T ret = 1.;
  for (int mul = 1; mul <= abs(expo); mul++)
    ret *= base;
  return (expo > 0) ? ret : 1./ret;
}


SingleReconstruction_WENO::SingleReconstruction_WENO(const Data &gdata, const CronosFluid &fluid, int dir, int qReconst, int substep):
	SingleReconstruction(gdata, fluid, dir, qReconst, substep)
{

	w_l.resize(Index::set(-2),Index::set(gdata.mx[dir]+1));
	w_c.resize(Index::set(-2),Index::set(gdata.mx[dir]+1));
	w_r.resize(Index::set(-2),Index::set(gdata.mx[dir]+1));

	dudx_l.resize(Index::set(-2),Index::set(gdata.mx[dir]+1));
	dudx_c.resize(Index::set(-2),Index::set(gdata.mx[dir]+1));
	dudx_r.resize(Index::set(-2),Index::set(gdata.mx[dir]+1));

	w_lP.resize(Index::set(-2),Index::set(gdata.mx[dir]+1));
	w_cP.resize(Index::set(-2),Index::set(gdata.mx[dir]+1));
	w_rP.resize(Index::set(-2),Index::set(gdata.mx[dir]+1));
	w_lM.resize(Index::set(-2),Index::set(gdata.mx[dir]+1));
	w_cM.resize(Index::set(-2),Index::set(gdata.mx[dir]+1));
	w_rM.resize(Index::set(-2),Index::set(gdata.mx[dir]+1));

	dudx_lP.resize(Index::set(-2),Index::set(gdata.mx[dir]+1));
	dudx_cP.resize(Index::set(-2),Index::set(gdata.mx[dir]+1));
	dudx_rP.resize(Index::set(-2),Index::set(gdata.mx[dir]+1));
	dudx_lM.resize(Index::set(-2),Index::set(gdata.mx[dir]+1));
	dudx_cM.resize(Index::set(-2),Index::set(gdata.mx[dir]+1));
	dudx_rM.resize(Index::set(-2),Index::set(gdata.mx[dir]+1));

	fac_cen = 13./3.;
	eps_WENO = 1.e-12;
	p_WENO = 2;
	use_Limiter = true;

	if (gdata.rank == 0) {
		cout << "  Using Reconst: " << "WENO" << " - " << fluid.get_fieldName(qReconst) << " - dir " << dir << " - substep " << substep << endl;
	}

}


SingleReconstruction_WENO::SingleReconstruction_WENO(const Data &gdata, int dir, int substep):
	SingleReconstruction(gdata, dir, substep)
{

	w_l.resize(Index::set(-2),Index::set(gdata.mx[dir]+1));
	w_c.resize(Index::set(-2),Index::set(gdata.mx[dir]+1));
	w_r.resize(Index::set(-2),Index::set(gdata.mx[dir]+1));

	dudx_l.resize(Index::set(-2),Index::set(gdata.mx[dir]+1));
	dudx_c.resize(Index::set(-2),Index::set(gdata.mx[dir]+1));
	dudx_r.resize(Index::set(-2),Index::set(gdata.mx[dir]+1));

	w_lP.resize(Index::set(-2),Index::set(gdata.mx[dir]+1));
	w_cP.resize(Index::set(-2),Index::set(gdata.mx[dir]+1));
	w_rP.resize(Index::set(-2),Index::set(gdata.mx[dir]+1));
	w_lM.resize(Index::set(-2),Index::set(gdata.mx[dir]+1));
	w_cM.resize(Index::set(-2),Index::set(gdata.mx[dir]+1));
	w_rM.resize(Index::set(-2),Index::set(gdata.mx[dir]+1));

	dudx_lP.resize(Index::set(-2),Index::set(gdata.mx[dir]+1));
	dudx_cP.resize(Index::set(-2),Index::set(gdata.mx[dir]+1));
	dudx_rP.resize(Index::set(-2),Index::set(gdata.mx[dir]+1));
	dudx_lM.resize(Index::set(-2),Index::set(gdata.mx[dir]+1));
	dudx_cM.resize(Index::set(-2),Index::set(gdata.mx[dir]+1));
	dudx_rM.resize(Index::set(-2),Index::set(gdata.mx[dir]+1));

	fac_cen = 13./3.;
	eps_WENO = 1.e-12;
	p_WENO = 2;
	use_Limiter = true;

	if (gdata.rank == 0) {
		cout << "  Using Reconst: " << "WENO" << " - " << "GENERIC" << " - dir " << dir << " - substep " << substep << endl;
	}

}


SingleReconstruction_WENO::~SingleReconstruction_WENO() {
}



void SingleReconstruction_WENO::get_Vals_EW(const Data &gdata, phys_fields_0D &xFieldsW,
		phys_fields_0D &xFieldsE, int ix, int iy, int iz)
{

	// Take into account shifted collocation points
	REAL shift(0.);
#if (GEOM != CARTESIAN)
#if (SHIFTED_COLLOCATION == TRUE)
	shift = sources->shift_Geom_WE(gdata, i);
#endif
#endif

#if (NON_LINEAR_GRID == CRONOS_ON)
		REAL delx = gdata.getCen_dx(0, ix );
#endif

	int q = qReconst;

#if (NON_LINEAR_GRID == CRONOS_ON)
	double AVal_x = gdata.om[q](ix,iy,iz) - wCx/12.*(dudxp_q - dudxm_q)*delx;
	double BVal_x = (wRx*(dudxp_q) + 0.5*wCx*(dudx0_q) + wLx*(dudxm_q))*delx;
	double CVal_x = 0.25*wCx*(dudxp_q - dudxm_q)*delx;
#else
	double AVal_x = gdata.om[q](ix,iy,iz) - wCx/12.*(dudxp_q - dudxm_q);
	double BVal_x = wRx*(dudxp_q) + 0.5*wCx*(dudx0_q) + wLx*(dudxm_q);
	double CVal_x = 0.25*wCx*(dudxp_q - dudxm_q);
#endif
	xFieldsW.uPri(q) = AVal_x - 0.5*BVal_x + CVal_x;
	xFieldsE.uPri(q) = xFieldsW.uPri(q) + BVal_x;

}

void SingleReconstruction_WENO::get_Vals_SN(const Data &gdata, phys_fields_0D &xFieldsS,
		phys_fields_0D &xFieldsN, int ix, int iy, int iz)
{

#if (NON_LINEAR_GRID == CRONOS_ON)
	REAL dely = gdata.getCen_dx(1, iy );	// iterate over all indices
#endif

	int q = qReconst;

#if (NON_LINEAR_GRID == CRONOS_ON)
	double AVal_y = gdata.om[q](ix,iy,iz) - wCy/12.*(dudyp_q - dudym_q)*dely;
	double BVal_y = (wRy*(dudyp_q) + 0.5*wCy*(dudy0_q) + wLy*(dudym_q))*dely;
	double CVal_y = 0.25*wCy*(dudyp_q - dudym_q)*dely;
#else
	double AVal_y = gdata.om[q](ix,iy,iz) - wCy/12.*(dudyp_q - dudym_q);
	double BVal_y = wRy*(dudyp_q) + 0.5*wCy*(dudy0_q) + wLy*(dudym_q);
	double CVal_y = 0.25*wCy*(dudyp_q - dudym_q);
#endif
	xFieldsS.uPri(q) = AVal_y - 0.5*BVal_y + CVal_y;
	xFieldsN.uPri(q) = xFieldsS.uPri(q) + BVal_y;

}

void SingleReconstruction_WENO::get_Vals_BT(const Data &gdata, phys_fields_0D &xFieldsB,
		phys_fields_0D &xFieldsT, int ix, int iy, int iz)
{

#if (NON_LINEAR_GRID == CRONOS_ON)
	REAL delz = gdata.getCen_dx(2, iz );	// iterate over all indices
#endif

	int q = qReconst;

#if (NON_LINEAR_GRID == CRONOS_ON)
	double AVal_z = gdata.om[q](ix,iy,iz) - wCz/12.*(dudzp_q - dudzm_q)*delz;
	double BVal_z = (wRz*(dudzp_q) + 0.5*wCz*(dudz0_q) + wLz*(dudzm_q))*delz;
	double CVal_z = 0.25*wCz*(dudzp_q - dudzm_q)*delz;
#else
	double AVal_z = gdata.om[q](ix,iy,iz) - wCz/12.*(dudzp_q - dudzm_q);
	double BVal_z = wRz*(dudzp_q) + 0.5*wCz*(dudz0_q) + wLz*(dudzm_q);
	double CVal_z = 0.25*wCz*(dudzp_q - dudzm_q);
#endif
	xFieldsB.uPri(q) = AVal_z - 0.5*BVal_z + CVal_z;
	xFieldsT.uPri(q) = xFieldsB.uPri(q) + BVal_z;

}


void SingleReconstruction_WENO::perpareDerivs(const Data &gdata, int ix, int iy, int iz) {
	//! Compute derivate from at given position

	// Compute local derivatives in all directions
	getDerivs(gdata, ix, iy, iz);

	// set default weights
	wLx = 0.25;
	wCx = 0.5;
	wRx = 0.25;

	wLy = 0.25;
	wCy = 0.5;
	wRy = 0.25;

	wLz = 0.25;
	wCz = 0.5;
	wRz = 0.25;

	if(use_Limiter) {
		// Compute smootheness indicators (Eq. (2.9) in KL2000)
		double IS_l_x = sqr(dudxm_q);
		double IS_r_x = sqr(dudxp_q);
		double IS_c_x = fac_cen*sqr(dudxp_q - dudxm_q) + sqr(dudx0_q);

		double IS_l_y = sqr(dudym_q);
		double IS_r_y = sqr(dudyp_q);
		double IS_c_y = fac_cen*sqr(dudyp_q - dudym_q) + sqr(dudy0_q);

		double IS_l_z = sqr(dudzm_q);
		double IS_r_z = sqr(dudzp_q);
		double IS_c_z = fac_cen*sqr(dudzp_q - dudzm_q) + sqr(dudz0_q);

		// Alpha factors (Eq. (2.8) in KL2000)
		double alpha_l_x = wLx / power(std::max(eps_WENO, IS_l_x), p_WENO);
		double alpha_c_x = wCx / power(std::max(eps_WENO, IS_c_x), p_WENO);
		double alpha_r_x = wRx / power(std::max(eps_WENO, IS_r_x), p_WENO);
		double alpha_x = alpha_l_x + alpha_c_x + alpha_r_x;

		double alpha_l_y = wLy/power(std::max(eps_WENO,IS_l_y),p_WENO);
		double alpha_c_y = wCy/power(std::max(eps_WENO,IS_c_y),p_WENO);
		double alpha_r_y = wRy/power(std::max(eps_WENO,IS_r_y),p_WENO);
		double alpha_y = alpha_l_y + alpha_c_y + alpha_r_y;

		double alpha_l_z = wLz/power(std::max(eps_WENO,IS_l_z),p_WENO);
		double alpha_c_z = wCz/power(std::max(eps_WENO,IS_c_z),p_WENO);
		double alpha_r_z = wRz/power(std::max(eps_WENO,IS_r_z),p_WENO);
		double alpha_z = alpha_l_z + alpha_c_z + alpha_r_z;

		// Compute resulting weights (Eq. (2.8) in KL2000)
		wLx = alpha_l_x/alpha_x;
		wCx = alpha_c_x/alpha_x;
		wRx = alpha_r_x/alpha_x;

		wLy = alpha_l_y/alpha_y;
		wCy = alpha_c_y/alpha_y;
		wRy = alpha_r_y/alpha_y;

		wLz = alpha_l_z/alpha_z;
		wCz = alpha_c_z/alpha_z;
		wRz = alpha_r_z/alpha_z;
	}

}





void SingleReconstruction_WENO::get_weights(const Data &gdata, const NumMatrix<REAL,1> &input,
		NumMatrix<REAL,1> &wL, NumMatrix<REAL,1> &wC, NumMatrix<REAL,1> &wR,
		NumMatrix<REAL,1> &dudxL, NumMatrix<REAL,1> &dudxC, NumMatrix<REAL,1> &dudxR, Buffer<REAL, 1>& inputSYCL) {
	for (int i = -2; i <= gdata.mx[dir]+1; ++i){
		// Compute local derivatives (left, centred, right)
		getDeriv(gdata, input, i, inputSYCL);
		dudxL(i) = dudxm;
		dudxC(i) = dudx0;
		dudxR(i) = dudxp;

		// set default weights
		wL(i) = 0.25;
		wC(i) = 0.5;
		wR(i) = 0.25;

		if(use_Limiter) {
			// Compute smootheness indicators (Eq. (2.9) in KL2000)
			double IS_l = sqr(dudxm);
			double IS_r = sqr(dudxp);
			double IS_c = fac_cen*sqr(dudxp - dudxm) + sqr(dudx0);

			// Alpha factors (Eq. (2.8) in KL2000)
			double alpha_l = wL(i)/power(std::max(eps_WENO,IS_l),p_WENO);
			double alpha_c = wC(i)/power(std::max(eps_WENO,IS_c),p_WENO);
			double alpha_r = wR(i)/power(std::max(eps_WENO,IS_r),p_WENO);
			double alpha = alpha_l + alpha_c + alpha_r;

			// Compute resulting weights (Eq. (2.8) in KL2000)
			wL(i) = alpha_l/alpha;
			wC(i) = alpha_c/alpha;
			wR(i) = alpha_r/alpha;

		}
	}

}


void SingleReconstruction_WENO::computeNormal(const Data &gdata,
		NumMatrix<REAL,1> &input, NumMatrix<REAL,1> &lhs, NumMatrix<REAL,1> &rhs, NumMatrix<REAL,1> &deriv, Buffer<REAL, 1>& inputSYCL)
{
	//std::cout << "SingleReconstruction_WENO" << std::endl << std::flush;
	// Compute weights for reconstruction (derivatives are stored as dudxm, dudx0, dudxp)
	get_weights(gdata, input, w_l, w_c, w_r, dudx_l, dudx_c, dudx_r, inputSYCL);

	for (int i = -2; i <= gdata.mx[dir]+1; ++i){

		// Take into account shifted collocation points
		REAL shift(0.);
#if (GEOM != CARTESIAN)
#if (SHIFTED_COLLOCATION == TRUE)
		if(dir == 0) {
			shift = sources->shift_Geom_WE(gdata, i);
		}
#endif
#endif


#if (NON_LINEAR_GRID == CRONOS_ON)
		REAL delx = gdata.getCen_dx(dir, i  );
		double AVal = input(i)-w_c(i)/12.*(dudxp - dudxm)*delx;
		double BVal = (w_r(i)*(dudxp) + 0.5*w_c(i)*(dudx0) + w_l(i)*(dudxm))*delx;
		double CVal = 0.25*w_c(i)*(dudxp - dudxm)*delx;
#else
		double AVal = input(i)-w_c(i)/12.*(dudx_r(i) - dudx_l(i));
		double BVal = w_r(i)*(dudx_r(i)) + 0.5*w_c(i)*(dudx_c(i)) + w_l(i)*(dudx_l(i));
		double CVal = 0.25*w_c(i)*(dudx_r(i) - dudx_l(i));
#endif
		lhs(i) = AVal - 0.5*BVal + CVal;
		rhs(i) = lhs(i) + BVal; // B can be stored as approximation for limited derivative
		deriv(i) = BVal;

	}

}

void SingleReconstruction_WENO::computePar(const Data &gdata,
		NumMatrix<REAL,1> &inputPar, NumMatrix<REAL,1> &lhs, NumMatrix<REAL,1> &rhs)
{
	for (int i = -2; i <= gdata.mx[dir]+1; ++i){
		lhs(i) = inputPar(i-1);
		rhs(i) = inputPar(i);
	}
}


void SingleReconstruction_WENO::computePerp(const Data &gdata,
		NumMatrix<REAL,1> &inputPerpPORIG, NumMatrix<REAL,1> &inputPerpMORIG,
		NumMatrix<REAL,1> &lhs, NumMatrix<REAL,1> &rhs, NumMatrix<REAL,1> &derivP, Buffer<REAL, 1>& inputPerpPSYCL, Buffer<REAL, 1>& inputPerpMSYCL)
{
	get_weights(gdata, inputPerpPORIG, w_lP, w_cP, w_rP, dudx_lP, dudx_cP, dudx_rP, inputPerpPSYCL);
	get_weights(gdata, inputPerpMORIG, w_lM, w_cM, w_rM, dudx_lM, dudx_cM, dudx_rM, inputPerpMSYCL);

	for (int i = -2; i <= gdata.mx[dir]+1; ++i){

		// Take into account shifted collocation points
		REAL shift(0.);
#if (GEOM != CARTESIAN)
#if (SHIFTED_COLLOCATION == TRUE)
		if(dir == 0) {
			shift = sources->shift_Geom_WE(gdata, i);
		}
#endif
#endif

#if (NON_LINEAR_GRID == CRONOS_ON)
		REAL delx = gdata.getCen_dx(dir, i  );
		double AValM = inputPerpMORIG(i)-w_cM(i)/12.*(dudx_rM(i) - dudx_lM(i))*delx;
		double BValM = (w_rM(i)*(dudx_rM(i)) + 0.5*w_cM(i)*(dudx_cM(i)) + w_lM(i)*(dudx_lM(i)))*delx;
		double CValM = 0.25*w_cM(i)*(dudx_rM(i) - dudx_lM(i))*delx;
		double AValP = inputPerpPORIG(i)-w_cP(i)/12.*(dudx_rP(i) - dudx_lP(i))*delx;
		double BValP = (w_rP(i)*(dudx_rP(i)) + 0.5*w_cP(i)*(dudx_cP(i)) + w_lP(i)*(dudx_lP(i)))*delx;
		double CValP = 0.25*w_cP(i)*(dudx_rP(i) - dudx_lP(i))*delx;
#else
		double AValM = inputPerpMORIG(i)-w_cM(i)/12.*(dudx_rM(i) - dudx_lM(i));
		double BValM = w_rM(i)*(dudx_rM(i)) + 0.5*w_cM(i)*(dudx_cM(i)) + w_lM(i)*(dudx_lM(i));
		double CValM = 0.25*w_cM(i)*(dudx_rM(i) - dudx_lM(i));
		double AValP = inputPerpPORIG(i)-w_cP(i)/12.*(dudx_rP(i) - dudx_lP(i));
		double BValP = w_rP(i)*(dudx_rP(i)) + 0.5*w_cP(i)*(dudx_cP(i)) + w_lP(i)*(dudx_lP(i));
		double CValP = 0.25*w_cP(i)*(dudx_rP(i) - dudx_lP(i));
#endif
		lhs(i) = 0.5*(AValM + AValP - 0.5*(BValM + BValP) + (CValM + CValP));
		rhs(i) = lhs(i) + 0.5*(BValM + BValP); // B can be stored as approximation for limited derivative
		derivP(i) = 0.5*(BValM + BValP);

	}  
}




//Reconstruction2D_WENO::Reconstruction2D_WENO(const Data &gdata, const int &dir)  : Reconstruction2D(gdata, dir){
//	fac_cen = 13./3.;
//	eps_WENO = 1.e-12;
//	p_WENO = 2;
//	use_Limiter = true;
//}
//
//
//
//void Reconstruction2D_2nd::computeNormal(const Data &gdata,
//                                     NumMatrix<REAL,2> &input,
//                                     NumMatrix<REAL,2> &deriv_dir0,
//                                     NumMatrix<REAL,2> &deriv_dir1,
//                                     NumMatrix<REAL,2> &recLL,
//                                     NumMatrix<REAL,2> &recLR,
//                                     NumMatrix<REAL,2> &recRL,
//                                     NumMatrix<REAL,2> &recRR)
//{
//
//#ifdef CRONOS_SAVEMEM
//	/* Will only compute derivatives if not saved from earlier
//	   Otherwise they are taken from input.
//	*/
//	getDerivs(gdata, input, deriv_dir0, deriv_dir1);
//#endif
//	for (int i = -1; i <= gdata.mx[dir0]+1; ++i){
//		for (int j = -1; j <= gdata.mx[dir1]+1; ++j){
//
//#if (NON_LINEAR_GRID == CRONOS_ON)
//			REAL delx0 = gdata.getCen_dx(dir0, i  );
//			REAL delx1 = gdata.getCen_dx(dir1, j  );
//
//			recLL(i,j) = input(i,j) - 0.5*(deriv_dir0(i,j)*delx0 +
//			                               deriv_dir1(i,j)*delx1);
//			recLR(i,j) = recLL(i,j) + deriv_dir1(i,j)*delx1;
//			recRL(i,j) = recLL(i,j) + deriv_dir0(i,j)*delx0;
//			recRR(i,j) = recRL(i,j) + deriv_dir1(i,j)*delx1;
//#else
//			recLL(i,j) = input(i,j) - 0.5*deriv_dir0(i,j) - 0.5*deriv_dir1(i,j);
//			recLR(i,j) = recLL(i,j) + deriv_dir1(i,j);
//			recRL(i,j) = recLL(i,j) + deriv_dir0(i,j);
//			recRR(i,j) = recRL(i,j) + deriv_dir1(i,j);
//#endif
//
//		}
//	}
//}
//
//
//void Reconstruction2D_2nd::computePerp(const Data &gdata,
//                                   NumMatrix<REAL,2> &input,
//                                   NumMatrix<REAL,2> &deriv_inp,
//                                   const int &dir,
//                                   NumMatrix<REAL,2> &recLL,
//                                   NumMatrix<REAL,2> &recLR,
//                                   NumMatrix<REAL,2> &recRL,
//                                   NumMatrix<REAL,2> &recRR)
//{
//	assert(dir == dir0 || dir == dir1);
//
//#ifdef CRONOS_SAVEMEM
//	getDeriv_limit(gdata, input, deriv_inp, dir);
//#endif
//	if(dir == dir0) {
//		for (int i = -1; i <= gdata.mx[dir0]+1; ++i){
//			for (int j = -1; j <= gdata.mx[dir1]+1; ++j){
//
//#if (NON_LINEAR_GRID == CRONOS_ON)
//				REAL delxm = gdata.getCen_dx(dir, i);
//				recLL(i,j) = input(i,j-1) - 0.5*deriv_inp(i,j-1)*delxm;
//				recRL(i,j) = recLL(i,j) + deriv_inp(i,j-1)*delxm;
//				REAL delxp = delxm; // For orthogonal grids
//				recLR(i,j) = input(i,  j) - 0.5*deriv_inp(i,j  )*delxp;
//				recRR(i,j) = recLR(i,j) + deriv_inp(i,j  )*delxp;
//#else
//				recLL(i,j) = input(i,j-1) - 0.5*deriv_inp(i,j-1);
//				recRL(i,j) = recLL(i,j) + deriv_inp(i,j-1);
//				recLR(i,j) = input(i,  j) - 0.5*deriv_inp(i,j  );
//				recRR(i,j) = recLR(i,j) + deriv_inp(i,j  );
//#endif
//
//			}
//		}
//	} else if (dir == dir1) {
//		for (int i = -1; i <= gdata.mx[dir0]+1; ++i){
//			for (int j = -1; j <= gdata.mx[dir1]+1; ++j){
//#if (NON_LINEAR_GRID == CRONOS_ON)
//				REAL delxm = gdata.getCen_dx(dir, j);
//				recLL(i,j) = input(i-1,j) - 0.5*deriv_inp(i-1,j)*delxm;
//				recLR(i,j) = recLL(i,j) + deriv_inp(i-1,j)*delxm;
//				REAL delxp = delxm;
//				recRL(i,j) = input(i  ,j) - 0.5*deriv_inp(i  ,j)*delxp;
//				recRR(i,j) = recRL(i,j) + deriv_inp(i  ,j)*delxp;
//#else
//				recLL(i,j) = input(i-1,j) - 0.5*deriv_inp(i-1,j);
//				recLR(i,j) = recLL(i,j) + deriv_inp(i-1,j);
//				recRL(i,j) = input(i  ,j) - 0.5*deriv_inp(i  ,j);
//				recRR(i,j) = recRL(i,j) + deriv_inp(i  ,j);
//#endif
//			}
//		}
//	}
//
//}
//
//
//
//void Reconstruction2D_2nd::getDerivs(const Data &gdata,
//                                 const NumMatrix<REAL,2> &input,
//                                 NumMatrix<REAL,2> &deriv_dir0,
//                                 NumMatrix<REAL,2> &deriv_dir1) {
//
//	for (int i = -2; i <= gdata.mx[dir0]+1; ++i){
//		for (int j = -2; j <= gdata.mx[dir1]+1; ++j){
//
//			// Compute derivatives and store in dudxp, dudx0, dudxm
//			getDeriv(gdata, input, dir0, i, j);
//
//			deriv_dir0(i,j) =  Limiter.compute(dudxp,dudx0,dudxm);
//		}
//	}
//
//
//	for (int i = -2; i <= gdata.mx[dir0]+1; ++i){
//		for (int j = -2; j <= gdata.mx[dir1]+1; ++j){
//
//			// Compute derivatives and store in dudyp, dudy0, dudym
//			getDeriv(gdata, input, dir1, i, j);
//
//			deriv_dir1(i,j) =  Limiter.compute(dudyp,dudy0,dudym);
//		}
//	}
//
//}
//
//
//void Reconstruction2D_WENO::get_weights(const Data &gdata, const NumMatrix<REAL,2> &input,
//		NumMatrix<REAL,1> &wL, NumMatrix<REAL,1> &wC, NumMatrix<REAL,1> &wR, int dir) {
//	if(dir == dir0) {
//		for (int i = -2; i <= gdata.mx[dir0]+1; ++i){
//			for (int j = -2; j <= gdata.mx[dir1]+1; ++j){
//
//				// Compute derivatives and store in dudxp, dudx0, dudxm
//				getDeriv(gdata, input, dir0, i, j);
//
//				// set default weights
//				wL(i) = 0.25;
//				wC(i) = 0.5;
//				wR(i) = 0.25;
//
//				if(use_Limiter) {
//					// Compute smootheness indicators (Eq. (2.9) in KL2000)
//					double IS_l = sqr(dudxm);
//					double IS_r = sqr(dudxp);
//					double IS_c = fac_cen*sqr(dudxp - dudxm) + sqr(dudx0);
//
//					// Alpha factors (Eq. (2.8) in KL2000)
//					double alpha_l = wL(i)/power(max(eps_WENO,IS_l),p_WENO);
//					double alpha_c = wC(i)/power(max(eps_WENO,IS_c),p_WENO);
//					double alpha_r = wR(i)/power(max(eps_WENO,IS_r),p_WENO);
//					double alpha = alpha_l + alpha_c + alpha_r;
//
//					// Compute resulting weights (Eq. (2.8) in KL2000)
//					wL(i) = alpha_l/alpha;
//					wC(i) = alpha_l/alpha;
//					wR(i) = alpha_l/alpha;
//
//				}
//
//				deriv(i,j) =  Limiter.compute(dudxp,dudx0,dudxm);
//			}
//		}
//	} else if (dir == dir1) {
//		for (int i = -2; i <= gdata.mx[dir0]+1; ++i){
//			for (int j = -2; j <= gdata.mx[dir1]+1; ++j){
//
//				// Compute derivatives and store in dudyp, dudy0, dudym
//				getDeriv(gdata, input, dir1, i, j);
//
//				deriv(i,j) =  Limiter.compute(dudyp,dudy0,dudym);
//			}
//		}
//	}
//}


