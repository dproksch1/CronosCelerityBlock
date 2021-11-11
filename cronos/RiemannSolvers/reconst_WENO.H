
#ifndef CRONOS_RECONST_WENO_H
#define CRONOS_RECONST_WENO_H = 1

/*===================================================================================*/
/*!
 * \file
 * \brief Third order WENO reconstruction
 *
 * The specific implementation is adapted from:
 * Kurganov & Levy (2000) SISC 22, 1461 --> KL2000
 *
 *  */
/*===================================================================================*/

#include "reconst.H"
#include "limiter.H"
#include "fields_local.H"
#include "data.H"
#include "matrix.H"
#include <vector>

class SingleReconstruction_WENO: public SingleReconstruction {
	/** Class for second order spatial reconstruction
	 * --> using different limiters
	*/
public:
	SingleReconstruction_WENO(){}
	SingleReconstruction_WENO(const Data &, int dir, int substep);
	SingleReconstruction_WENO(const Data &, const CronosFluid &fluid, int dir, int qReconst, int substep);
	virtual ~SingleReconstruction_WENO();


	/**
	 * Compute limited version of derivatives for block-structured version of code
	 * */
	void prepareDerivs(const Data &, int ix, int iy, int iz);
	void computeNormal(const Data &, NumMatrix<REAL,1> &, NumMatrix<REAL,1> &,
	                   NumMatrix<REAL,1> &, NumMatrix<REAL,1> &, Buffer<REAL, 1>&);
	void computePar(const Data &, NumMatrix<REAL,1> &, NumMatrix<REAL,1> &,
	                NumMatrix<REAL,1> &);
	void computePerp(const Data &, NumMatrix<REAL,1> &, NumMatrix<REAL,1> &,
	                 NumMatrix<REAL,1> &, NumMatrix<REAL,1> &,
	                 NumMatrix<REAL,1> &, Buffer<REAL, 1>&, Buffer<REAL, 1>&);
	void get_Vals_EW(const Data &, phys_fields_0D &xFieldsW,
			phys_fields_0D &xFieldsE, int, int, int);
	void get_Vals_SN(const Data &, phys_fields_0D &xFieldsS,
			phys_fields_0D &xFieldsN, int, int, int);
	void get_Vals_BT(const Data &, phys_fields_0D &xFieldsB,
			phys_fields_0D &xFieldsT, int, int, int);
private:
	/**
	 * Compute smootheness indicators and weights for WENO reconstruction
	 *  */
	void get_weights(const Data &, const NumMatrix<REAL,1> &input,
			NumMatrix<REAL,1> &wL, NumMatrix<REAL,1> &wC, NumMatrix<REAL,1> &wR,
			NumMatrix<REAL,1> &dudxL, NumMatrix<REAL,1> &dudxC, NumMatrix<REAL,1> &dudxR, Buffer<REAL, 1>&);
//	template<typename T>
//	T Reconstruction_WENO::power(T base, int expo);
	/** Weights
	 * */
	NumMatrix<REAL, 1>  w_l, w_c, w_r;
	NumMatrix<REAL, 1>  w_lP, w_cP, w_rP;
	NumMatrix<REAL, 1>  w_lM, w_cM, w_rM;

	NumMatrix<REAL,1> dudx_l, dudx_c, dudx_r;
	NumMatrix<REAL,1> dudx_lP, dudx_cP, dudx_rP;
	NumMatrix<REAL,1> dudx_lM, dudx_cM, dudx_rM;

	REAL wLx, wCx, wRx;
	REAL wLy, wCy, wRy;
	REAL wLz, wCz, wRz;

	double eps_WENO, fac_cen;
	int p_WENO;
	bool use_Limiter;
};

//class Reconstruction2D_WENO: public Reconstruction2D {
//public:
//	Reconstruction2D_WENO(const Data &, const int &);
//	virtual ~Reconstruction2D_WENO(){}
//	void computeNormal(const Data &,
//	                   NumMatrix<REAL,2> &, NumMatrix<REAL,2> &,
//	                   NumMatrix<REAL,2> &, NumMatrix<REAL,2> &,
//	                   NumMatrix<REAL,2> &, NumMatrix<REAL,2> &,
//	                   NumMatrix<REAL,2> &);
//	void computePerp(const Data &, NumMatrix<REAL,2> &, NumMatrix<REAL,2> &,
//	                 const int &,
//	                 NumMatrix<REAL,2> &, NumMatrix<REAL,2> &,
//	                 NumMatrix<REAL,2> &, NumMatrix<REAL,2> &);
//	void getDerivs(const Data &, const NumMatrix<REAL,2> &,
//	               NumMatrix<REAL,2> &, NumMatrix<REAL,2> &);
//	void get_weights(const Data &, const NumMatrix<REAL,2> &,
//			NumMatrix<REAL,1> &wL, NumMatrix<REAL,1> &wC, NumMatrix<REAL,1> &wR, int);
//private:
//	double eps_WENO, fac_cen;
//	int p_WENO;
//	bool use_Limiter;
//};

#endif
