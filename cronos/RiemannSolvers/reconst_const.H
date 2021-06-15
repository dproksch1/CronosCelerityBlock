
#ifndef CRONOS_RECONST_CONST_H
#define CRONOS_RECONST_CONST_H = 1

#include "reconst.H"
#include "limiter.H"
#include "fields_local.H"
#include "data.H"
#include "matrix.H"
#include <vector>

class SingleReconstruction_constant: public SingleReconstruction {
	/** Class for second order spatial reconstruction
	 * --> using different limiters
	*/
public:
	SingleReconstruction_constant(){}
	SingleReconstruction_constant(const Data &, int dir, int substep);
	SingleReconstruction_constant(const Data &, const CronosFluid &fluid, int dir, int qReconst, int substep);

	virtual ~SingleReconstruction_constant();

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
};

class Reconstruction2D_constant: public Reconstruction2D {
public:
	Reconstruction2D_constant(const Data &, const int &);
	virtual ~Reconstruction2D_constant(){}
	void computeNormal(const Data &,
	                   NumMatrix<REAL,2> &, NumMatrix<REAL,2> &,
	                   NumMatrix<REAL,2> &, NumMatrix<REAL,2> &,
	                   NumMatrix<REAL,2> &, NumMatrix<REAL,2> &,
	                   NumMatrix<REAL,2> &);
	void computePerp(const Data &, NumMatrix<REAL,2> &, NumMatrix<REAL,2> &,
	                 const int &,
	                 NumMatrix<REAL,2> &, NumMatrix<REAL,2> &,
	                 NumMatrix<REAL,2> &, NumMatrix<REAL,2> &);
	void getDerivs(const Data &, const NumMatrix<REAL,2> &,
	               NumMatrix<REAL,2> &, NumMatrix<REAL,2> &);
	void getDeriv_limit(const Data &, const NumMatrix<REAL,2> &, NumMatrix<REAL,2> &,
	              int);
};

#endif
