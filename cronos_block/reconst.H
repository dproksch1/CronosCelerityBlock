#ifndef CRONOS_RECONST_H
#define CRONOS_RECONST_H = 1

#include "matrix.H"
#include "data.H"
#include "limiter.H"
#include "fields_local.H"
#include <vector>

class SingleReconstruction{
public:

	enum DudIndex {
		_p,
		_0, 
		_m,
		Max_Index,
	};

	enum DudDir {
		_x,
		_y,
		_z,
		Max_Dir,
	};

	SingleReconstruction() {}
	SingleReconstruction(const Data &, int dir, int substep);
	SingleReconstruction(const Data &, const CronosFluid &fluid, int dir, int qReconst, int substep);
	virtual ~SingleReconstruction();

	/**
	 * Compute reconstruction of cell-centered variable
	 * */
	void computeNormal(Queue& /*queue*/, const Data& /*gdata*/, NumMatrix<REAL, 1>& /*inputORIG*/, NumMatrix<REAL, 1>& /*lhs*/,
		NumMatrix<REAL, 1>& /*rhs*/, NumMatrix<REAL, 1>& /*derivORIG*/, Buffer<REAL, 1> /*lhsSYCL*/, Buffer<REAL, 1> /*rhsSYCL*/, Buffer<REAL, 1> /*inputSYCL*/, Buffer<REAL, 1> /*derivSYCL*/) = 0;

	/**
	 * Compute reconstruction of face-centered variable in direction of respective cell face
	 * */
	void computePar(const Data &, NumMatrix<REAL,1> &, NumMatrix<REAL,1> &,
	                NumMatrix<REAL,1> &) = 0;

	/**
	 * Compute reconstruction of face-centered variable in perpendicular direction
	 * */
	void computePerp(Queue& /*queue*/, const Data& /*gdata*/, NumMatrix<REAL, 1>& /*inputPerpPORIG*/, NumMatrix<REAL, 1>& /*inputPerpMORIG*/,
							 NumMatrix<REAL, 1>& /*lhs*/, NumMatrix<REAL, 1>& /*rhs*/, NumMatrix<REAL, 1>& /*derivP*/,
							 Buffer<REAL, 1> /*inputPerpPSYCL*/, Buffer<REAL, 1> /*inputPerpMSYCL*/, Buffer<REAL, 1> /*derivSYCL*/) = 0;

	void prepareDerivs(const Data& gdata, int ix, int iy, int iz) {};

	void get_Vals_EW(const Data &, phys_fields_0D &xFieldsW,
			phys_fields_0D &xFieldsE, int, int, int) = 0;
	void get_Vals_SN(const Data &, phys_fields_0D &xFieldsS,
			phys_fields_0D &xFieldsN, int, int, int) = 0;
	void get_Vals_BT(const Data &, phys_fields_0D &xFieldsB,
			phys_fields_0D &xFieldsT, int, int, int) = 0;
            
protected:
	/**
	 * Compute left-handed, centred, and right-handed derivatives in given direction
	 * */
	void getDeriv(Queue& /*queue*/, const Data& /*gdata*/, const NumMatrix<REAL, 1>& /*inputORIG*/, int /*iPos*/, NumMatrix<REAL, 1>& /*dudxORIG*/, Buffer<REAL, 1> /*inputSYCL*/, Buffer<REAL, 1> /*dudxSYCL*/);

	/**
	 * Compute left-handed, centred, and right-handed derivatives in all spatial directions
	 * --> block-structured version
	 * */
	void getDerivs(const Data &, int ix, int iy, int iz);
	int dir;
	int qReconst;
	int substep;

	NumMatrix<REAL, 1> derivPerp, derivM;
	REAL deriv_x, deriv_y, deriv_z;
	REAL dudxp_q, dudx0_q, dudxm_q;
	REAL dudyp_q, dudy0_q, dudym_q;
	REAL dudzp_q, dudz0_q, dudzm_q;
	//std::vector<Buffer<REAL, 1>> dud_q = std::vector<Buffer<REAL, 1>>(DudDir::Max_Dir, Range<1>(DudIndex::Max_Index));

	//REAL dudxp, dudx0, dudxm;
	//Buffer<REAL, 1> dudx;
    
private:

 	limiter Limiter;
}

class Reconstruction{
public:
	Reconstruction(){}
	Reconstruction(const Data &, int dir, const CronosFluid &fluid, int substep=-1);
	Reconstruction(const Data &, int dir, int num, int substep=-1);
	virtual ~Reconstruction();

	/**
	 * Compute reconstruction for block-structured code
	 * */
	void compute(const Data &, std::vector<phys_fields_0D> &allFields, int, int, int, Direction dir = (Direction)-1);

protected:
	void set_singleReconstructions(const Data &);
	void set_singleReconstructions(const Data &, const CronosFluid &fluid);

	int dir;
	int substep;
	std::vector<int>::iterator iter;
	std::vector<int> ListNormal;
	std::vector<int> ListParallel;
	std::vector<int> ListPerp;

	std::vector<SingleReconstruction*> ListReconstructionNormal;
	std::vector<SingleReconstruction*> ListReconstructionPar;
	std::vector<SingleReconstruction*> ListReconstructionPerp;

};

#endif