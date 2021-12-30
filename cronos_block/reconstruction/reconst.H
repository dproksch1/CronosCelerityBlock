#ifndef CRONOS_RECONST_H
#define CRONOS_RECONST_H = 1

#include "matrix.H"
#include "data.H"
#include "limiter.H"
#include "fields_local.H"
#include "direction.H"
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

	void prepareDerivs(const Data& gdata, int ix, int iy, int iz);

	void get_Vals_EW(const Data &, phys_fields_0D &xFieldsW,
			phys_fields_0D &xFieldsE, int, int, int);
	void get_Vals_SN(const Data &, phys_fields_0D &xFieldsS,
			phys_fields_0D &xFieldsN, int, int, int);
	void get_Vals_BT(const Data &, phys_fields_0D &xFieldsB,
			phys_fields_0D &xFieldsT, int, int, int);
      int qReconst;      
protected:
	/**
	 * Compute left-handed, centred, and right-handed derivatives in all spatial directions
	 * --> block-structured version
	 * */
	void getDerivs(const Data &, int ix, int iy, int iz);
	int dir;
	
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
};

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