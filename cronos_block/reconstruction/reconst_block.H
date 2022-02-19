#ifndef CRONOS_RECONST_BLOCK_H
#define CRONOS_RECONST_BLOCK_H = 1

#include <vector>

#include "matrix.H"
#include "data.H"
#include "limiter.H"
#include "fields_local.H"
#include "standalone_usage.H"

#if(STANDALONE_USAGE == TRUE)
#include "direction.H"
#endif

class SingleReconstruction_Block{
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

	SingleReconstruction_Block() {}
	SingleReconstruction_Block(const Data &, int dir, int substep);
	SingleReconstruction_Block(const Data &, const CronosFluid &fluid, int dir, int qReconst, int substep);
	virtual ~SingleReconstruction_Block();

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

	NumMatrix<double, 1> derivPerp, derivM;
	double deriv_x, deriv_y, deriv_z;
	double dudxp_q, dudx0_q, dudxm_q;
	double dudyp_q, dudy0_q, dudym_q;
	double dudzp_q, dudz0_q, dudzm_q;
	//std::vector<Buffer<double, 1>> dud_q = std::vector<Buffer<double, 1>>(DudDir::Max_Dir, Range<1>(DudIndex::Max_Index));

	//double dudxp, dudx0, dudxm;
	//Buffer<double, 1> dudx;
    
private:

 	limiter Limiter;
};

class Reconstruction_Block{
public:
	Reconstruction_Block(){}
	Reconstruction_Block(const Data &, int dir, const CronosFluid &fluid, int substep=-1);
	Reconstruction_Block(const Data &, int dir, int num, int substep=-1);
	virtual ~Reconstruction_Block();

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

	std::vector<SingleReconstruction_Block*> ListReconstructionNormal;
	std::vector<SingleReconstruction_Block*> ListReconstructionPar;
	std::vector<SingleReconstruction_Block*> ListReconstructionPerp;

};

#endif