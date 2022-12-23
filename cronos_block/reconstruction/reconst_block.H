#ifndef CRONOS_RECONST_BLOCK_H
#define CRONOS_RECONST_BLOCK_H = 1

#include <vector>

#include "matrix.H"
#include "data.H"
#include "limiter.H"
#include "fields_local.H"
#include "standalone_usage.H"
#include "utils.H"

#if(STANDALONE_USAGE == TRUE)
#include "direction.H"
#endif

namespace gpu {

	struct Dud {
		double _m;
		double _0;
		double _p;
	};
	
	inline double minmodCompute(double Theta, double& up, double& uc, double& um) {
		up *= Theta;
		um *= Theta;
		if (up * um > 0.) {
			if (up > 0) {
				return cl::sycl::fmin(up, (cl::sycl::fmin(uc, um)));
			} else if (up < 0) {
				return cl::sycl::fmax(up, (cl::sycl::fmax(uc, um)));
			} else {
				return 0.0;
			}
		} else {
			return 0.;
		}
	}

	template<typename AccessorType>
	auto getDerivs(AccessorType& omAcc, size_t ix, size_t iy, size_t iz) {
		//size_t q = qReconst;
		Dud x;
		Dud y;
		Dud z;

		x._m = omAcc[ix][iy][iz] - omAcc[ix - 1][iy][iz];
		x._p = omAcc[ix + 1][iy][iz] - omAcc[ix][iy][iz];
		x._0 = 0.5 * (x._m + x._p);

		y._m = omAcc[ix][iy][iz] - omAcc[ix][iy - 1][iz];
		y._p = omAcc[ix][iy + 1][iz] - omAcc[ix][iy][iz];
		y._0 = 0.5 * (y._m + y._p);

		z._m = omAcc[ix][iy][iz] - omAcc[ix][iy][iz - 1];
		z._p = omAcc[ix][iy][iz + 1] - omAcc[ix][iy][iz];
		z._0 = 0.5 * (z._m + z._p);

		return std::tuple{x, y, z};
	}
	
	template<typename AccessorType>
	auto prepareDerivs(AccessorType& omAcc, size_t ix, size_t iy, size_t iz) {

		auto [dudx, dudy, dudz] = getDerivs(omAcc, ix, iy, iz);

		//dproksch TODO: read from project parameters
		double Theta = 1;

		double deriv_x = minmodCompute(Theta, dudx._p, dudx._0, dudx._m);
		double deriv_y = minmodCompute(Theta, dudy._p, dudy._0, dudy._m);
		double deriv_z = minmodCompute(Theta, dudz._p, dudz._0, dudz._m);

		return std::tuple{ deriv_x, deriv_y, deriv_z };
	}
	
	template<typename AccessorType>
	auto prepareDerivsOffset(AccessorType& omAcc, size_t ix, size_t iy, size_t iz) {

		Dud x;
		Dud y;
		Dud z;

		x._m = omAcc[ix - 1][iy][iz] - omAcc[ix - 2][iy][iz];
		x._p = omAcc[ix][iy][iz] - omAcc[ix - 1][iy][iz];
		x._0 = 0.5 * (x._m + x._p);

		y._m = omAcc[ix][iy - 1][iz] - omAcc[ix][iy - 2][iz];
		y._p = omAcc[ix][iy][iz] - omAcc[ix][iy - 1][iz];
		y._0 = 0.5 * (y._m + y._p);

		z._m = omAcc[ix][iy][iz - 1] - omAcc[ix][iy][iz - 2];
		z._p = omAcc[ix][iy][iz] - omAcc[ix][iy][iz - 1];
		z._0 = 0.5 * (z._m + z._p);

		//dproksch TODO: read from project parameters
		double Theta = 1;

		double deriv_x = minmodCompute(Theta, x._p, x._0, x._m);
		double deriv_y = minmodCompute(Theta, y._p, y._0, y._m);
		double deriv_z = minmodCompute(Theta, z._p, z._0, z._m);

		return std::tuple{ deriv_x, deriv_y, deriv_z };
	}

	template<typename AccessorType>
	void get_Vals_EW(AccessorType& omAcc, double uPri[][N_OMINT], size_t ix, size_t iy, size_t iz, int q, double deriv_x) {
		double shift(0.);
		uPri[FaceWest][q] = omAcc[ix][iy][iz] - (0.5 + shift) * deriv_x;
		uPri[FaceEast][q] = uPri[FaceWest][q] + deriv_x;
	}

	template<typename AccessorType>
	void get_Vals_SN(AccessorType& omAcc, double uPri[][N_OMINT], size_t ix, size_t iy, size_t iz, int q, double deriv_y) {
		double shift(0.);
		uPri[FaceSouth][q] = omAcc[ix][iy][iz] - 0.5 * deriv_y;
		uPri[FaceNorth][q] = uPri[FaceSouth][q] + deriv_y;
	}

	template<typename AccessorType>
	void get_Vals_BT(AccessorType& omAcc, double uPri[][N_OMINT], size_t ix, size_t iy, size_t iz, int q, double deriv_z) {
		double shift(0.);
		uPri[FaceBottom][q] = omAcc[ix][iy][iz] - 0.5 * deriv_z;
		uPri[FaceTop][q] = uPri[FaceBottom][q] + deriv_z;
	}
	
	template<typename AccessorType>
	void compute(AccessorType& omAcc, double physValsAcc[][N_OMINT], size_t ix, size_t iy, size_t iz, int q, Direction dir) {

		auto [deriv_x, deriv_y, deriv_z] = prepareDerivs(omAcc, ix, iy, iz);
		
		if (dir == DirAll || dir == DirX) {
			get_Vals_EW(omAcc, physValsAcc, ix, iy, iz, q, deriv_x);
		}

		if (dir == DirAll || dir == DirY) {
			get_Vals_SN(omAcc, physValsAcc, ix, iy, iz, q, deriv_y);
		}

		if (dir == DirAll || dir == DirZ) {
			get_Vals_BT(omAcc, physValsAcc, ix, iy, iz, q, deriv_z);
		}
		
	}

	template<typename AccessorType>
	auto computeWSB(AccessorType& omAcc, size_t ix, size_t iy, size_t iz) {

		auto [deriv_x, deriv_y, deriv_z] = prepareDerivs(omAcc, ix, iy, iz);
		double shift(0.);

		double uPriWest = omAcc[ix][iy][iz] - (0.5 + shift) * deriv_x;
		double uPriSouth = omAcc[ix][iy][iz] - 0.5 * deriv_y;
		double uPriBottom = omAcc[ix][iy][iz] - 0.5 * deriv_z;
		
		return std::tuple{uPriWest, uPriSouth, uPriBottom};
	}
	
	template<typename AccessorType>
	auto computeENT(AccessorType& omAcc, size_t ix, size_t iy, size_t iz) {

		auto [deriv_x, deriv_y, deriv_z] = prepareDerivsOffset(omAcc, ix, iy, iz);
		double shift(0.);

		double uPriEast = omAcc[ix - 1][iy][iz] + (0.5 - shift) * deriv_x;
		double uPriNorth = omAcc[ix][iy - 1][iz] + 0.5 * deriv_y;
		double uPriTop = omAcc[ix][iy][iz - 1] + 0.5 * deriv_z;
		
		return std::tuple{uPriEast, uPriNorth, uPriTop};
	}

	template<typename AccessorType>
	auto computeLowerCorner(AccessorType& omAcc, size_t ix, size_t iy, size_t iz) {

		auto [deriv_x, deriv_y, deriv_z] = prepareDerivs(omAcc, ix, iy, iz);
		double shift(0.);

		double uPriWest = omAcc[ix][iy][iz] - (0.5 + shift) * deriv_x;
		double uPriSouth = omAcc[ix][iy][iz] - 0.5 * deriv_y;
		double uPriBottom = omAcc[ix][iy][iz] - 0.5 * deriv_z;
		double uPriEast = uPriWest + deriv_x;
		double uPriNorth = uPriSouth + deriv_y;
		double uPriTop = uPriBottom + deriv_z;
		
		return std::tuple{uPriWest, uPriEast, uPriSouth, uPriNorth, uPriBottom, uPriTop};
	}

}

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