#ifndef CRONOS_RECONST_BLOCK_H
#define CRONOS_RECONST_BLOCK_H = 1

#include <vector>

#include "matrix.H"
#include "data.H"
#include "limiter.H"
#include "fields_local.H"
#include "standalone_usage.H"
#include "utils.H"
#include "direction.H"

// #if(STANDALONE_USAGE == TRUE)
// #include "direction.H"
// #endif

namespace gpu {

	struct Dud {
		double _m;
		double _0;
		double _p;
	};
	
	inline double minmodCompute(const double Theta, double& up, double& uc, double& um) {
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
	auto prepareDerivs(AccessorType& omAcc, const double Theta, size_t ix, size_t iy, size_t iz) {

		auto [dudx, dudy, dudz] = getDerivs(omAcc, ix, iy, iz);

		double deriv_x = minmodCompute(Theta, dudx._p, dudx._0, dudx._m);
		double deriv_y = minmodCompute(Theta, dudy._p, dudy._0, dudy._m);
		double deriv_z = minmodCompute(Theta, dudz._p, dudz._0, dudz._m);

		return std::tuple{ deriv_x, deriv_y, deriv_z };
	}
	
	template<typename AccessorType>
	auto prepareDerivsOffset(AccessorType& omAcc, const double Theta, size_t ix, size_t iy, size_t iz) {

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
	void compute(AccessorType& omAcc, double physValsAcc[][N_OMINT], const double Theta, size_t ix, size_t iy, size_t iz, int q, Direction dir) {

		auto [deriv_x, deriv_y, deriv_z] = prepareDerivs(omAcc, Theta, ix, iy, iz);
		
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
	auto computeWSB(AccessorType& omAcc, const double Theta, size_t ix, size_t iy, size_t iz) {

		auto [deriv_x, deriv_y, deriv_z] = prepareDerivs(omAcc, Theta, ix, iy, iz);
		double shift(0.);

		double uPriWest = omAcc[ix][iy][iz] - (0.5 + shift) * deriv_x;
		double uPriSouth = omAcc[ix][iy][iz] - 0.5 * deriv_y;
		double uPriBottom = omAcc[ix][iy][iz] - 0.5 * deriv_z;
		
		return std::tuple{uPriWest, uPriSouth, uPriBottom};
	}
	
	template<typename AccessorType>
	auto computeENT(AccessorType& omAcc, const double Theta, size_t ix, size_t iy, size_t iz) {

		auto [deriv_x, deriv_y, deriv_z] = prepareDerivsOffset(omAcc, Theta, ix, iy, iz);
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

#endif