/**
 * @file
 * @brief GPU-compatible second-order reconstruction for the main kernel
 *
 * This file provides a GPU-compatible version of second-order reconstruction specifically
 *  scoped for the Celerity Cronos main kernel.
 */

#ifndef CRONOS_RECONST_BLOCK_H
#define CRONOS_RECONST_BLOCK_H = 1

#include <vector>

#include "matrix.H"
#include "data.H"
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
	
	//! @brief Applies a minmod limiter to the given values
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

	//! @brief Computes the reconstruction derivatives
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
	
	//! @brief Computes the reconstruction derivatives and applies a minmod limiter
	template<typename AccessorType>
	auto prepareDerivs(AccessorType& omAcc, const double Theta, size_t ix, size_t iy, size_t iz) {

		auto [dudx, dudy, dudz] = getDerivs(omAcc, ix, iy, iz);

		double deriv_x = minmodCompute(Theta, dudx._p, dudx._0, dudx._m);
		double deriv_y = minmodCompute(Theta, dudy._p, dudy._0, dudy._m);
		double deriv_z = minmodCompute(Theta, dudz._p, dudz._0, dudz._m);

		return std::tuple{ deriv_x, deriv_y, deriv_z };
	}
	
	//! @brief Computes the reconstruction derivatives with an offset of -1 and applies a minmod limiter
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

	//! @brief Reconstructs the vertices for the directions West, South and Bottom
	template<typename AccessorType>
	auto computeWSB(AccessorType& omAcc, const double Theta, size_t ix, size_t iy, size_t iz) {

		auto [deriv_x, deriv_y, deriv_z] = prepareDerivs(omAcc, Theta, ix, iy, iz);
		double shift(0.);

		double uPriWest = omAcc[ix][iy][iz] - (0.5 + shift) * deriv_x;
		double uPriSouth = omAcc[ix][iy][iz] - 0.5 * deriv_y;
		double uPriBottom = omAcc[ix][iy][iz] - 0.5 * deriv_z;
		
		return std::tuple{uPriWest, uPriSouth, uPriBottom};
	}
	
	//! @brief Reconstructs the vertices for the directions East, North and Top
	template<typename AccessorType>
	auto computeENT(AccessorType& omAcc, const double Theta, size_t ix, size_t iy, size_t iz) {

		auto [deriv_x, deriv_y, deriv_z] = prepareDerivsOffset(omAcc, Theta, ix, iy, iz);
		double shift(0.);

		double uPriEast = omAcc[ix - 1][iy][iz] + (0.5 - shift) * deriv_x;
		double uPriNorth = omAcc[ix][iy - 1][iz] + 0.5 * deriv_y;
		double uPriTop = omAcc[ix][iy][iz - 1] + 0.5 * deriv_z;
		
		return std::tuple{uPriEast, uPriNorth, uPriTop};
	}

}

#endif