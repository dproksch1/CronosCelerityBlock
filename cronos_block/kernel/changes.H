#ifndef CHANGES_H
#define CHANGES_H = 1

//#include "save_data.H"
#include "fields_local.H"
#include "fluid.H"

namespace gpu {

	//! @brief Calculates the cell update modifiers for the given field
	template<typename FluxAccessorType, typename IdxConstType>
	double get_Changes1DCart(const FluxAccessorType& numfM, const FluxAccessorType &numfP,
			int dir, IdxConstType idx, int q) {

		return (numfP[q] - numfM[q])*idx[dir];
	}

	//! @brief Converts the numerical flux to cell update modifiers and applies them to the accumulation buffers
	template<typename NumAccessorType, typename NumFluxAccessorType, typename SinglePtotalAccessorType, typename IdxConstType>
	void get_Changes(NumAccessorType& nomAcc, int ix, int iy, int iz, int dir, const NumFluxAccessorType& numfM,
					const SinglePtotalAccessorType &numPtotal_M, const NumFluxAccessorType &numfP, const SinglePtotalAccessorType &numPtotal_P,
					int num_len, size_t z_max, IdxConstType idx) {
		
		int iPos;
		if(dir==0) {
			iPos = ix;
		} else if (dir==1) {
			iPos = iy;
		} else {
			iPos = iz;
		}
		
		for (int q = 0; q < num_len; ++q){

			int q_global = q;

			// Cartesian case
			nomAcc[ix][iy][iz].mat[q] += get_Changes1DCart(numfM, numfP, dir, idx, q);

			// Handle possibly extracted pressure
#if EXTRACT_PRESSURE == TRUE
			nom[q_global](ix, iy, iz) += get_Changes1DCart(gdata, numPtotal_P, numPtotal_M, dir, iPos, q);
#endif

		}
	}

}

#endif