#ifndef CHANGES_H
#define CHANGES_H = 1

//#include "save_data.H"
#include "fields_local.H"
#include "fluid.H"

namespace gpu {

	inline size_t nom_iyiz(size_t iy, size_t iz, size_t z_max) {
		return (iy * z_max) + iz;
	}

	template<typename NumFluxAccessorType, typename IdxConstType>
	double get_Changes1DCart(const NumFluxAccessorType& numfM, const NumFluxAccessorType &numfP,
			int dir, IdxConstType idx, int q) {

/*#if (NON_LINEAR_GRID == CRONOS_ON)
		//double idx = gdata.getCen_idx(dir, iPos);
#else
		double idx_dir = idx[dir];
#endif*/

		return (numfP[q] - numfM[q])*idx[dir];

	}

	template<typename NumAccessorType, typename NumFluxAccessorType, typename IdxConstType>
	void apply_change(NumAccessorType& nomAcc, int q, int ix, int iy, int iz, int dir, const NumFluxAccessorType& numfM, const NumFluxAccessorType &numfP, int num_len, IdxConstType idx) {
		switch(q) {
			case 0: nomAcc[ix][iy][iz].rho  += get_Changes1DCart(numfM, numfP, dir, idx, 0); break;
			case 1: nomAcc[ix][iy][iz].sx   += get_Changes1DCart(numfM, numfP, dir, idx, 1); break;
			case 2: nomAcc[ix][iy][iz].sy   += get_Changes1DCart(numfM, numfP, dir, idx, 2); break;
			case 3: nomAcc[ix][iy][iz].sz   += get_Changes1DCart(numfM, numfP, dir, idx, 3); break;
			case 4: nomAcc[ix][iy][iz].Eges += get_Changes1DCart(numfM, numfP, dir, idx, 4); break;
			default: nomAcc[ix][iy][iz].rho += get_Changes1DCart(numfM, numfP, dir, idx, 0);
		}
	}


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

#if (GEOM == CARTESIAN)
			// Cartesian case
			//nomAcc[ix][nom_iyiz(iy,iz,z_max)][q_global] += get_Changes1DCart(numfM, numfP, dir, idx, q);
			apply_change(nomAcc, q, ix, iy, iz, dir, numfM, numfP, num_len, idx);
#else

			/*// general case
			double change = (numfP[q]*AreaP - numfM[q]*AreaM)*iVol;

			nom[ix][iy][iz] += change;*/
			nomAcc[ix][nom_iyiz(iy,iz,z_max)][q_global] += get_Changes1DCart(numfM, numfP, dir, iPos, q);

			// Compute mass flux in case of spherical coordinates
#if (GEOM == SPHERICAL)
			if(dir==0) {
				gdata.massFlux(ix) += numfM.flux_num[0]*AreaM;
			}
#endif

#endif // IF GEOM

			// Handle possibly extracted pressure (unimplemented for now)
/*#if EXTRACT_PRESSURE == TRUE
			// Get corresponding indices
#if(FLUID_TYPE == CRONOS_MULTIFLUID)
			int q_sx_global = gdata.fluids->fluids[iFluid].get_q_sx_global();
			int q_sy_global = gdata.fluids->fluids[iFluid].get_q_sy_global();
#else
			int q_sx_global = q_sx;
			int q_sy_global = q_sy;
#endif

#if (GEOM == CARTESIAN)

			nom[q_global](ix, iy, iz) += get_Changes1DCart(gdata, numfP.ptotal_num, numfM.ptotal_num, dir, iPos, q);

#else

			// general case
			double change = (numfP.ptotal_num*AreaP - numfM.ptotal_num*AreaM)*iVol;

			nom[q_global](ix, iy, iz) += change;

#endif // IF GEOM
#endif // IF EXTRACT PRESSURE
*/

		}
	}

	template<typename NumAccessorType, typename NumFluxAccessorType, typename SinglePtotalAccessorType, typename IdxConstType>
	void get_Changes(NumAccessorType& nomAcc, int q, int ix, int iy, int iz, int dir, const NumFluxAccessorType& numfM,
					const SinglePtotalAccessorType &numPtotal_M, const NumFluxAccessorType &numfP, const SinglePtotalAccessorType &numPtotal_P,
					int num_len, size_t z_max, IdxConstType idx) {
		

		int q_global = q;

#if (GEOM == CARTESIAN)
		// Cartesian case
		nomAcc[ix][iy][iz] += get_Changes1DCart(numfM, numfP, dir, idx, q);
#else

		/*// general case
		double change = (numfP[q]*AreaP - numfM[q]*AreaM)*iVol;

		nom[ix][iy][iz] += change;*/
		nomAcc[ix][nom_iyiz(iy,iz,z_max)][q_global] += get_Changes1DCart(numfM, numfP, dir, iPos, q);

		// Compute mass flux in case of spherical coordinates
#if (GEOM == SPHERICAL)
		if(dir==0) {
			gdata.massFlux(ix) += numfM.flux_num[0]*AreaM;
		}
#endif

#endif 

	}

	template<typename NumAccessorType, typename NumFluxAccessorType, typename IdxConstType>
	void get_Changes_Local(NumAccessorType& nomAcc, int q, int ix, int iy, int iz, int dir,
							const NumFluxAccessorType& numfM, const NumFluxAccessorType &numfP, IdxConstType idx) {
		

#if (GEOM == CARTESIAN)
		// Cartesian case
		nomAcc[ix][iy][iz] += get_Changes1DCart(numfM, numfP, dir, idx, q);
#endif // IF GEOM

	}
}

void get_Changes(const Data &, const num_fields_0D &, const num_fields_0D &,
			NumMatrix<double,3> [], int, int, int, int dir, int, int iFluid=0);

double get_Changes1DCart(const Data &, const NumArray<double> &, const NumArray<double> &,
			int, int, int);

#endif