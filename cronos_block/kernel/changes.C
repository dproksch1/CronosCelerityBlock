#include "changes.H"
#include <iomanip>


void get_Changes(const Data &gdata,
		const num_fields_0D &numfM, const num_fields_0D &numfP, NumMatrix<REAL, 3> nom[],
		int ix, int iy, int iz, int dir, int Fluid_Type, int iFluid) {

	// Length of first dimension contains number of fields
	int num = numfP.flux_num.getLength();
	int iPos;
	if(dir==0) {
		iPos = ix;
	} else if (dir==1) {
		iPos = iy;
	} else {
		iPos = iz;
	}

#if (GEOM != CARTESIAN)
	REAL AreaM, AreaP;
	if(dir==0) {
		AreaM = gdata.get_CellArea_x(ix,iy,iz);
		AreaP = gdata.get_CellArea_x(ix+1,iy,iz);
	} else if (dir==1) {
		AreaM = gdata.get_CellArea_y(ix,iy,iz);
		AreaP = gdata.get_CellArea_y(ix,iy+1,iz);
	} else {
		AreaM = gdata.get_CellArea_z(ix,iy,iz);
		AreaP = gdata.get_CellArea_z(ix,iy,iz+1);
	}
	REAL iVol = 1./gdata.get_CellVolume(ix,iy,iz);
#endif


	for (int q = 0; q < num; ++q){
		
		int q_global = q;

#if (GEOM == CARTESIAN)
		// Cartesian case
		nom[q_global](ix, iy, iz) += get_Changes1DCart(gdata, numfP.flux_num, numfM.flux_num, dir, iPos, q);
#else

		// general case
		double change = (numfP.flux_num[q]*AreaP - numfM.flux_num[q]*AreaM)*iVol;

		nom[q_global](ix, iy, iz) += change;

		// Compute mass flux in case of spherical coordinates
#if (GEOM == SPHERICAL)
		if(dir==0) {
			gdata.massFlux(ix) += numfM.flux_num[0]*AreaM;
		}
#endif

#endif // IF GEOM

		// Handle possibly extracted pressure
#if EXTRACT_PRESSURE == TRUE
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

	}

}

REAL get_Changes1DCart(const Data &gdata,
		const NumArray<REAL> &fluxP, const NumArray<REAL> &fluxM,
		int dir, int iPos, int q) {

#if (NON_LINEAR_GRID == CRONOS_ON)
	REAL idx = gdata.getCen_idx(dir, iPos);
#else
	REAL idx = gdata.idx[dir];
#endif

	return (fluxP[q] - fluxM[q])*idx;

}
