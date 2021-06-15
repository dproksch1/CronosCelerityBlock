#include "gridgen.H"
#include <iomanip>


void HyperbolicSolver::get_Changes(Data &gdata,
                                   NumMatrix<REAL, 1> flux[],
                                   NumMatrix<REAL, 3> nom[],
                                   NumMatrix<REAL, 1> &ptotal,
                                   int dir, int j, int k, int num,
                                   int Fluid_Type, int iFluid) {

#if GEOM == CARTESIAN
	for (int q = 0; q < num; ++q){ 
#if(FLUID_TYPE == CRONOS_MULTIFLUID)
		int q_global = gdata.fluids->fluids[iFluid].get_IndexGlobal(q);
#else
		int q_global = q;
#endif
//		if((Fluid_Type == CRONOS_MHD && (q < 4 || q > 6)) ||
//		   Fluid_Type != CRONOS_MHD) {
		if((Fluid_Type != CRONOS_MHD) || (q!=q_Bx && q!=q_By && q!=q_Bz)) {
#if (NON_LINEAR_GRID == CRONOS_OFF)
			REAL idx = gdata.idx[dir];
#endif

			if(dir == 0) {
				for (int i = 0; i <= gdata.mx[dir]; ++i){
#if (NON_LINEAR_GRID == CRONOS_ON)
					REAL idx = gdata.getCen_idx(dir, i);
#endif
					nom[q_global](i,j,k) += (flux[q](i+1) - flux[q](i  ))*idx;
				}
			} else if (dir == 1) {
				for (int i = 0; i <= gdata.mx[dir]; ++i){
#if (NON_LINEAR_GRID == CRONOS_ON)
					REAL idx = gdata.getCen_idx(dir, i);
#endif
					nom[q_global](j,i,k) += (flux[q](i+1) - flux[q](i  ))*idx;
				}
			} else {
				for (int i = 0; i <= gdata.mx[dir]; ++i){
#if (NON_LINEAR_GRID == CRONOS_ON)
					REAL idx = gdata.getCen_idx(dir, i);
#endif
					nom[q_global](j,k,i) += (flux[q](i+1) - flux[q](i  ))*idx;
				}
			}
		}
	}
#else
	// Compute changes including geometrical factors -- be aware that
	// gflux is not a flux anymore (due to the fact that an area is
	// included
	for (int q = 0; q < num; ++q){ 
#if(FLUID_TYPE == CRONOS_MULTIFLUID)
		int q_global = gdata.fluids->fluids[iFluid].get_IndexGlobal(q);
#else
		int q_global = q;
#endif

//		if((Fluid_Type == CRONOS_MHD && (q < 4 || q > 6)) ||
//		   Fluid_Type != CRONOS_MHD) {
		if((Fluid_Type != CRONOS_MHD) || (q!=q_Bx && q!=q_By && q!=q_Bz)) {

#if (NON_LINEAR_GRID == CRONOS_OFF)
			REAL idx = gdata.idx[dir];
#endif

			if(dir == 0) {

				for (int i = 0; i <= gdata.mx[dir]+1; ++i){

					// Get geometrical factors
#if (NON_LINEAR_GRID == CRONOS_OFF)
					REAL f_geom_x = gdata.h1(i-0.5,j,k)*gdata.h2(i-0.5,j,k);
#else
					REAL f_geom_x = gdata.h1(i,j,k,-1,0,0)*gdata.h2(i,j,k,-1,0,0);
#endif
					gflux[dir](i) = flux[q](i)*f_geom_x;
				}

				for (int i = 0; i <= gdata.mx[dir]; ++i){
#if (NON_LINEAR_GRID == CRONOS_ON)
					REAL idx = gdata.getCen_idx(dir, i);
#endif
					REAL iVol = idx/gdata.get_CellGeomTrafo(i,j,k);
					nom[q_global](i,j,k) += (gflux[dir](i+1) - gflux[dir](i  ))*iVol;
				}

				// Comute mass flux for spherical geometry:
#if (GEOM == SPHERICAL)
				if(q==0) {
					for (int i = 0; i <= gdata.mx[dir]+1; ++i){
						gdata.massFlux(i) += gflux[dir](i)*gdata.getCen_dx(1,j)*gdata.getCen_dx(2,k);
					}
				}
#endif

			} else if (dir == 1) {

				for (int i = 0; i <= gdata.mx[dir]+1; ++i){
#if (NON_LINEAR_GRID == CRONOS_OFF)
					REAL f_geom_y = gdata.h0(j,i-0.5,k)*gdata.h2(j,i-0.5,k);
#else
					REAL f_geom_y = gdata.h0(j,i,k,0,-1,0)*gdata.h2(j,i,k,0,-1,0);
#endif
					gflux[dir](i) = flux[q](i)*f_geom_y;
				}

				for (int i = 0; i <= gdata.mx[dir]; ++i){
#if (NON_LINEAR_GRID == CRONOS_ON)
					REAL idx = gdata.getCen_idx(dir, i);
#endif
					REAL iVol = idx/gdata.get_CellGeomTrafo(j,i,k);
					nom[q_global](j,i,k) += (gflux[dir](i+1) -gflux[dir](i  ))*iVol;

				}

			} else {

				for (int i = 0; i <= gdata.mx[dir]+1; ++i){
#if (NON_LINEAR_GRID == CRONOS_OFF)
					REAL f_geom_z = gdata.h0(j,k,i-0.5)*gdata.h1(j,k,i-0.5);
#else
					REAL f_geom_z = gdata.h0(j,k,i,0,0,-1)*gdata.h1(j,k,i,0,0,-1);
#endif
					gflux[dir](i) = flux[q](i)*f_geom_z;
				}

				for (int i = 0; i <= gdata.mx[dir]; ++i){
#if (NON_LINEAR_GRID == CRONOS_ON)
					REAL idx = gdata.getCen_idx(dir, i);
#endif
					REAL iVol = idx/gdata.get_CellGeomTrafo(j,k,i);
					nom[q_global](j,k,i) += (gflux[dir](i+1) - gflux[dir](i  ))*iVol;
				}

			}
		}
	}
#endif


#if EXTRACT_PRESSURE == TRUE

	if(Fluid_Type != CRONOS_USER) {

#if GEOM == CARTESIAN

#if (NON_LINEAR_GRID == CRONOS_OFF)
		REAL idx = gdata.idx[dir];
#endif

		if(dir == 0) {
#if(FLUID_TYPE == CRONOS_MULTIFLUID)
			int q_sx_global = gdata.fluids->fluids[iFluid].get_q_sx_global();
#else
			int q_sx_global = q_sx;
#endif
			for (int i = 0; i <= gdata.mx[dir]; ++i){
#if (NON_LINEAR_GRID == CRONOS_ON)
				REAL idx = gdata.getCen_idx(dir, i);
#endif
				nom[q_sx_global](i,j,k) += (ptotal(i+1) - ptotal(i  ))*idx;
			}
		} else if (dir == 1) {
#if(FLUID_TYPE == CRONOS_MULTIFLUID)
			int q_sy_global = gdata.fluids->fluids[iFluid].get_q_sy_global();
#else
			int q_sy_global = q_sy;
#endif
			for (int i = 0; i <= gdata.mx[dir]; ++i){		
#if (NON_LINEAR_GRID == CRONOS_ON)
				REAL idx = gdata.getCen_idx(dir, i);
#endif
				nom[q_sy_global](j,i,k) += (ptotal(i+1) - ptotal(i  ))*idx;
			}
		} else {
#if(FLUID_TYPE == CRONOS_MULTIFLUID)
			int q_sz_global = gdata.fluids->fluids[iFluid].get_q_sz_global();
#else
			int q_sz_global = q_sz;
#endif
			for (int i = 0; i <= gdata.mx[dir]; ++i){
#if (NON_LINEAR_GRID == CRONOS_ON)
				REAL idx = gdata.getCen_idx(dir, i);
#endif
				nom[q_sz_global](j,k,i) += (ptotal(i+1) - ptotal(i  ))*idx;
			}
		}

#else

#if (NON_LINEAR_GRID == CRONOS_OFF)
		REAL idx = gdata.idx[dir];
#endif
	
		if(dir == 0) {

#if(FLUID_TYPE == CRONOS_MULTIFLUID)
			int q_sx_global = gdata.fluids->fluids[iFluid].get_q_sx_global();
#else
			int q_sx_global = q_sx;
#endif

			NumMatrix<double,1> gptotal(ptotal);
			for (int i = 0; i <= gdata.mx[dir]+1; ++i){
#if (NON_LINEAR_GRID == CRONOS_OFF)
				REAL f_geom_x = gdata.h1(i-0.5,j,k)*gdata.h2(i-0.5,j,k);
#else
				REAL f_geom_x = gdata.h1(i,j,k,-1,0,0)*gdata.h2(i,j,k,-1,0,0);
#endif
				gptotal(i) = ptotal(i)*f_geom_x;
			}
			
			for (int i = 0; i <= gdata.mx[dir]; ++i){
#if (NON_LINEAR_GRID == CRONOS_ON)
				REAL idx = gdata.getCen_idx(dir, i);
#endif
				REAL iVol = idx/gdata.get_CellGeomTrafo(i,j,k);
				nom[q_sx_global](i,j,k) += (gptotal(i+1) - gptotal(i  ))*iVol;
			}
		} else if (dir == 1) {
			
#if(FLUID_TYPE == CRONOS_MULTIFLUID)
			int q_sy_global = gdata.fluids->fluids[iFluid].get_q_sy_global();
#else
			int q_sy_global = q_sy;
#endif

			NumMatrix<double,1> gptotal(ptotal);
			for (int i = 0; i <= gdata.mx[dir]+1; ++i){
#if (NON_LINEAR_GRID == CRONOS_OFF)
				REAL f_geom_y = gdata.h0(j,i-0.5,k)*gdata.h2(j,i-0.5,k);
#else
				REAL f_geom_y = gdata.h0(j,i,k,0,-1,0)*gdata.h2(j,i,k,0,-1,0);
#endif
				gptotal(i) = ptotal(i)*f_geom_y;
			}
			for (int i = 0; i <= gdata.mx[dir]; ++i){
#if (NON_LINEAR_GRID == CRONOS_ON)
				REAL idx = gdata.getCen_idx(dir, i);
#endif
				REAL iVol = idx/gdata.get_CellGeomTrafo(j,i,k);
				nom[q_sy_global](j,i,k) += (gptotal(i+1) - gptotal(i  ))*iVol;
			}
			
		} else {

#if(FLUID_TYPE == CRONOS_MULTIFLUID)
			int q_sz_global = gdata.fluids->fluids[iFluid].get_q_sz_global();
#else
			int q_sz_global = q_sz;
#endif

			NumMatrix<double,1> gptotal(ptotal);
			for (int i = 0; i <= gdata.mx[dir]+1; ++i){
#if (NON_LINEAR_GRID == CRONOS_OFF)
				REAL f_geom_z = gdata.h0(j,k,i-0.5)*gdata.h1(j,k,i-0.5);
#else
				REAL f_geom_z = gdata.h0(j,k,i,0,0,-1)*gdata.h1(j,k,i,0,0,-1);
#endif
				gptotal(i) = ptotal(i)*f_geom_z;
			}
			
			for (int i = 0; i <= gdata.mx[dir]; ++i){
#if (NON_LINEAR_GRID == CRONOS_ON)
				REAL idx = gdata.getCen_idx(dir, i);
#endif
				REAL iVol = idx/gdata.get_CellGeomTrafo(j,k,i);
				nom[q_sz_global](j,k,i) += (gptotal(i+1) - gptotal(i  ))*iVol;
			}
		}
	
#endif

	}

#endif

}

void HyperbolicSolver::get_Changes(Data &gdata,
		NumArray<REAL> &fluxM_x, NumArray<REAL> &fluxP_x,
		NumArray<REAL> &fluxM_y, NumArray<REAL> &fluxP_y,
		NumArray<REAL> &fluxM_z, NumArray<REAL> &fluxP_z,
		NumMatrix<REAL, 3> nom[], double ptotal,
		int ix, int iy, int iz, int Fluid_Type, int iFluid) {

	// Length of first dimension contains number of fields
	int num = fluxP_x.getLength();

#if (GEOM != CARTESIAN)
	REAL AreaM_x = gdata.get_CellArea_x(ix,iy,iz);
	REAL AreaP_x = gdata.get_CellArea_x(ix+1,iy,iz);
	REAL AreaM_y = gdata.get_CellArea_y(ix,iy,iz);
	REAL AreaP_y = gdata.get_CellArea_y(ix,iy+1,iz);
	REAL AreaM_z = gdata.get_CellArea_z(ix,iy,iz);
	REAL AreaP_z = gdata.get_CellArea_z(ix,iy,iz+1);
	REAL iVol = 1./gdata.get_CellVolume(ix,iy,iz);
#endif


	for (int q = 0; q < num; ++q){
#if(FLUID_TYPE == CRONOS_MULTIFLUID)
		int q_global = gdata.fluids->fluids[iFluid].get_IndexGlobal(q);
#else
		int q_global = q;
#endif

		// Cartesian case
#if (GEOM == CARTESIAN)

		// Flux in x-dimension
		nom[q_global](ix, iy, iz) = get_Changes1DCart(gdata, fluxP_x, fluxM_x, 0, ix, q);
		// Flux in y-dimension
		nom[q_global](ix, iy, iz) += get_Changes1DCart(gdata, fluxP_y, fluxM_y, 1, iy, q);
		// Flux in z-dimension
		nom[q_global](ix, iy, iz) += get_Changes1DCart(gdata, fluxP_z, fluxM_z, 2, iz, q);

#else

		// general case
		double change = (fluxP_x[q]*AreaP_x - fluxM_x[q]*AreaM_x);
		change += fluxP_y[q]*AreaP_y - fluxM_y[q]*AreaM_y;
		change += fluxP_z[q]*AreaP_z - fluxM_z[q]*AreaM_z;
		change *= iVol;

		nom[q_global](ix, iy, iz) += change;

		// Compute mass flux in case of spherical coordinates
#if (GEOM == SPHERICAL)
		gdata.massFlux(ix) += fluxM_x[0]*AreaM_x;
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

		nom[q_global](ix, iy, iz) += get_Changes1DCart(gdata, ptotal(1), ptotal(0), 0, ix, q);
		nom[q_global](ix, iy, iz) += get_Changes1DCart(gdata, ptotal(3), ptotal(2), 1, iy, q);
		nom[q_global](ix, iy, iz) += get_Changes1DCart(gdata, ptotal(5), ptotal(3), 2, iz, q);

#else

		// general case
		double change = (ptotal(1)*AreaP_x - ptotal(0)*AreaM_x);
		change = (ptotal(3)*AreaP_y - ptotal(2)*AreaM_y);
		change = (ptotal(5)*AreaP_z - ptotal(4)*AreaM_z);
		change *= iVol;

		nom[q_global](ix, iy, iz) += change;

#endif // IF GEOM


#endif // IF EXTRACT PRESSURE

	}

}



void HyperbolicSolver::get_Changes(const Data &gdata,
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
#if(FLUID_TYPE == CRONOS_MULTIFLUID)
		int q_global = gdata.fluids->fluids[iFluid].get_IndexGlobal(q);
#else
		int q_global = q;
#endif

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

REAL HyperbolicSolver::get_Changes1DCart(const Data &gdata,
		const NumArray<REAL> &fluxP, const NumArray<REAL> &fluxM,
		int dir, int iPos, int q) {

#if (NON_LINEAR_GRID == CRONOS_ON)
	REAL idx = gdata.getCen_idx(dir, iPos);
#else
	REAL idx = gdata.idx[dir];
#endif

	return (fluxP[q] - fluxM[q])*idx;

}

void HyperbolicSolver::get_ChangesX(Data &gdata, NumMatrix<REAL, 1> &fluxP, NumMatrix<REAL,1> &fluxM,
		NumMatrix<REAL, 3> nom[], NumMatrix<REAL, 1> &ptotal,
		int i, int j, int k, int Fluid_Type, int iFluid) {

}

void HyperbolicSolver::get_Changes(Data &gdata, NumMatrix<REAL, 1> &fluxP, NumMatrix<REAL,1> &fluxM,
		NumMatrix<REAL, 3> nom[], NumMatrix<REAL, 1> &ptotal,
		int i, int j, int k, int Fluid_Type, int iFluid) {

	// Length of first dimension contains number of fields
	int num = fluxM.getHigh(0);

	for(int dir=0; dir<DIM; ++dir) {

		for (int q = 0; q < num; ++q){
#if(FLUID_TYPE == CRONOS_MULTIFLUID)
			int q_global = gdata.fluids->fluids[iFluid].get_IndexGlobal(q);
#else
			int q_global = q;
#endif

			if((Fluid_Type != CRONOS_MHD) || (q!=q_Bx && q!=q_By && q!=q_Bz)) {

#if GEOM == CARTESIAN

#if (NON_LINEAR_GRID == CRONOS_OFF)
				REAL idx = gdata.idx[dir];
#else
				REAL idx = gdata.getCen_idx(dir, i);
#endif
				nom[q_global](i,j,k) +=  (fluxP(q) - fluxM(q))*idx;


#else // non-cartesian grids

//				nom[q_global](i,j,k) +=


#endif

			} // Check if magnetic field
		} // q
	} // dir
}


void HyperbolicSolver::get_ChangesEmf(Data &gdata,
                                      NumMatrix<REAL, 3> nom [],
                                      NumMatrix<REAL, 2> &emf,
                                      int dir, int layer) {

	if(IntegrateA) {

		if(layer >= 0 && layer <= gdata.mx[dir]) {

			if(dir == 0) {

				for (int k = 0; k <= gdata.mx[2]+1; ++k){
					for (int j = 0; j <= gdata.mx[1]+1; ++j){
						
						nom[q_Bx](layer,j,k) = emf(j,k);
						
					}
				}
				
			} else if (dir == 1) {

				// Beware dimensions 1 and 3 are swapped(!)
				for (int k = 0; k <= gdata.mx[2]+1; ++k){
					for (int i = 0; i <= gdata.mx[0]+1; ++i){
						
						nom[q_By](i,layer,k) = emf(k,i);
						
					}
				}
				
			} else {
			
				for (int j = 0; j <= gdata.mx[1]+1; ++j){
					for (int i = 0; i <= gdata.mx[0]+1; ++i){
						
						nom[q_Bz](i,j,layer) = emf(i,j);

					}
				}

			}

		}

	} else {

		if(layer >= -1 && layer <= gdata.mx[dir]) {
			
#if (NON_LINEAR_GRID == CRONOS_OFF)
			REAL idx = gdata.idx[0];
			REAL idy = gdata.idx[1];
			REAL idz = gdata.idx[2];
#endif

			if(dir == 0) {

				// Emf = Ex in this case:
				for (int k = 0; k <= gdata.mx[2]; ++k){
					for (int j = 0; j <= gdata.mx[1]; ++j){
#if (NON_LINEAR_GRID == CRONOS_ON)
						REAL idz = gdata.getCen_idx(2, k);
#endif
						nom[q_By](layer,j,k) +=  (emf(j+1,k+1) - emf(j+1,k  ))*idz;
						
					}
				}

				for (int k = 0; k <= gdata.mx[2]; ++k){
					for (int j = 0; j <= gdata.mx[1]; ++j){
#if (NON_LINEAR_GRID == CRONOS_ON)
						REAL idy = gdata.getCen_idx(1, j);
#endif
						nom[q_Bz](layer,j,k) += -(emf(j+1,k+1) - emf(j  ,k+1))*idy;

					}
				}

			} else if (dir == 1) {
				
				// Emf = Ey in this case:

				// Beware dimensions 1 and 3 are swapped(!)
				for (int k = 0; k <= gdata.mx[2]; ++k){
					for (int i = 0; i <= gdata.mx[0]; ++i){
#if (NON_LINEAR_GRID == CRONOS_ON)
						REAL idz = gdata.getCen_idx(2, k);
#endif						
						nom[q_Bx](i,layer,k) += -(emf(k+1,i+1) - emf(k  ,i+1))*idz;
						
					}
				}

				// Beware dimensions 1 and 3 are swapped(!)
				for (int k = 0; k <= gdata.mx[2]; ++k){
					for (int i = 0; i <= gdata.mx[0]; ++i){
#if (NON_LINEAR_GRID == CRONOS_ON)
						REAL idx = gdata.getCen_idx(0, i);
#endif	
						nom[q_Bz](i,layer,k) +=  (emf(k+1,i+1) - emf(k+1,i  ))*idx;
						
					}
				}

			} else {

				// Emf = Ez in this case

				for (int j = 0; j <= gdata.mx[1]; ++j){
					for (int i = 0; i <= gdata.mx[0]; ++i){
#if (NON_LINEAR_GRID == CRONOS_ON)
						REAL idy = gdata.getCen_idx(1, j);
#endif							
						nom[q_Bx](i,j,layer) +=  (emf(i+1,j+1) - emf(i+1,j  ))*idy;

					}
				}

				for (int j = 0; j <= gdata.mx[1]; ++j){
					for (int i = 0; i <= gdata.mx[0]; ++i){
#if (NON_LINEAR_GRID == CRONOS_ON)
						REAL idx = gdata.getCen_idx(0, i);
#endif						
						nom[q_By](i,j,layer) += -(emf(i+1,j+1) - emf(i  ,j+1))*idx;

					}
				}

			}
		}

	}
}
