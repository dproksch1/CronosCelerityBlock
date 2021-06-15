#include "DissipationMHD.H"

void DissipationMHD::modifyEmf(Data &gdata, ProblemType &Problem, fields_2D &fl,
                               int dir, int layer) {
	assert(dir >= 0 && dir < DIM);
 
	if(dir == 0) {

#if (NON_LINEAR_GRID == CRONOS_OFF)
		REAL idy = gdata.idx[1];
		REAL idz = gdata.idx[2];
#endif

		int i(layer);
		for (int k = -1; k <= gdata.mx[2]+1; ++k){
#if (NON_LINEAR_GRID == CRONOS_ON)
			REAL idz = gdata.getCen_idx(2,k);
#endif
			for (int j = -1; j <= gdata.mx[1]+1; ++j){
#if (NON_LINEAR_GRID == CRONOS_ON)
				REAL idy = gdata.getCen_idx(1,j);
#endif

				REAL eta_val = Problem.eta(gdata, 1.*layer,j-0.5,k-0.5);

#if (GEOM != CARTESIAN)

#if (NON_LINEAR_GRID == CRONOS_OFF)
				REAL idA = 1./(gdata.h1(i,j-0.5,k-0.5)*gdata.h2(i,j-0.5,k-0.5));

				fl.emf(j,k) += eta_val*((fl.B_z(j  ,k-1)*gdata.h2(i,j  ,k-0.5) -
				                         fl.B_z(j-1,k-1)*gdata.h2(i,j-1,k-0.5))*idy -
				                        (fl.B_y(j-1,k  )*gdata.h1(i,j-0.5,k  ) -
				                         fl.B_y(j-1,k-1)*gdata.h1(i,j-0.5,k-1))*idz)*idA;

#else
				REAL idA = 1./(gdata.h1(i,j,k,0,-1,-1)*gdata.h2(i,j,k,0,-1,-1));

				fl.emf(j,k) += eta_val*((fl.B_z(j  ,k-1)*gdata.h2(i,j  ,k,0,0,-1) -
				                         fl.B_z(j-1,k-1)*gdata.h2(i,j-1,k,0,0,-1))*idy -
				                        (fl.B_y(j-1,k  )*gdata.h1(i,j,k  ,0,-1,0  ) -
				                         fl.B_y(j-1,k-1)*gdata.h1(i,j,k-1,0,-1,0))*idz)*idA;
	
				
#endif

#else

				fl.emf(j,k) += eta_val*((fl.B_z(j  ,k-1) -
				                         fl.B_z(j-1,k-1))*idy -
				                        (fl.B_y(j-1,k  ) -
				                         fl.B_y(j-1,k-1))*idz);

#endif

				}
		}
		
	} else if (dir == 1) {

#if (NON_LINEAR_GRID == CRONOS_OFF)
		REAL idx = gdata.idx[0];
		REAL idz = gdata.idx[2];
#endif

		// Beware dimensions 1 and 3 are swapped(!)
		REAL j(1.*layer);
		for (int k = -1; k <= gdata.mx[2]+1; ++k){
#if (NON_LINEAR_GRID == CRONOS_ON)
			REAL idz = gdata.getCen_idx(2,k);
#endif
			for (int i = -1; i <= gdata.mx[0]+1; ++i){
#if (NON_LINEAR_GRID == CRONOS_ON)
				REAL idx = gdata.getCen_idx(0,i);
#endif

				REAL eta_val = Problem.eta(gdata, i-0.5,1.*layer,k-0.5);

#if (GEOM != CARTESIAN)

#if (NON_LINEAR_GRID == CRONOS_OFF)
				REAL idA = 1./(gdata.h0(i-0.5,j,k-0.5)*gdata.h2(i-0.5,j,k-0.5));

				fl.emf(k,i) += eta_val*((fl.B_x(k  ,i-1)*gdata.h0(i-0.5,j,k  ) - 
				                         fl.B_x(k-1,i-1)*gdata.h0(i-0.5,j,k-1))*idz -
				                        (fl.B_z(k-1,i  )*gdata.h2(i  ,j,k-0.5) -
				                         fl.B_z(k-1,i-1)*gdata.h2(i-1,j,k-0.5))*idx)*idA;

#else
				REAL idA = 1./(gdata.h0(i,j,k,-1,0,-1)*
				               gdata.h2(i,j,k,-1,0,-1));

				fl.emf(k,i) += eta_val*((fl.B_x(k  ,i-1)*gdata.h0(i,j,k  ,-1,0,0) - 
				                         fl.B_x(k-1,i-1)*gdata.h0(i,j,k-1,-1,0,0))*idz -
				                        (fl.B_z(k-1,i  )*gdata.h2(i  ,j,k,0,0,-1) -
				                         fl.B_z(k-1,i-1)*gdata.h2(i-1,j,k,0,0,-1))*idx)*idA;		
				
#endif

#else

				fl.emf(k,i) += eta_val*((fl.B_x(k  ,i-1) - 
				                         fl.B_x(k-1,i-1))*idz -
				                        (fl.B_z(k-1,i  ) -
				                         fl.B_z(k-1,i-1))*idx);

#endif

			}
		}

	} else {
		
#if (NON_LINEAR_GRID == CRONOS_OFF)
		REAL idx = gdata.idx[0];
		REAL idy = gdata.idx[1];
#endif

		REAL k(1.*layer);
		for (int j = -1; j <= gdata.mx[1]+1; ++j){
#if (NON_LINEAR_GRID == CRONOS_ON)
			REAL idy = gdata.getCen_idx(1,j);
#endif
			for (int i = -1; i <= gdata.mx[0]+1; ++i){
#if (NON_LINEAR_GRID == CRONOS_ON)
				REAL idx = gdata.getCen_idx(0,i);
#endif

				REAL eta_val = Problem.eta(gdata, i-0.5,j-0.5,1.*layer);

#if (GEOM != CARTESIAN)

#if (NON_LINEAR_GRID == CRONOS_OFF)
				REAL idA = 1./(gdata.h0(i-0.5,j-0.5,k)*gdata.h1(i-0.5,j-0.5,k));

				fl.emf(i,j) += eta_val*((fl.B_y(i  ,j-1)*gdata.h1(i  ,j-0.5,k) - 
				                         fl.B_y(i-1,j-1)*gdata.h1(i-1,j-0.5,k))*idx -
				                        (fl.B_x(i-1,j  )*gdata.h0(i-0.5,j  ,k) -
				                         fl.B_x(i-1,j-1)*gdata.h0(i-0.5,j-1,k))*idy)*idA;
#else
				REAL idA = 1./(gdata.h0(i,j,k,-1,-1,0)*
				               gdata.h1(i,j,k,-1,-1,0));

				fl.emf(i,j) += eta_val*((fl.B_y(i  ,j-1)*gdata.h1(i  ,j,k,0,-1,0) - 
				                         fl.B_y(i-1,j-1)*gdata.h1(i-1,j,k,0,-1,0))*idx -
				                        (fl.B_x(i-1,j  )*gdata.h0(i,j  ,k,-1,0,0) -
				                         fl.B_x(i-1,j-1)*gdata.h0(i,j-1,k,-1,0,0))*idy)*idA;		

#endif

#else

				fl.emf(i,j) += eta_val*((fl.B_y(i  ,j-1) - 
				                         fl.B_y(i-1,j-1))*idx -
				                        (fl.B_x(i-1,j  ) -
				                         fl.B_x(i-1,j-1))*idy);

#endif

			}
		}
	}

}
