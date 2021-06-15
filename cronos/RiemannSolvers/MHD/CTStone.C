#include <iomanip>
#include "CTStone.H"

CTStone::CTStone(const Data &gdata, int dir,
                 bool IntegrateA, int i_magFluid) : ConstrainedTransport(gdata, dir, IntegrateA, i_magFluid)
{
}

void CTStone::get_NumEmf2D(const Data &gdata,
                           const phys_fields_2D &pfLL,
                           const phys_fields_2D &pfLR,
                           const phys_fields_2D &pfRL,
                           const phys_fields_2D &pfRR,
                           fields_2D &fl) {
	
#if (FLUID_TYPE == CRONOS_MHD)
	this->eps_CentreCT = 1.e-6;
	int dir0, dir1;
	assert(dir >= 0 && dir < DIM);
  
	if(dir == 0) {
		dir0 = 1;
		dir1 = 2;
	} else if (dir == 1) {
		dir0 = 2;
		dir1 = 0;
	} else {
		dir0 = 0;
		dir1 = 1;
	}

	if(dir == 0) {

		for (int k = -1; k <= gdata.mx[2]+1; ++k){
			for (int j = -1; j <= gdata.mx[1]+1; ++j){
				fl.Ex(j,k) = -(fl.v_y(j,k)*0.5*(fl.B_z(j,k) + fl.B_z(j,k-1)) -
				               fl.v_z(j,k)*0.5*(fl.B_y(j,k) + fl.B_y(j-1,k)));
			}
		}

#if (STONE_TYPE == STONE_SIMPLE)
		for (int k = 0; k <= gdata.mx[2]+1; ++k){
			for (int j = 0; j <= gdata.mx[1]+1; ++j){
				fl.emf(j,k) = (0.5*((fl.fluxBT[q_By](j,k) + fl.fluxBT[q_By](j-1,k)) -
				                    (fl.fluxSN[q_Bz](j,k) + fl.fluxSN[q_Bz](j,k-1))) -
				               0.25*(fl.Ex(j-1,k-1) + fl.Ex(j  ,k-1) +
				                     fl.Ex(j-1,k  ) + fl.Ex(j  ,k  )));
			}
		}

#elif (STONE_TYPE == STONE_CENTRE) 

		for (int k = 0; k <= gdata.mx[2]+1; ++k){
			for (int j = 0; j <= gdata.mx[1]+1; ++j){

				// Equivalent to Gardiner & Stones derivative at j-1/4 where
				// -fl.fluxSN[q_Bz](j+1,k) is Ex(j+1/2,k)

				REAL dExdyjm1_4P = 2.*(fl.Ex(j,k  ) + fl.fluxSN[q_Bz](j,k  ));
				REAL dExdyjm1_4M = 2.*(fl.Ex(j,k-1) + fl.fluxSN[q_Bz](j,k-1));

				// for the different cases see, e.g., Eq (50) in Gardiner & Stone
				if(fl.fluxBT[0](j,k) < -eps_CentreCT) {
					fl.dExdyjm1_4(j,k) = dExdyjm1_4P;
				} else if(fl.fluxBT[0](j,k) > eps_CentreCT) {
					fl.dExdyjm1_4(j,k) = dExdyjm1_4M;
				} else {
					fl.dExdyjm1_4(j,k) = 0.5*(dExdyjm1_4M + dExdyjm1_4P);
				}

			}
		}


		for (int k = 0; k <= gdata.mx[2]+1; ++k){
			for (int j = 0; j <= gdata.mx[1]+1; ++j){

				// Equivalent to Gardiner & Stones derivative at j-1/4 where
				// -fl.fluxSN[q_Bz](j+1,k) is Ex(j+1/2,k)

				REAL dExdyjm3_4P = 2.*(-fl.fluxSN[q_Bz](j,k  ) - fl.Ex(j-1,k  ));
				REAL dExdyjm3_4M = 2.*(-fl.fluxSN[q_Bz](j,k-1) - fl.Ex(j-1,k-1));

				// for the different cases see, e.g., Eq (50) in Gardiner & Stone
				if(fl.fluxBT[0](j-1,k) < -eps_CentreCT) {
					fl.dExdyjm3_4(j,k) = dExdyjm3_4P;
				} else if(fl.fluxBT[0](j-1,k) > eps_CentreCT) {
					fl.dExdyjm3_4(j,k) = dExdyjm3_4M;
				} else {
					fl.dExdyjm3_4(j,k) = 0.5*(dExdyjm3_4M + dExdyjm3_4P);
				}

			}
		}

		for (int k = 0; k <= gdata.mx[2]+1; ++k){
			for (int j = 0; j <= gdata.mx[1]+1; ++j){

				// Equivalent to Gardiner & Stones derivative at k+1/4 where
				//  fluxBT[5](j,k+1) is Ex(j,k+1/2)
				REAL dExdzkm1_4P = 2.*( fl.Ex(j  ,k) - fl.fluxBT[5](j  ,k));
				REAL dExdzkm1_4M = 2.*( fl.Ex(j-1,k) - fl.fluxBT[5](j-1,k));

				// for different cases see, e.g., Eq (50) in Gardiner & Stone

				if(fl.fluxSN[0](j,k) < -eps_CentreCT) {
					fl.dExdzkm1_4(j,k) = dExdzkm1_4P;
				} else if(fl.fluxSN[0](j,k) > eps_CentreCT) {
					fl.dExdzkm1_4(j,k) = dExdzkm1_4M;
				} else {
					fl.dExdzkm1_4(j,k) = 0.5*(dExdzkm1_4M + dExdzkm1_4P);
				}

			}
		}

		for (int k = 0; k <= gdata.mx[2]+1; ++k){
			for (int j = 0; j <= gdata.mx[1]+1; ++j){

				// Equivalent to Gardiner & Stones derivative at j+1/4 where
				//  fluxBT[q_By](j,k+1) is Ex(j,k+1/2)
				REAL dExdzkm3_4P = 2.*( fl.fluxBT[q_By](j  ,k) - fl.Ex(j  ,k-1));
				REAL dExdzkm3_4M = 2.*( fl.fluxBT[q_By](j-1,k) - fl.Ex(j-1,k-1));

				// for different cases see, e.g., Eq (50) in Gardiner & Stone
				if(fl.fluxSN[0](j,k-1) < -eps_CentreCT) {
					fl.dExdzkm3_4(j,k) = dExdzkm3_4P;
				} else if(fl.fluxSN[0](j,k-1) > eps_CentreCT) {
					fl.dExdzkm3_4(j,k) = dExdzkm3_4M;
				} else {
					fl.dExdzkm3_4(j,k) = 0.5*(dExdzkm3_4M + dExdzkm3_4P);
				}

			}
		}




		for (int k = 0; k <= gdata.mx[2]+1; ++k){
			for (int j = 0; j <= gdata.mx[1]+1; ++j){
				fl.emf(j,k) = (0.25*((fl.fluxBT[q_By](j,k) + fl.fluxBT[q_By](j-1,k)) -
				                     (fl.fluxSN[q_Bz](j,k) + fl.fluxSN[q_Bz](j,k-1))) +
				               0.125*(fl.dExdyjm3_4(j,k) - fl.dExdyjm1_4(j,k) + 
				                      fl.dExdzkm3_4(j,k) - fl.dExdzkm1_4(j,k)));
			}
		}	


#endif

	} else if (dir == 1) {

		for (int i = -1; i <= gdata.mx[0]+1; ++i){
			for (int k = -1; k <= gdata.mx[2]+1; ++k){
				fl.Ey(k,i) = -(fl.v_z(k,i)*0.5*(fl.B_x(k,i) + fl.B_x(k,i-1)) - 
				               fl.v_x(k,i)*0.5*(fl.B_z(k,i) + fl.B_z(k-1,i)));
			}
		}

#if (STONE_TYPE == STONE_SIMPLE)

		for (int i = 0; i <= gdata.mx[0]+1; ++i){
			for (int k = 0; k <= gdata.mx[2]+1; ++k){
				fl.emf(k,i) = (0.5*((fl.fluxWE[q_Bz](k,i) + fl.fluxWE[q_Bz](k-1,i)) -
				                    (fl.fluxBT[q_Bx](k,i) + fl.fluxBT[q_Bx](k,i-1))) -
				               0.25*(fl.Ey(k-1,i-1) + fl.Ey(k  ,i-1) +
				                     fl.Ey(k-1,i  ) + fl.Ey(k  ,i  )));


			}
		}

	
#elif (STONE_TYPE == STONE_CENTRE) 

		for (int i = 0; i <= gdata.mx[0]+1; ++i){
			for (int k = 0; k <= gdata.mx[2]+1; ++k){

				// Equivalent to Gardiner & Stones derivative at i-1/4 where
				// fl.fluxWE[q_Bz](k,i+1) is Ey(i+1/2,k)
				REAL dEydxim1_4P = 2.*(fl.Ey(k  ,i) - fl.fluxWE[q_Bz](k  ,i));
				REAL dEydxim1_4M = 2.*(fl.Ey(k-1,i) - fl.fluxWE[q_Bz](k-1,i));

				// for different cases see, e.g., Eq (50) in Gardiner & Stone

				if(fl.fluxBT[0](k,i) < -eps_CentreCT) {
					fl.dEydxim1_4(k,i) = dEydxim1_4P;
				} else if(fl.fluxBT[0](k,i) > eps_CentreCT) {
					fl.dEydxim1_4(k,i) = dEydxim1_4M;
				} else {
					fl.dEydxim1_4(k,i) = 0.5*(dEydxim1_4M + dEydxim1_4P);
				}

			}
		}

		for (int i = 0; i <= gdata.mx[0]+1; ++i){
			for (int k = 0; k <= gdata.mx[2]+1; ++k){

				// Equivalent to Gardiner & Stones derivative at i-1/4 where
				// fl.fluxWE[q_Bz](k,i+1) is Ey(i+1/2,k)
				REAL dEydxim3_4P = 2.*(fl.fluxWE[q_Bz](k  ,i) - fl.Ey(k  ,i-1));
				REAL dEydxim3_4M = 2.*(fl.fluxWE[q_Bz](k-1,i) - fl.Ey(k-1,i-1));

				// for different cases see, e.g., Eq (50) in Gardiner & Stone
				if(fl.fluxBT[0](k,i-1) < -eps_CentreCT) {
					fl.dEydxim3_4(k,i) = dEydxim3_4P;
				} else if(fl.fluxBT[0](k,i-1) > eps_CentreCT) {
					fl.dEydxim3_4(k,i) = dEydxim3_4M;
				} else {
					fl.dEydxim3_4(k,i) = 0.5*(dEydxim3_4M + dEydxim3_4P);
				}

			}
		}

		for (int i = 0; i <= gdata.mx[0]+1; ++i){
			for (int k = 0; k <= gdata.mx[2]+1; ++k){

				// Equivalent to Gardiner & Stones derivative at j+1/4 where
				// -fl.fluxBT[q_Bx](k,i  ) is Ey(i,k+1/2)
				REAL dEydzkm1_4P = 2.*( fl.Ey(k,i  ) + fl.fluxBT[q_Bx](k,i  ));
				REAL dEydzkm1_4M = 2.*( fl.Ey(k,i-1) + fl.fluxBT[q_Bx](k,i-1));

				// for different cases see, e.g., Eq (50) in Gardiner & Stone

				if(fl.fluxWE[0](k,i) < -eps_CentreCT) {
					fl.dEydzkm1_4(k,i) = dEydzkm1_4P;
				} else if(fl.fluxWE[0](k,i) > eps_CentreCT) {
					fl.dEydzkm1_4(k,i) = dEydzkm1_4M;
				} else {
					fl.dEydzkm1_4(k,i) = 0.5*(dEydzkm1_4M + dEydzkm1_4P);
				}

			}
		}


		for (int i = 0; i <= gdata.mx[0]+1; ++i){
			for (int k = 0; k <= gdata.mx[2]+1; ++k){

				// Equivalent to Gardiner & Stones derivative at j+1/4 where
				//  fl.fluxSN[q_Bx](i,j+1) is Ez(i,j+1/2)
				REAL dEydzkm3_4P = 2.*(-fl.fluxBT[q_Bx](k,i  ) - fl.Ey(k-1,i  ));
				REAL dEydzkm3_4M = 2.*(-fl.fluxBT[q_Bx](k,i-1) - fl.Ey(k-1,i-1));

				// for different cases see, e.g., Eq (50) in Gardiner & Stone

				if(fl.fluxWE[0](k-1,i) < -eps_CentreCT) {
					fl.dEydzkm3_4(k,i) = dEydzkm3_4P;
				} else if(fl.fluxWE[0](k-1,i) > eps_CentreCT) {
					fl.dEydzkm3_4(k,i) = dEydzkm3_4M;
				} else {
					fl.dEydzkm3_4(k,i) = 0.5*(dEydzkm3_4M + dEydzkm3_4P);
				}

			}
		}


		for (int i = 0; i <= gdata.mx[0]+1; ++i){
			for (int k = 0; k <= gdata.mx[2]+1; ++k){
				fl.emf(k,i) = (0.25*((fl.fluxWE[q_Bz](k,i) + fl.fluxWE[q_Bz](k-1,i)) -
				                     (fl.fluxBT[q_Bx](k,i) + fl.fluxBT[q_Bx](k,i-1))) +
				               0.125*(fl.dEydxim3_4(k,i) - fl.dEydxim1_4(k,i) + 
				                      fl.dEydzkm3_4(k,i) - fl.dEydzkm1_4(k,i)));

			}
		}

#endif

	} else {

		for (int j = -1; j <= gdata.mx[1]+1; ++j){
			for (int i = -1; i <= gdata.mx[0]+1; ++i){
				fl.Ez(i,j) = -(fl.v_x(i,j)*0.5*(fl.B_y(i,j) + fl.B_y(i,j-1)) -
				               fl.v_y(i,j)*0.5*(fl.B_x(i,j) + fl.B_x(i-1,j)));
			}
		}

#if (STONE_TYPE == STONE_SIMPLE)

		for (int j = 0; j <= gdata.mx[1]+1; ++j){
			for (int i = 0; i <= gdata.mx[0]+1; ++i){
				fl.emf(i,j) = (0.5*((fl.fluxSN[q_Bx](i,j) + fl.fluxSN[q_Bx](i-1,j)) -
				                    (fl.fluxWE[q_By](i,j) + fl.fluxWE[q_By](i,j-1))) -
				               0.25*(fl.Ez(i-1,j-1) + fl.Ez(i  ,j-1) +
				                     fl.Ez(i-1,j  ) + fl.Ez(i  ,j  )));

			}
		}

#elif (STONE_TYPE == STONE_CENTRE) 

		for (int j = 0; j <= gdata.mx[1]+1; ++j){
			for (int i = 0; i <= gdata.mx[0]+1; ++i){

				// Equivalent to Gardiner & Stones derivative at i-1/4 where
				// -fl.fluxWE[q_By](i+1,j) is Ez(i+1/2,j)
				REAL dEzdxim1_4P = 2.*(fl.Ez(i,j  ) + fl.fluxWE[q_By](i,j  ));
				REAL dEzdxim1_4M = 2.*(fl.Ez(i,j-1) + fl.fluxWE[q_By](i,j-1));

				// for different cases see, e.g., Eq (50) in Gardiner & Stone

				if(fl.fluxSN[0](i,j) < -eps_CentreCT) {
					fl.dEzdxim1_4(i,j) = dEzdxim1_4P;
				} else if(fl.fluxSN[0](i,j) > eps_CentreCT) {
					fl.dEzdxim1_4(i,j) = dEzdxim1_4M;
				} else {
					fl.dEzdxim1_4(i,j) = 0.5*(dEzdxim1_4M + dEzdxim1_4P);
				}

			}
		}

		for (int j = 0; j <= gdata.mx[1]+1; ++j){
			for (int i = 0; i <= gdata.mx[0]+1; ++i){

				// Equivalent to Gardiner & Stones derivative at i-1/4 where
				// -fl.fluxWE[q_By](i+1,j) is Ez(i+1/2,j)
				REAL dEzdxim3_4P = 2.*(-fl.fluxWE[q_By](i,j  ) - fl.Ez(i-1,j  ));
				REAL dEzdxim3_4M = 2.*(-fl.fluxWE[q_By](i,j-1) - fl.Ez(i-1,j-1));

				// for different cases see, e.g., Eq (50) in Gardiner & Stone

				if(fl.fluxSN[0](i-1,j) < -eps_CentreCT) {
					fl.dEzdxim3_4(i,j) = dEzdxim3_4P;
				} else if(fl.fluxSN[0](i-1,j) > eps_CentreCT) {
					fl.dEzdxim3_4(i,j) = dEzdxim3_4M;
				} else {
					fl.dEzdxim3_4(i,j) = 0.5*(dEzdxim3_4M + dEzdxim3_4P);
				}

			}
		}


		for (int j = 0; j <= gdata.mx[1]+1; ++j){
			for (int i = 0; i <= gdata.mx[0]+1; ++i){

				// Equivalent to Gardiner & Stones derivative at j+1/4 where
				//  fl.fluxSN[q_Bx](i,j+1) is Ez(i,j+1/2)
				REAL dEzdyjm1_4P = 2.*( fl.Ez(i  ,j) - fl.fluxSN[q_Bx](i  ,j));
				REAL dEzdyjm1_4M = 2.*( fl.Ez(i-1,j) - fl.fluxSN[q_Bx](i-1,j));

				// for different cases see, e.g., Eq (50) in Gardiner & Stone

				if(fl.fluxWE[0](i,j) < -eps_CentreCT) {
					fl.dEzdyjm1_4(i,j) = dEzdyjm1_4P;
				} else if(fl.fluxWE[0](i,j) > eps_CentreCT) {
					fl.dEzdyjm1_4(i,j) = dEzdyjm1_4M;
				} else {
					fl.dEzdyjm1_4(i,j) = 0.5*(dEzdyjm1_4M + dEzdyjm1_4P);
				}

			}
		}


		for (int j = 0; j <= gdata.mx[1]+1; ++j){
			for (int i = 0; i <= gdata.mx[0]+1; ++i){

				// Equivalent to Gardiner & Stones derivative at j+1/4 where
				//  fl.fluxSN[q_Bx](i,j+1) is Ez(i,j+1/2)
				REAL dEzdyjm3_4P = 2.*( fl.fluxSN[q_Bx](i  ,j) - fl.Ez(i  ,j-1));
				REAL dEzdyjm3_4M = 2.*( fl.fluxSN[q_Bx](i-1,j) - fl.Ez(i-1,j-1));

				// for different cases see, e.g., Eq (50) in Gardiner & Stone

				if(fl.fluxWE[0](i,j-1) < -eps_CentreCT) {
					fl.dEzdyjm3_4(i,j) = dEzdyjm3_4P;
				} else if(fl.fluxWE[0](i,j-1) > eps_CentreCT) {
					fl.dEzdyjm3_4(i,j) = dEzdyjm3_4M;
				} else {
					fl.dEzdyjm3_4(i,j) = 0.5*(dEzdyjm3_4M + dEzdyjm3_4P);
				}

			}
		}




		for (int j = 0; j <= gdata.mx[1]+1; ++j){
			for (int i = 0; i <= gdata.mx[0]+1; ++i){
				fl.emf(i,j) = (0.25*((fl.fluxSN[q_Bx](i,j) + fl.fluxSN[q_Bx](i-1,j)) -
				                     (fl.fluxWE[q_By](i,j) + fl.fluxWE[q_By](i,j-1))) +
				               0.125*(fl.dEzdyjm3_4(i,j) - fl.dEzdyjm1_4(i,j) +
				                      fl.dEzdxim3_4(i,j) - fl.dEzdxim1_4(i,j)));

			}
		}
#endif

	}
#endif
}
