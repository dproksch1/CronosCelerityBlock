#include "RiemannSolverMHD.H"
#include <iomanip>

HLLSolver::HLLSolver(const Data &gdata, int dir, int Fluid_Type) : RiemannSolverMHD(gdata, dir, Fluid_Type) {
	veps = 1.e-120;
	norm = 1./sqrt(2.);
}

void HLLSolver::get_NumFlux(Queue /*queue*/, Data &gdata,
                            phys_fields_1D &pfL,
                            phys_fields_1D &pfR,
                            fields_1D &fl, const int &dir, int iFluid) const {
	for(int q=0; q<fl.get_num(); ++q) {
		for (int i = -1; i <= gdata.mx[dir]+1; ++i){
      
			REAL fac = 1./(fl.v_ch_pORIG(i)+fl.v_ch_mORIG(i)+veps);
      
			fl.fluxORIG[q](i) = (fl.v_ch_mORIG(i)*pfL.fluxORIG[q](i) +
			                 fl.v_ch_pORIG(i)*pfR.fluxORIG[q](i-1) -
			                 fl.v_ch_mORIG(i)*fl.v_ch_pORIG(i)*(pfL.uConORIG[q](i  )-
			                                            pfR.uConORIG[q](i-1)))*fac;
		}
	}

#if EXTRACT_PRESSURE == TRUE
	if(fl.isgeneric() == true) {
		for (int i = -1; i <= gdata.mx[dir]+1; ++i){
			
			REAL fac = 1./(fl.v_ch_p(i)+fl.v_ch_m(i)+veps);
			
#if(FLUID_TYPE==CRONOS_MULTIFLUID)
			fl.ptotal[iFluid](i) = (fl.v_ch_m(i)*pfL.ptotal(i) +
					fl.v_ch_p(i)*pfR.ptotal(i-1))*fac;
#else
			fl.ptotal(i) = (fl.v_ch_m(i)*pfL.ptotal(i) +
			                fl.v_ch_p(i)*pfR.ptotal(i-1))*fac;
#endif
		}
	}

#endif

}


void HLLSolver::get_NumFlux(const Data &gdata, const phys_fields_0D &pfM,
	const phys_fields_0D &pfP, num_fields_0D &f_num, int dir, int iFluid) const {
	for(int q=0; q<pfM.get_num(); ++q) {

		REAL fac = 1./(f_num.v_ch_p+f_num.v_ch_m+veps);

		f_num.flux_num[q] = f_num.flux_num[q] = (f_num.v_ch_m*pfP.flux_phys[q] +
				f_num.v_ch_p*pfM.flux_phys[q] -
				f_num.v_ch_m*f_num.v_ch_p*(pfP.uCon[q] - pfM.uCon[q]))*fac;
	}

#if EXTRACT_PRESSURE == TRUE
	if(pfM.isgeneric() == true) {

		REAL fac = 1./(f_num.v_ch_p+f_num.v_ch_m+veps);

#if(FLUID_TYPE==CRONOS_MULTIFLUID)
		cerr << " Not possible in hll " << endl;
		exit(3);
#else
		f_num.ptotal_num = (f_num.v_ch_m*pfP.ptotal +
				f_num.v_ch_p*pfM.ptotal)*fac;
#endif
	}

#endif

}


void HLLSolver::get_NumEmf2D(Data &gdata,
                             phys_fields_2D &pfLL,
                             phys_fields_2D &pfLR,
                             phys_fields_2D &pfRL,
                             phys_fields_2D &pfRR,
                             fields_2D &fl,
                             const int &dir) {

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
		for (int j = 0; j <= gdata.mx[dir1]+1; ++j){
			for (int i = 0; i <= gdata.mx[dir0]+1; ++i){

				REAL fac0 = 1./(fl.v_cor_p[0](i,j)+fl.v_cor_m[0](i,j)+veps);
				REAL fac1 = 1./(fl.v_cor_p[1](i,j)+fl.v_cor_m[1](i,j)+veps);
				REAL fac = fac0*fac1;

				fl.emf(i,j) = (fl.v_cor_m[0](i,j)*fl.v_cor_m[1](i,j)*pfLL.emf(i,j) +
				               fl.v_cor_p[0](i,j)*fl.v_cor_m[1](i,j)*pfRL.emf(i-1,j) +
				               fl.v_cor_m[0](i,j)*fl.v_cor_p[1](i,j)*pfLR.emf(i,j-1) +
				               fl.v_cor_p[0](i,j)*fl.v_cor_p[1](i,j)*pfRR.emf(i-1,j-1))*fac +
					fl.v_cor_m[0](i,j)*fl.v_cor_p[0](i,j)*(pfLL.B_z(i  ,j) -
					                                       pfRL.B_z(i-1,j))*fac0 -
					fl.v_cor_m[1](i,j)*fl.v_cor_p[1](i,j)*(pfLL.B_y(i  ,j) -
					                                       pfLR.B_y(i,j-1))*fac1;
			}
		}
	} else if(dir == 1) {
		for (int j = 0; j <= gdata.mx[dir1]+1; ++j){
			for (int i = 0; i <= gdata.mx[dir0]+1; ++i){

				REAL fac0 = 1./(fl.v_cor_p[0](i,j)+fl.v_cor_m[0](i,j)+veps);
				REAL fac1 = 1./(fl.v_cor_p[1](i,j)+fl.v_cor_m[1](i,j)+veps);
				REAL fac = fac0*fac1;
	      
				fl.emf(i,j) = (fl.v_cor_m[0](i,j)*fl.v_cor_m[1](i,j)*pfLL.emf(i,j) +
				               fl.v_cor_p[0](i,j)*fl.v_cor_m[1](i,j)*pfRL.emf(i-1,j) +
				               fl.v_cor_m[0](i,j)*fl.v_cor_p[1](i,j)*pfLR.emf(i,j-1) +
				               fl.v_cor_p[0](i,j)*fl.v_cor_p[1](i,j)*pfRR.emf(i-1,j-1))*fac +
					fl.v_cor_m[0](i,j)*fl.v_cor_p[0](i,j)*(pfLL.B_x(i  ,j) -
					                                       pfRL.B_x(i-1,j))*fac0 -
					fl.v_cor_m[1](i,j)*fl.v_cor_p[1](i,j)*(pfLL.B_z(i  ,j) -
					                                       pfLR.B_z(i,j-1))*fac1;

			}
		}
	} else {
		for (int j = 0; j <= gdata.mx[dir1]+1; ++j){
			for (int i = 0; i <= gdata.mx[dir0]+1; ++i){

				REAL fac0 = 1./(fl.v_cor_p[0](i,j)+fl.v_cor_m[0](i,j)+veps);
				REAL fac1 = 1./(fl.v_cor_p[1](i,j)+fl.v_cor_m[1](i,j)+veps);
				REAL fac = fac0*fac1;

				fl.emf(i,j) = (fl.v_cor_m[0](i,j)*fl.v_cor_m[1](i,j)*pfLL.emf(i,j) +
				               fl.v_cor_p[0](i,j)*fl.v_cor_m[1](i,j)*pfRL.emf(i-1,j) +
				               fl.v_cor_m[0](i,j)*fl.v_cor_p[1](i,j)*pfLR.emf(i,j-1) +
				               fl.v_cor_p[0](i,j)*fl.v_cor_p[1](i,j)*pfRR.emf(i-1,j-1))*fac +
					fl.v_cor_m[0](i,j)*fl.v_cor_p[0](i,j)*(pfLL.B_y(i  ,j) -
					                                       pfRL.B_y(i-1,j))*fac0 -
					fl.v_cor_m[1](i,j)*fl.v_cor_p[1](i,j)*(pfLL.B_x(i  ,j) -
					                                       pfLR.B_x(i,j-1))*fac1;

			}
		}
	}


}

