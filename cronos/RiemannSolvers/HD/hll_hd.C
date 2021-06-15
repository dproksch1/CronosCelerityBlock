#include "RiemannSolverHD.H"

//HLLSolver_gen::HLLSolver_gen(const Data &gdata, int dir, int Fluid_Type) : RiemannSolverHD(gdata, dir, Fluid_Type) {
//	veps = 1.e-120;
//	norm = 1./sqrt(2.);
//}
//
//void HLLSolver_gen::get_NumFlux(Queue /*queue*/, Data &gdata,
//                                phys_fields_1D &pfL,
//                                phys_fields_1D &pfR,
//                                fields_1D &fl, const int &dir, int iFluid) const {
//	for(int q=0; q<fl.get_num(); ++q) {
//		for (int i = -1; i <= gdata.mx[dir]+1; ++i){
//      
//			REAL fac = 1./(fl.v_ch_pORIG(i)+fl.v_ch_mORIG(i)+veps);
//      
//			fl.fluxORIG[q](i) = (fl.v_ch_mORIG(i)*pfL.fluxORIG[q](i) +
//			                 fl.v_ch_pORIG(i)*pfR.fluxORIG[q](i-1) -
//			                 fl.v_ch_mORIG(i)*fl.v_ch_pORIG(i)*(pfL.uConORIG[q](i  )-
//			                                            pfR.uConORIG[q](i-1)))*fac;
//		}
//	}
//
//#if EXTRACT_PRESSURE == TRUE
//	if(fl.isgeneric() == true) {
//		for (int i = -1; i <= gdata.mx[dir]+1; ++i){
//			
//			REAL fac = 1./(fl.v_ch_p(i)+fl.v_ch_m(i)+veps);
//			
//#if(FLUID_TYPE==CRONOS_MULTIFLUID)
//			fl.ptotal[iFluid](i) = (fl.v_ch_m(i)*pfL.ptotal(i) +
//			                fl.v_ch_p(i)*pfR.ptotal(i-1))*fac;
//#else
//			fl.ptotal(i) = (fl.v_ch_m(i)*pfL.ptotal(i) +
//			                fl.v_ch_p(i)*pfR.ptotal(i-1))*fac;
//#endif
//		}
//	}
//
//#endif
//
//}
//
//
//
//void HLLSolver_gen::get_NumFlux(Data &gdata, phys_fields_0D &pfM,
//		phys_fields_0D &pfP, num_fields_0D &f_num, int dir, int iFluid) const {
//	for(int q=0; q<pfM.get_num(); ++q) {
//		REAL fac = 1./(f_num.v_ch_p+f_num.v_ch_m+veps);
//
//		f_num.flux_num[q] = (f_num.v_ch_m*pfP.flux_phys[q] +
//			                 f_num.v_ch_p*pfM.flux_phys[q] -
//			                 f_num.v_ch_m*f_num.v_ch_p*(pfP.uCon[q]- pfM.uCon[q]))*fac;
//	}
//
//#if EXTRACT_PRESSURE == TRUE
//	if(pfL.isgeneric() == true) {
//
//		REAL fac = 1./(f_num.v_ch_p+f_num.v_ch_m+veps);
//
//#if(FLUID_TYPE==CRONOS_MULTIFLUID)
//		cerr << " Not possible in hll_hd " << endl;
//		exit(3);
//		f_num.ptotal[iFluid] = (f_num.v_ch_m*pfP.ptotal +
//				f_num.v_ch_p*pfM.ptotal)*fac;
//#else
//		f_num.ptotal_num = (f_num.v_ch_m*pfP.ptotal +
//				f_num.v_ch_p*pfM.ptotal)*fac;
//#endif
//	}
//
//#endif
//
//}



