#include "CTLondrilloDelZanna.H"

CTLondrilloDelZanna::CTLondrilloDelZanna(const Data &gdata, int dir,
                                         bool IntegrateA, int i_magFluid) : ConstrainedTransport(gdata, dir, IntegrateA, i_magFluid)
{
	this->veps = 1.e-120;
}


void CTLondrilloDelZanna::Reconstruct(const Data &gdata) {

	// Doing the reconstruction:
	ReconstEmf->compute(gdata, *fieldsEmf, *physValEmfLL, *physValEmfLR,
	                    *physValEmfRL, *physValEmfRR);

}


void CTLondrilloDelZanna::get_PhysEmfs(const Data &gdata,
                                       const ProblemType &Problem,
                                       phys_fields_2D &pf,
                                       cronos::vector<REAL> &ipos) 
{
	if(dir == 0) {
		
		pf.emf = pf.v_z*pf.B_y - pf.v_y*pf.B_z;
		
	} else if (dir == 1) {
		
		pf.emf = pf.v_x*pf.B_z - pf.v_z*pf.B_x; 

	} else if (dir == 2) {
		
		pf.emf = pf.v_y*pf.B_x - pf.v_x*pf.B_y;
		
	}

}



void CTLondrilloDelZanna::get_vChar2D(const Data &gdata,
                                      const ProblemType &Problem,
                                      cronos::vector<REAL> &ipos)
{
	for (int i = -1; i <= gdata.mx[dir0]+1; ++i){
		for (int j = -1; j <= gdata.mx[dir1]+1; ++j){
			fieldsEmf->v_cor_m[0](i,j) = std::max(fieldsEmf->v_ch_m[0](i,j  ),
			                                      fieldsEmf->v_ch_m[0](i,j-1));
			fieldsEmf->v_cor_p[0](i,j) = std::max(fieldsEmf->v_ch_p[0](i,j  ),
			                                      fieldsEmf->v_ch_p[0](i,j-1));
			fieldsEmf->v_cor_m[1](i,j) = std::max(fieldsEmf->v_ch_m[1](i  ,j),
			                                      fieldsEmf->v_ch_m[1](i-1,j));
			fieldsEmf->v_cor_p[1](i,j) = std::max(fieldsEmf->v_ch_p[1](i  ,j),
			                                      fieldsEmf->v_ch_p[1](i-1,j));
		}
	}
}


void CTLondrilloDelZanna::get_NumEmf2D(const Data &gdata,
                                       const phys_fields_2D &pfLL,
                                       const phys_fields_2D &pfLR,
                                       const phys_fields_2D &pfRL,
                                       const phys_fields_2D &pfRR,
                                       fields_2D &fl)
{
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

