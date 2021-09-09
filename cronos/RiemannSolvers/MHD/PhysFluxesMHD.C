#include "PhysFluxesMHD.H"

PhysFluxesMHD::PhysFluxesMHD(const Data &gdata,
		const CronosFluid &fluid) : PhysFluxes(gdata, fluid) {
	DQ_SB = 3;
	q_rho = fluid.get_q_rho();
	q_sx = fluid.get_q_sx();
	q_sy = fluid.get_q_sy();
	q_sz = fluid.get_q_sz();
	q_Bx = fluid.get_q_Bx();
	q_By = fluid.get_q_By();
	q_Bz = fluid.get_q_Bz();
	q_Eges = fluid.get_q_Eges();
	q_Eadd = fluid.get_q_Eadd();
}

void PhysFluxesMHD::get_PhysFlux(Queue queue, Data &gdata,
                                 ProblemType &Problem,
                                 cronos::vector<REAL> &iPos,
                                 phys_fields_1D &pf,
                                 const int &dir, REAL shift, int iFluid){
  
//#if(FLUID_TYPE == CRONOS_MULTIFLUID)
////	int q_rho = Problem.q_rho;
////	int q_sx = Problem.q_sx;
////	int q_sy = Problem.q_sy;
////	int q_sz = Problem.q_sz;
////	int q_Eges = Problem.q_Eges;
////	int q_Eadd = Problem.q_Eadd;
//	int q_rho = gdata.fluids->get_q_rho(iFluid);
//	int q_sx = gdata.fluids->get_q_sx(iFluid);
//	int q_sy = gdata.fluids->get_q_sy(iFluid);
//	int q_sz = gdata.fluids->get_q_sz(iFluid);
//	int q_Bx = gdata.fluids->get_q_Bx(iFluid);
//	int q_By = gdata.fluids->get_q_By(iFluid);
//	int q_Bz = gdata.fluids->get_q_Bz(iFluid);
//	int q_Eges = gdata.fluids->get_q_Eges(iFluid);
//	int q_Eadd = gdata.fluids->get_q_Eadd(iFluid);
//#endif
	// shift between magnetic field and velocity components
	DQ_SB = q_Bx - q_sx;

	pf.fluxORIG[q_sx].clear();
	pf.fluxORIG[q_sy].clear();
	pf.fluxORIG[q_sz].clear();

	// flux for mass density
	for (int i = -1; i <= gdata.mx[dir]+1; ++i){

		// pf.flux[q_rho](i) = pf.uCon[q_sx+dir](i);
		pf.fluxORIG[q_rho](i) = pf.uPriORIG[q_sx+dir](i)*pf.uPriORIG[q_rho](i);

	}

	// flux for momentum vector
#if EXTRACT_PRESSURE == TRUE
	for(int q=q_sx; q<=q_sz; ++q) {
		for (int i = -1; i <= gdata.mx[dir]+1; ++i){
			pf.flux[q](i) = (pf.uCon[q  ](i)*pf.uPri[dir+q_sx](i) -
					pf.uPri[q+DQ_SB](i)*pf.uPri[dir+DQ_SB+q_sx](i));
//			pf.flux[q](i) = (pf.uPri[q  ](i)*pf.uPri[dir+q_sx](i)*pf.uPri[q_rho](i) -
//			                 pf.uPri[q+DQ_SB](i)*pf.uPri[dir+DQ_SB+q_sx](i));
		}
	}
#else
	for(int q=q_sx; q<=q_sz; ++q) {
		if(q == dir+1) {
			for (int i = -1; i <= gdata.mx[dir]+1; ++i){
				pf.fluxORIG[q](i) += (pf.uPriORIG[q](i)*pf.uConORIG[q](i) +
						pf.pthermORIG(i) -
						0.5*sqr(pf.uPriORIG[dir+DQ_SB+q_sx](i)));
//				pf.flux[q](i) += (pf.uPri[q](i)*pf.uPri[q](i)*pf.uPri[q_rho](i)+
//				                  pf.ptherm(i) -
//				                  0.5*sqr(pf.uPri[dir+DQ_SB+q_sx](i)));
			}
		} else {
			for (int i = -1; i <= gdata.mx[dir]+1; ++i){
				pf.fluxORIG[q](i) = (pf.uConORIG[q](i)*pf.uPriORIG[dir+q_sx](i) -
						pf.uPriORIG[q+DQ_SB](i)*pf.uPriORIG[dir+DQ_SB+q_sx](i));
//				pf.flux[q](i) = (pf.uPri[q](i)*pf.uPri[dir+q_sx](i)*pf.uPri[q_rho](i) -
//				                 pf.uPri[q+DQ_SB](i)*pf.uPri[dir+DQ_SB+q_sx](i));
				pf.fluxORIG[dir+q_sx](i) += 0.5*sqr(pf.uPriORIG[q+DQ_SB](i));
			}
		}
	}
#endif

#ifdef FLUX_CT
	if(Problem.mag) {

		for(int dirVar=0; dirVar<3; dirVar++) {
			for (int i = -1; i <= gdata.mx[dir]+1; ++i){
//#if (USE_COROTATION == CRONOS_ON)
//				pf.flux[dirVar+q_Bx](i) = (pf.uInertial[dir](i)*pf.uPri[dirVar+q_Bx](i) -
//						pf.uInertial[dirVar](i)*pf.uPri[dir   +q_Bx](i));
//#else
				pf.flux[dirVar+q_Bx](i) = (pf.uPri[dir   +q_sx](i)*pf.uPri[dirVar+q_Bx](i) -
				                           pf.uPri[dirVar+q_sx](i)*pf.uPri[dir   +q_Bx](i));
//#endif
			}
		}

	}
#endif

    
	if(ENERGETICS == FULL){
		for (int i = -1; i <= gdata.mx[dir]+1; ++i){
      
#if (USE_COROTATION == CRONOS_ON)
			pf.flux[q_Eges](i) = (pf.uCon[q_Eges](i)*pf.uPri[q_sx+dir](i) + pf.ptotal(i)*pf.uInertial[dir](i) -
					(pf.uInertial[0](i)*pf.uPri[q_Bx](i) +
							pf.uInertial[1](i)*pf.uPri[q_By](i) +
							pf.uInertial[2](i)*pf.uPri[q_Bz](i))*pf.uPri[q_Bx+dir](i));
#else
			pf.fluxORIG[q_Eges](i) = ((pf.uConORIG[q_Eges](i) + pf.ptotalORIG(i))*pf.uPriORIG[q_sx+dir](i) -
			                      (pf.uPriORIG[q_sx](i)*pf.uPriORIG[q_Bx](i) +
			                       pf.uPriORIG[q_sy](i)*pf.uPriORIG[q_By](i) +
			                       pf.uPriORIG[q_sz](i)*pf.uPriORIG[q_Bz](i))*pf.uPriORIG[q_Bx+dir](i));
#endif

		}
	}

#if(CRSWITCH_DUAL_ENERGY == CRONOS_ON)
	// Compute physical flux for entropy evolution -- flux is
	// Entropy*v_par
	for (int i = -1; i <= gdata.mx[dir]+1; ++i){
		pf.flux[q_Eadd](i) = pf.uCon[q_Eadd](i)*pf.uPri[q_sx+dir](i);
	}
#endif

#if (USE_ANGULAR_MOMENTUM == TRUE)

	// int ishift = static_cast<int>(2.*shift + 1.e-4);
	int ishift = static_cast<int>((2. + 1.e-4)*shift );
	int shift_vec[3] = {0,0,0};
	// cout << " shift: " << shift << " " << ishift << " " << dir << endl;
	// cout << "        ";
	// shift_vec[dir] = ishift;
	// cout << shift_vec[0] << " ";
	// cout << shift_vec[1] << " ";
	// cout << shift_vec[2] << " ";
	// cout << endl;

	for (int i = -1; i <= gdata.mx[dir]+1; ++i){
		iPos.set(dir, i);

		pf.flux[q_sy](i) *= gdata.h1(iPos.get(0), iPos.get(1), iPos.get(2),
		                          shift_vec[0],shift_vec[1],shift_vec[2]);
		
	}

	// int ishift = static_cast<int>(2.*shift + 1.e-4);
	// if(dir == 0) {
	// 	pf.flux[2] *= gdata.h1(i,iPos.get(1),iPos.get(2),ishift,0,0);
	// } else if (dir == 1) {
	// 	pf.flux[2] *= gdata.h1(i,j,k,0,ishift,0);
	// } else {
	// 	pf.flux[2] *= gdata.h1(i,j,k,0,0,ishift);
	// }
		
#endif

  
}



void PhysFluxesMHD::get_PhysFlux(/*Queue queue, */const Data &gdata, const ProblemType &Problem, phys_fields_0D &fields,
		int ix, int iy, int iz, int face, int iFluid) {

	int dir = face/2;

	fields.flux_phys[q_sx] = 0.;
	fields.flux_phys[q_sy] = 0.;
	fields.flux_phys[q_sz] = 0.;


	// flux for mass density
	fields.flux_phys[q_rho] = fields.uPri[q_sx+dir]*fields.uPri[q_rho];


	// flux for momentum vector
#if EXTRACT_PRESSURE == TRUE
	for(int q=q_sx; q<=q_sz; ++q) {
		fields.flux_phys[q] = (fields.uCon[q  ]*fields.uPri[dir+q_sx] -
				fields.uPri[q+DQ_SB]*fields.uPri[dir+DQ_SB+q_sx]);
	}
#else
	for(int q=q_sx; q<=q_sz; ++q) {
		if(q == dir+1) {
			fields.flux_phys[q] += (fields.uPri[q]*fields.uCon[q] +
					fields.ptherm - 0.5*sqr(fields.uPri[dir+DQ_SB+q_sx]));
		} else {
			fields.flux_phys[q] = (fields.uCon[q]*fields.uPri[dir+q_sx] -
					fields.uPri[q+DQ_SB]*fields.uPri[dir+DQ_SB+q_sx]);
			fields.flux_phys[dir+q_sx] += 0.5*sqr(fields.uPri[q+DQ_SB]);
		}
	}
#endif

#ifdef FLUX_CT
	if(Problem.mag) {

		for(int dirVar=0; dirVar<3; dirVar++) {
			//#if (USE_COROTATION == CRONOS_ON)
			//				pf.flux[dirVar+q_Bx](i) = (pf.uInertial[dir](i)*pf.uPri[dirVar+q_Bx](i) -
			//						pf.uInertial[dirVar](i)*pf.uPri[dir   +q_Bx](i));
			//#else
			fields.flux_phys[dirVar+q_Bx] = (fields.uPri[dir   +q_sx]*fields.uPri[dirVar+q_Bx] -
					fields.uPri[dirVar+q_sx]*fields.uPri[dir   +q_Bx]);
			//#endif
		}

	}
#endif


	if(ENERGETICS == FULL){

#if (USE_COROTATION == CRONOS_ON)
		fields.flux_phys[q_Eges] = (fields.uCon[q_Eges]*fields.uPri[q_sx+dir] + fields.ptotal*fields.uInertial[dir] -
				(fields.uInertial[0]*fields.uPri[q_Bx] +
						fields.uInertial[1]*fields.uPri[q_By] +
						fields.uInertial[2]*fields.uPri[q_Bz])*fields.uPri[q_Bx+dir]);
#else
		fields.flux_phys[q_Eges] = ((fields.uCon[q_Eges] + fields.ptotal)*fields.uPri[q_sx+dir] -
				(fields.uPri[q_sx]*fields.uPri[q_Bx] +
						fields.uPri[q_sy]*fields.uPri[q_By] +
						fields.uPri[q_sz]*fields.uPri[q_Bz])*fields.uPri[q_Bx+dir]);
#endif

	}

#if(CRSWITCH_DUAL_ENERGY == CRONOS_ON)
	// Compute physical flux for entropy evolution -- flux is
	// Entropy*v_par
	fields.flux_phys[q_Eadd] = fields.uCon[q_Eadd]*fields.uPri[q_sx+dir];
#endif

#if (USE_ANGULAR_MOMENTUM == TRUE)

	int ishift = 0;
	if(face % 2 > 0) ishift = 1;
	int shift_vec[3] = {0,0,0};
	shift_vec[dir] = ishift;

	fields.flux_phys[q_sy] *= gdata.h1(ix, iy, iz,
			shift_vec[0],shift_vec[1],shift_vec[2]);

#endif


}




void PhysFluxesMHD::get_PhysEmfs(Data &gdata,
                                 ProblemType &Problem,
                                 cronos::vector<REAL> &iPos,
                                 phys_fields_2D &pf,
                                 const int &dir){
	if(dir == 0) {
		pf.emf = pf.v_z*pf.B_y - pf.v_y*pf.B_z;
	} else if (dir == 1) {
		pf.emf = pf.v_x*pf.B_z - pf.v_z*pf.B_x; 
	} else if (dir == 2) {
		pf.emf = pf.v_y*pf.B_x - pf.v_x*pf.B_y;
	}
}

