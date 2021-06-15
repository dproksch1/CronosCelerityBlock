#include "PhysFluxesHD.H"

PhysFluxesHD::PhysFluxesHD(const Data &gdata,
		const CronosFluid &fluid) : PhysFluxes(gdata, fluid) {

}

void PhysFluxesHD::get_PhysFlux(Queue queue, Data &gdata,
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
//	int q_Eges = gdata.fluids->get_q_Eges(iFluid);
//	int q_Eadd = gdata.fluids->get_q_Eadd(iFluid);
//	cout << " Getting q_Eges " << gdata.fluids->get_q_Eges(iFluid) << " " << iFluid << endl;
//#endif

	pf.fluxORIG[q_sx].clear();
	pf.fluxORIG[q_sy].clear();
	pf.fluxORIG[q_sz].clear();

	clearBuffer(queue, pf.fluxSYCL[q_sx]);
	clearBuffer(queue, pf.fluxSYCL[q_sy]);
	clearBuffer(queue, pf.fluxSYCL[q_sz]);

	//assert_sycl_eq(pf.fluxORIG[q_sx], pf.fluxSYCL[q_sx]);
	//assert_sycl_eq(pf.fluxORIG[q_sy], pf.fluxSYCL[q_sy]);
	//assert_sycl_eq(pf.fluxORIG[q_sz], pf.fluxSYCL[q_sz]);

	// flux for the density
	for (int i = -1; i <= gdata.mx[dir]+1; ++i){

		// pf.flux[q_rho](i) = pf.uCon[q_sx+dir](i);
		pf.fluxORIG[q_rho](i) = pf.uPriORIG[q_sx+dir](i)*pf.uPriORIG[q_rho](i) ;

	}

	int size = gdata.mx[dir];

#ifdef USE_SYCL
	queue.submit([&](sycl::handler& cgh) {
		auto flux_acc = pf.fluxSYCL[q_rho].get_access<cl::sycl::access::mode::write>(cgh);
		auto uPriRho_acc = pf.uPriSYCL[q_rho].get_access<cl::sycl::access::mode::read>(cgh);
		auto uPriOther_acc = pf.uPriSYCL[q_sx + dir].get_access<cl::sycl::access::mode::read>(cgh);
		cgh.parallel_for<class FluxDensity>(Range<1>(size+3), Id<1>(1), [=](Item<1> item) {
			size_t i = item.get(0);
			flux_acc[i] = uPriRho_acc[i] * uPriOther_acc[i];
		});
	});
#endif // USE_SYCL

	//assert_sycl_eq(pf.fluxSYCL[q_rho], pf.fluxORIG[q_rho]);

	// flux for momentum vector
#if EXTRACT_PRESSURE == TRUE
	for(int q=q_sx; q<=q_sz; ++q) {
		for (int i = -1; i <= gdata.mx[dir]+1; ++i){
			pf.flux[q](i) = (pf.uCon[q  ](i)*pf.uPri[dir+q_sx](i));
//			pf.flux[q](i) = (pf.uPri[q  ](i)*pf.uPri[dir+q_sx](i)*pf.uPri[q_rho](i));
		}
	}
#else
	for(int q=q_sx; q<=q_sz; ++q) {
		if(q == dir+1) {
			for (int i = -1; i <= gdata.mx[dir] + 1; ++i) {
				pf.fluxORIG[q](i) += (pf.uPriORIG[q](i)*pf.uConORIG[q](i) + pf.pthermORIG(i));
//				pf.flux[q](i) += (pf.uPri[q](i)*pf.uPri[q](i)*pf.uPri[q_rho](i)
//				                  + pf.ptherm(i));
			}
		} else {
			for (int i = -1; i <= gdata.mx[dir] + 1; ++i) {
				pf.fluxORIG[q](i) = (pf.uConORIG[q](i)*pf.uPriORIG[dir+q_sx](i));
//				pf.flux[q](i) = (pf.uPri[q](i)*pf.uPri[dir+q_sx](i)*pf.uPri[q_rho](i));
			}
		}
	}

#ifdef USE_SYCL
	for (int q = q_sx; q <= q_sz; ++q) {
		if (q == dir + 1) {
			queue.submit([&](sycl::handler& cgh) {
				auto flux_acc = pf.fluxSYCL[q].get_access<cl::sycl::access::mode::read_write>(cgh);
				auto uCon_acc = pf.uConSYCL[q].get_access<cl::sycl::access::mode::read>(cgh);
				auto ptherm_acc = pf.pthermSYCL.get_access<cl::sycl::access::mode::read>(cgh);
				auto uPri_acc = pf.uPriSYCL[q].get_access<cl::sycl::access::mode::read>(cgh);
				cgh.parallel_for<class FluxMomentumA>(Range<1>(size + 3), Id<1>(1), [=](Item<1> item) {
					size_t i = item.get(0);
					flux_acc[i] += (uPri_acc[i] * uCon_acc[i] + ptherm_acc[i]);
				});
			});
		} else {
			queue.submit([&](sycl::handler& cgh) {
				auto flux_acc = pf.fluxSYCL[q].get_access<cl::sycl::access::mode::read_write>(cgh);
				auto uCon_acc = pf.uConSYCL[q].get_access<cl::sycl::access::mode::read>(cgh);
				auto ptherm_acc = pf.pthermSYCL.get_access<cl::sycl::access::mode::read>(cgh);
				auto uPri_acc = pf.uPriSYCL[dir + q_sx].get_access<cl::sycl::access::mode::read>(cgh);
				cgh.parallel_for<class FluxMomentumB>(Range<1>(size + 3), Id<1>(1), [=](Item<1> item) {
					size_t i = item.get(0);
					flux_acc[i] += (uPri_acc[i] * uCon_acc[i]);
				});
			});
		}
	}
#endif // USE_SYCL

#endif


	//for (int q = 0; q <= 4; ++q) {
	//	assert_sycl_eq(pf.fluxORIG[q], pf.fluxSYCL[q]);
	//	assert_sycl_eq(pf.uConORIG[q], pf.uConSYCL[q]);
	//	assert_sycl_eq(pf.uPriORIG[q], pf.uPriSYCL[q]);
	//	assert_sycl_eq(pf.pthermORIG, pf.pthermSYCL);
	//}
   
	// flux for total energy
	if(ENERGETICS == FULL){
		for (int i = -1; i <= gdata.mx[dir]+1; ++i){

#if (USE_COROTATION == CRONOS_ON)
			pf.flux[q_Eges](i) = pf.uCon[q_Eges](i)*pf.uPri[q_sx+dir](i) +
					pf.ptherm(i)*pf.uInertial[dir](i);
#else
			pf.fluxORIG[q_Eges](i) = (pf.uConORIG[q_Eges](i) +
			                      pf.pthermORIG(i))*pf.uPriORIG[q_sx+dir](i);
#endif

		}


#ifdef USE_SYCL
#if (USE_COROTATION == CRONOS_ON)
#else
		queue.submit([&](sycl::handler& cgh) {
			auto flux_acc = pf.fluxSYCL[q_Eges].get_access<cl::sycl::access::mode::read_write>(cgh);
			auto uCon_acc = pf.uConSYCL[q_Eges].get_access<cl::sycl::access::mode::read>(cgh);
			auto ptherm_acc = pf.pthermSYCL.get_access<cl::sycl::access::mode::read>(cgh);
			auto uPri_acc = pf.uPriSYCL[q_sx + dir].get_access<cl::sycl::access::mode::read>(cgh);
			cgh.parallel_for<class FluxEnergy>(Range<1>(size + 3), Id<1>(1), [=](Item<1> item) {
				size_t i = item.get(0);
				flux_acc[i] = (uCon_acc[i] + ptherm_acc[i]) * uPri_acc[i];
			});
		});
#endif

#endif // USE_SYCL
	}

	//assert_sycl_eq(pf.fluxORIG[q_Eges], pf.fluxSYCL[q_Eges]);


#if(CRSWITCH_DUAL_ENERGY == CRONOS_ON)
	// Compute physical flux for entropy evolution -- flux is
	// Entropy*v_par
	for (int i = -1; i <= gdata.mx[dir]+1; ++i){
		pf.flux[q_Eadd](i) = pf.uCon[q_Eadd](i)*pf.uPri[q_sx+dir](i);
	}
#endif

#if (USE_ANGULAR_MOMENTUM == TRUE)

	int ishift = static_cast<int>((2. + 1.e-4)*shift );
	int shift_vec[3] = {0,0,0};

	for (int i = -1; i <= gdata.mx[dir]+1; ++i){
		iPos.set(dir, i);

		pf.flux[2](i) *= gdata.h1(iPos.get(0), iPos.get(1), iPos.get(2),
		                          shift_vec[0],shift_vec[1],shift_vec[2]);
		
	}
		
#endif

  
}


void PhysFluxesHD::get_PhysFlux(const Data &gdata,
	const ProblemType &Problem, phys_fields_0D &fields,
        int /*ix*/, int /*iy*/, int /*iz*/, int face, int iFluid) {

	int dir = face/2;

	fields.flux_phys[q_sx] = 0.;
	fields.flux_phys[q_sy] = 0.;
	fields.flux_phys[q_sz] = 0.;

	// flux for the density
	fields.flux_phys[q_rho] = fields.uPri[q_sx+dir]*fields.uPri[q_rho] ;


	// flux for momentum vector
#if EXTRACT_PRESSURE == TRUE
	for(int q=q_sx; q<=q_sz; ++q) {
		fields.flux_phys[q] = (fields.uCon[q  ]*fields.uPri[dir+q_sx]);
	}
#else
	for(int q=q_sx; q<=q_sz; ++q) {
		if(q == dir+1) {
			fields.flux_phys[q] += (fields.uPri[q]*fields.uCon[q] + fields.ptherm);
		} else {
			fields.flux_phys[q] = (fields.uCon[q]*fields.uPri[dir+q_sx]);
		}
	}
#endif


	// flux for total energy
	if(ENERGETICS == FULL){

#if (USE_COROTATION == CRONOS_ON)
		fields.flux_phys[q_Eges] = fields.uCon[q_Eges]*fields.uPri[q_sx+dir] +
				fields.ptherm*fields.uInertial[dir];
#else
		fields.flux_phys[q_Eges] = (fields.uCon[q_Eges] +
				fields.ptherm)*fields.uPri[q_sx+dir];
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

	fields.flux_phys[2] *= gdata.h1(ix, iy, iz,
			shift_vec[0],shift_vec[1],shift_vec[2]);


#endif

}




