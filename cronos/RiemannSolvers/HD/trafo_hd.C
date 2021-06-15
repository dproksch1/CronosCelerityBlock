#include "transformations.H"

inline REAL Transformations::TransEth2E_HD(REAL rhoinv, REAL psq, REAL Bsq, REAL ETherm)
{
	REAL Energy(ETherm + 0.5*psq*rhoinv);
	return Energy;
}



inline REAL Transformations::TransT2E_HD(ProblemType &Problem,
		REAL rhoinv, REAL psq, REAL Bsq, REAL Temp)
{
	REAL Energy(1./(rhoinv*(Problem.gamma-1.))*Temp);
	return TransEth2E(rhoinv, psq, Bsq, Energy);
}

namespace transformations_free_functions {

	REAL TransEth2E_MHD(REAL rhoinv, REAL psq, REAL Bsq, REAL ETherm) {
		REAL Energy(ETherm + 0.5 * psq * rhoinv + 0.5 * Bsq);
		return Energy;
	}

	REAL TransT2E_MHD(REAL gamma, REAL rhoinv, REAL psq, REAL Bsq, REAL Temp) {
		REAL Energy(1. / (rhoinv * (gamma - 1.)) * Temp);
		return TransEth2E_MHD(rhoinv, psq, Bsq, Energy);
	}

} // end namespace transformations_free_functions

void Transformations::get_Cons_HD(Queue queue, Data &gdata,
                               ProblemType &Problem,
                               EquationOfState  &eos,
                               cronos::vector<REAL> &Pos,
                               phys_fields_1D &fields,
                               int dir, REAL shift) {

#if (USE_ANGULAR_MOMENTUM == TRUE)
	// I will omit this for the time being - let's see if it is
	// important indeed...
	cronos::vector<int> iPos(static_cast<int>(Pos.get(0)*(1. + 1.e-5)),
	                         static_cast<int>(Pos.get(1)*(1. + 1.e-5)),
	                         static_cast<int>(Pos.get(2)*(1. + 1.e-5)));
#endif

	// Compute thermal pressure:
	if(ENERGETICS == FULL) {
		if(thermal) { // Compute thermal pressure from Etherm
			// ORIG
			for (int i = -2; i <= gdata.mx[dir]+1; ++i){
				fields.pthermORIG(i) = (Problem.gamma-1)*fields.uPriORIG[q_Eges_loc](i);
			}

			// SYCL
			int size = gdata.mx[dir];
			double gamma = Problem.gamma;
			queue.submit([&](cl::sycl::handler& cgh) {
				auto ptherm_acc = fields.pthermSYCL.get_access<cl::sycl::access::mode::write>(cgh);
				auto uPri_acc = fields.uPriSYCL[q_Eges_loc].get_access<cl::sycl::access::mode::read>(cgh);

				cgh.parallel_for<class GetDeriv>(Range<1>(size + 4), Id<1>(0), [=](Item<1> item) {
					size_t i = item.get_range(0);
					ptherm_acc[i] = (gamma - 1) * uPri_acc[i];
				});
			});
			//assert_sycl_eq(fields.pthermSYCL, fields.pthermORIG);
			//assert_sycl_eq(fields.uPriSYCL[q_Eges_loc], fields.uPriORIG[q_Eges_loc]);
		} else { // Compute thermal pressure from Temperature
			// ORIG
			for (int i = -2; i <= gdata.mx[dir]+1; ++i){
				fields.pthermORIG(i) = fields.uPriORIG[q_rho_loc](i)*fields.uPriORIG[q_Eges_loc](i);
			}

			// SYCL
			int size = gdata.mx[dir];
			double gamma = Problem.gamma;
			queue.submit([&](cl::sycl::handler& cgh) {
				auto ptherm_acc = fields.pthermSYCL.get_access<cl::sycl::access::mode::write>(cgh);
				auto uPriEges_acc = fields.uPriSYCL[q_Eges_loc].get_access<cl::sycl::access::mode::read>(cgh);
				auto uPriRho_acc = fields.uPriSYCL[q_rho_loc].get_access<cl::sycl::access::mode::read>(cgh);

				cgh.parallel_for<class GetDeriv>(Range<1>(size + 4), Id<1>(0), [=](Item<1> item) {
					size_t i = item.get_range(0);
					ptherm_acc[i] = uPriRho_acc[i] * uPriEges_acc[i];
				});
			});
			//assert_sycl_eq(fields.pthermSYCL, fields.pthermORIG);
			//assert_sycl_eq(fields.uPriSYCL[q_Eges_loc], fields.uPriORIG[q_Eges_loc]);
		}
	} else {
		for (int i = -2; i <= gdata.mx[dir]+1; ++i){
// 			Pos.set(dir, gdata.getCen_x(dir, i)+shift);

// 			fields.ptherm(i) = eos.pressure(gdata, Problem, fields.uPri[0](i),
// 			                                1., 1., 1.);
			Pos.set(dir, i+shift);
						
			fields.pthermORIG(i) = eos.pressure(gdata, Problem, fields.uPriORIG[q_rho_loc](i),
			                                Pos.get(0), Pos.get(1), Pos.get(2));
		}
	}
  
	// Set / compute conservative variables:
  
	fields.uConORIG[q_rho_loc] = fields.uPriORIG[q_rho_loc]; // Density  
	fields.uConORIG[q_sx_loc] = fields.uPriORIG[q_sx_loc]; // v_x -> s_x
	fields.uConORIG[q_sx_loc] *= fields.uPriORIG[q_rho_loc];
	fields.uConORIG[q_sy_loc] = fields.uPriORIG[q_sy_loc]; // v_y -> s_y
	fields.uConORIG[q_sy_loc] *= fields.uPriORIG[q_rho_loc];
	fields.uConORIG[q_sz_loc] = fields.uPriORIG[q_sz_loc]; // v_z -> s_z
	fields.uConORIG[q_sz_loc] *= fields.uPriORIG[q_rho_loc];
//	fields.uCon[q_sx_loc] = fields.uPriORIG[q_sx_loc]*fields.uPriORIG[q_rho_loc]; // v_x -> s_x
//	fields.uCon[q_sy_loc] = fields.uPriORIG[q_sy_loc]*fields.uPriORIG[q_rho_loc]; // v_y -> s_y
//	fields.uCon[q_sz_loc] = fields.uPriORIG[q_sz_loc]*fields.uPriORIG[q_rho_loc]; // v_z -> s_z

	queue.submit([&](cl::sycl::handler& cgh) {
		auto pri_acc = fields.uPriSYCL[q_rho_loc].get_access<cl::sycl::access::mode::read>(cgh);
		auto con_acc = fields.uConSYCL[q_rho_loc].get_access<cl::sycl::access::mode::write>(cgh);
		cgh.copy(pri_acc, con_acc);
	});

	const auto& accumulate = [](Queue queue, Buffer<REAL, 1> priCur, Buffer<REAL, 1> priRho, Buffer<REAL, 1> con) {
		size_t size = priCur.get_count();
		queue.submit([&](cl::sycl::handler& cgh) {
			auto priCur_acc = priCur.get_access<cl::sycl::access::mode::read_write>(cgh);
			auto priRho_acc = priRho.get_access<cl::sycl::access::mode::read_write>(cgh);
			auto con_acc = con.get_access<cl::sycl::access::mode::write>(cgh);
			cgh.parallel_for<class GetDeriv>(Range<1>(size), [=](Item<1> item) {
				size_t i = item.get_range(0);
				con_acc[i] = priCur_acc[i] * priRho_acc[i];
			});
		});
	};

	accumulate(queue, fields.uPriSYCL[q_sx_loc], fields.uPriSYCL[q_rho_loc], fields.uConSYCL[q_sx_loc]);
	accumulate(queue, fields.uPriSYCL[q_sy_loc], fields.uPriSYCL[q_rho_loc], fields.uConSYCL[q_sy_loc]);
	accumulate(queue, fields.uPriSYCL[q_sz_loc], fields.uPriSYCL[q_rho_loc], fields.uConSYCL[q_sz_loc]);

	//assert_sycl_eq(fields.uConSYCL[q_rho_loc], fields.uConORIG[q_rho_loc]);
	//assert_sycl_eq(fields.uConSYCL[q_sx_loc], fields.uConORIG[q_sx_loc]);
	//assert_sycl_eq(fields.uConSYCL[q_sy_loc], fields.uConORIG[q_sy_loc]);
	//assert_sycl_eq(fields.uConSYCL[q_sz_loc], fields.uConORIG[q_sz_loc]);

	// When working with the angular momentum: do trafo.
#if (USE_ANGULAR_MOMENTUM == TRUE)
	int ishift = static_cast<int>((2. + 1.e-4)*shift );
	// cout << " shift: " << shift << " " << ishift << " " << dir << endl;
	int shift_vec[3] = {0,0,0};
	shift_vec[dir] = ishift;
	// fields.uCon[q_sy_loc] *= gdata.h1(iPos.get(0), iPos.get(1), iPos.get(2),
	//                               shift_vec[0],shift_vec[1],shift_vec[2]);
	for (int i = -2; i <= gdata.mx[dir]+1; ++i){
		// iPos.set(0,i); // That was stupid!
		iPos.set(dir,i);
		fields.uCon[q_sy_loc](i) *= gdata.h1(iPos.get(0), iPos.get(1), iPos.get(2),
		                                 shift_vec[0],shift_vec[1],shift_vec[2]);
	}

	// if(iPos.get(1) == 0 && iPos.get(2)==0 && dir==0) {
	// 	cout << " fac: ";
	// 	cout << gdata.h1(iPos.get(0), iPos.get(1), iPos.get(2),
	// 	                 shift_vec[0],shift_vec[1],shift_vec[2]);
	// 	cout << " ";
	// 	cout << Pos.get(1) << " " << Pos.get(2) << " ";
	// 	cout << iPos.get(1) << " " << iPos.get(2) << " ";
	// 	cout << shift << " " << ishift << " ";
	// 	cout << endl;
	// }
#endif

	REAL Bsq(0.);
	// Thermal energy / temperature / overall energy
	if(ENERGETICS == FULL) {
		for (int i = -2; i <= gdata.mx[dir]+1; ++i){
			REAL rhoinv = 1./fields.uPriORIG[q_rho_loc](i);
  
			REAL psq = (sqr(fields.uPriORIG[q_sx_loc](i)) + sqr(fields.uPriORIG[q_sy_loc](i)) +
			            sqr(fields.uPriORIG[q_sz_loc](i)))*sqr(fields.uPriORIG[q_rho_loc](i));
      
			if(thermal) {
				fields.uConORIG[q_Eges_loc](i) = TransEth2E(rhoinv,psq,Bsq,
				                                    fields.uPriORIG[q_Eges_loc](i));
			} else {
				fields.uConORIG[q_Eges_loc](i) = TransT2E(Problem, rhoinv,psq,Bsq,
				                                  fields.uPriORIG[q_Eges_loc](i));
			}
		}


		int size = gdata.mx[dir];
		double gamma = Problem.gamma;
		queue.submit([&](cl::sycl::handler& cgh) {
			auto uCon_acc = fields.uConSYCL[q_Eges_loc].get_access<cl::sycl::access::mode::write>(cgh);
			auto uPriRho_acc = fields.uPriSYCL[q_rho_loc].get_access<cl::sycl::access::mode::read>(cgh);
			auto uPriEges_acc = fields.uPriSYCL[q_Eges_loc].get_access<cl::sycl::access::mode::read>(cgh);
			auto uPriX_acc = fields.uPriSYCL[q_sx_loc].get_access<cl::sycl::access::mode::read>(cgh);
			auto uPriY_acc = fields.uPriSYCL[q_sy_loc].get_access<cl::sycl::access::mode::read>(cgh);
			auto uPriZ_acc = fields.uPriSYCL[q_sz_loc].get_access<cl::sycl::access::mode::read>(cgh);

			cgh.parallel_for<class ThermalEnergy>(Range<1>(size + 4), Id<1>(0), [=](Item<1> item) {
				size_t i = item.get_range(0);
				REAL rhoinv = 1. / uPriRho_acc[i];

				REAL psq = (sqr(uPriX_acc[i]) + sqr(uPriY_acc[i]) + sqr(uPriZ_acc[i])) * sqr(uPriRho_acc[i]);

				if (thermal) {
					uCon_acc[i] = transformations_free_functions::TransEth2E_MHD(rhoinv, psq, Bsq, uPriEges_acc[i]);
				} else {
					uCon_acc[i] = transformations_free_functions::TransT2E_MHD(gamma, rhoinv, psq, Bsq, uPriEges_acc[i]);
				}
			});
		});

		//assert_sycl_eq(fields.uConORIG[q_Eges_loc], fields.uConSYCL[q_Eges_loc]);

#if(CRSWITCH_DUAL_ENERGY == CRONOS_ON)
		// Computation of entropy for dual energy corrections
		for (int i = -2; i <= gdata.mx[dir]+1; ++i){
			double Etherm = fields.ptherm(i)/(Problem.gamma - 1.);
			fields.uPriORIG[q_Eadd_loc](i) = Etherm*pow(fields.uCon[q_rho_loc](i),
			                                    1.-Problem.gamma);
			fields.uCon[q_Eadd_loc](i) = fields.uPriORIG[q_Eadd_loc](i);
		}
#endif

	}

	// Total pressure:
	fields.ptotalORIG = fields.pthermORIG;

	queue.submit([&](cl::sycl::handler& cgh) {
		auto ptherm_acc = fields.pthermSYCL.get_access<cl::sycl::access::mode::read>(cgh);
		auto ptotal_acc = fields.ptotalSYCL.get_access<cl::sycl::access::mode::write>(cgh);
		cgh.copy(ptherm_acc, ptotal_acc);
	});

	//assert_sycl_eq(fields.ptotalORIG, fields.ptotalSYCL);

	// store inertial fram velocity and compute co-rotating frame velocity
#if (USE_COROTATION == CRONOS_ON)
	store_uInert(gdata, Pos, fields, dir);
#endif

}



void Transformations::get_Cons_HD(const Data &gdata, const ProblemType &Problem,
		const EquationOfState  &eos, phys_fields_0D &fields, int /*ix*/, int /*iy*/, int /*iz*/, int face) const {

	NumArray<REAL> Pos(3);
	//Pos(0) = gdata.getCen_x(ix);
	//Pos(1) = gdata.getCen_y(iy);
	//Pos(2) = gdata.getCen_z(iz);

	//NumArray<REAL> iPos(3);
	//iPos(0) = ix;
	//iPos(1) = iy;
	//iPos(2) = iz;

	int dir = face / 2;

	// Take care of specific face
	//switch (face) {
	//case 0:
	//	Pos(0) = gdata.getEdgL_x(ix);
	//	break;
	//case 1:
	//	Pos(0) = gdata.getEdgL_x(ix+1);
	//	break;
	//case 2:
	//	Pos(1) = gdata.getEdgL_y(iy);
	//	break;
	//case 3:
	//	Pos(1) = gdata.getEdgL_y(iy+1);
	//	break;
	//case 4:
	//	Pos(2) = gdata.getEdgL_z(iz);
	//	break;
	//case 5:
	//	Pos(2) = gdata.getEdgL_z(iz+1);
	//	break;
	//}

	// Compute thermal pressure:
	if(ENERGETICS == FULL) {
		if(thermal) { // Compute thermal pressure from Etherm
			fields.ptherm = (Problem.gamma-1.)*fields.uPri(q_Eges_loc);
		} else { // Compute thermal pressure from Temperature
			fields.ptherm = fields.uPri(q_rho_loc)*fields.uPri(q_Eges_loc);
		}
	} else {
		fields.ptherm = eos.pressure(gdata, Problem, fields.uPri(q_rho_loc),
				Pos(0), Pos(1), Pos(2));
	}

	// Set / compute conservative variables:

	fields.uCon[q_rho_loc] = fields.uPri[q_rho_loc]; // Density

	fields.uCon[q_sx_loc] = fields.uPri[q_sx_loc]*fields.uPri[q_rho_loc]; // v_x -> s_x
	fields.uCon[q_sy_loc] = fields.uPri[q_sy_loc]*fields.uPri[q_rho_loc]; // v_y -> s_y
	fields.uCon[q_sz_loc] = fields.uPri[q_sz_loc]*fields.uPri[q_rho_loc]; // v_z -> s_z

	// When working with the angular momentum: do trafo.
#if (USE_ANGULAR_MOMENTUM == TRUE)

	int ishift = 0;
	if(face % 2 > 0) ishift = 1;
	int shift_vec[3] = {0,0,0};
	shift_vec[dir] = ishift;

	fields.uCon[q_sy_loc] *= gdata.h1(ix, iy, iz,
			shift_vec[0],shift_vec[1],shift_vec[2]);

#endif

	REAL Bsq(0.);
	// Thermal energy / temperature / overall energy
	if(ENERGETICS == FULL) {
		REAL rhoinv = 1./fields.uPri[q_rho_loc];

		REAL psq = (sqr(fields.uPri[q_sx_loc]) + sqr(fields.uPri[q_sy_loc]) +
				sqr(fields.uPri[q_sz_loc]))*sqr(fields.uPri[q_rho_loc]);

		if(thermal) {
			fields.uCon[q_Eges_loc] = TransEth2E(rhoinv,psq,Bsq,
					fields.uPri[q_Eges_loc]);
		} else {
			fields.uCon[q_Eges_loc] = TransT2E(Problem, rhoinv,psq,Bsq,
					fields.uPri[q_Eges_loc]);
		}

#if(CRSWITCH_DUAL_ENERGY == CRONOS_ON)
		// Computation of entropy for dual energy corrections
		double Etherm = fields.ptherm/(Problem.gamma - 1.);
		fields.uPri[q_Eadd_loc] = Etherm*pow(fields.uCon[q_rho_loc],
				1.-Problem.gamma);
		fields.uCon[q_Eadd_loc] = fields.uPri[q_Eadd_loc];
#endif

	}

	// Total pressure:
	fields.ptotal = fields.ptherm;

	// store inertial fram velocity and compute co-rotating frame velocity
#if (USE_COROTATION == CRONOS_ON)
	store_uInert(gdata, fields, ix, iy, iz);
#endif

}
