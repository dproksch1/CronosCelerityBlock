#include "transformations.H"

inline REAL Transformations::TransEth2E_MHD(REAL rhoinv, REAL psq, REAL Bsq, REAL ETherm) const
{
	REAL Energy(ETherm + 0.5*psq*rhoinv + 0.5*Bsq);
	return Energy;
}



inline REAL Transformations::TransT2E_MHD(const ProblemType &Problem,
		REAL rhoinv, REAL psq, REAL Bsq, REAL Temp) const
{
	REAL Energy(1./(rhoinv*(Problem.gamma-1.))*Temp);
	return TransEth2E(rhoinv, psq, Bsq, Energy);
}



void Transformations::get_Cons_MHD(Data &gdata,
                               ProblemType &Problem,
                               EquationOfState  &eos,
                               //cronos::vector<int> &iPos,
                               cronos::vector<REAL> &Pos,
                               phys_fields_1D &fields,
                               int dir, REAL shift) {

#if (USE_ANGULAR_MOMENTUM == TRUE)
	cronos::vector<int> iPos(static_cast<int>(Pos.get(0)*(1. + 1.e-5)),
	                         static_cast<int>(Pos.get(1)*(1. + 1.e-5)),
	                         static_cast<int>(Pos.get(2)*(1. + 1.e-5)));
#endif
	// Compute thermal pressure:
	if(ENERGETICS == FULL) {
		if(thermal) { // Compute thermal pressure from Etherm
			for (int i = -2; i <= gdata.mx[dir]+1; ++i){
				fields.pthermORIG(i) = (Problem.gamma-1)*fields.uPriORIG[q_Eges_loc](i);
			}
		} else { // Compute thermal pressure from Temperature
			for (int i = -2; i <= gdata.mx[dir]+1; ++i){
				fields.pthermORIG(i) = fields.uPriORIG[q_rho_loc](i)*fields.uPriORIG[q_Eges_loc](i);
			}
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
  
	// fields.uCon[q_sx_loc] = fields.uPri[q_sx_loc]*fields.uPri[q_rho_loc]; // v_x -> s_x
	// fields.uCon[q_sy_loc] = fields.uPri[q_sy_loc]*fields.uPri[q_rho_loc]; // v_y -> s_y
	// fields.uCon[q_sz_loc] = fields.uPri[q_sz_loc]*fields.uPri[q_rho_loc]; // v_z -> s_z
	fields.uConORIG[q_sx_loc] = fields.uPriORIG[q_sx_loc]; // save v_x
	fields.uConORIG[q_sy_loc] = fields.uPriORIG[q_sy_loc]; // save v_y
	fields.uConORIG[q_sz_loc] = fields.uPriORIG[q_sz_loc]; // save v_z
	fields.uConORIG[q_sx_loc] *= fields.uPriORIG[q_rho_loc]; // v_x -> s_x
	fields.uConORIG[q_sy_loc] *= fields.uPriORIG[q_rho_loc]; // v_y -> s_y
	fields.uConORIG[q_sz_loc] *= fields.uPriORIG[q_rho_loc]; // v_z -> s_z
  
	fields.uConORIG[q_Bx_loc] = fields.uPriORIG[q_Bx_loc]; // B_x
	fields.uConORIG[q_By_loc] = fields.uPriORIG[q_By_loc]; // B_y
	fields.uConORIG[q_Bz_loc] = fields.uPriORIG[q_Bz_loc]; // B_y

	// When working with the angular momentum: do trafo.
#if (USE_ANGULAR_MOMENTUM == TRUE)
	// I will omit this for the time being - let's see if it is
	// important indeed...
	int ishift = static_cast<int>((2. + 1.e-4)*shift );
	// cout << " shift: " << shift << " " << ishift << " " << dir << endl;
	int shift_vec[3] = {0,0,0};
	shift_vec[dir] = ishift;

	for (int i = -2; i <= gdata.mx[dir]+1; ++i){
		iPos.set(dir,i);
		fields.uCon[q_sy_loc](i) *= gdata.h1(iPos.get(0), iPos.get(1), iPos.get(2),
		                                 shift_vec[0],shift_vec[1],shift_vec[2]);
	}

#endif

	// Thermal energy / temperature / overall energy
	if(ENERGETICS == FULL) {
		for (int i = -2; i <= gdata.mx[dir]+1; ++i){
			REAL rhoinv = 1./fields.uPriORIG[q_rho_loc](i);
      
			REAL Bsq = (sqr(fields.uPriORIG[q_Bx_loc](i)) +
			            sqr(fields.uPriORIG[q_By_loc](i)) +
			            sqr(fields.uPriORIG[q_Bz_loc](i)));
  
			REAL psq = (sqr(fields.uPriORIG[q_sx_loc](i)) + sqr(fields.uPriORIG[q_sy_loc](i)) +
			            sqr(fields.uPriORIG[q_sz_loc](i)))*sqr(fields.uPriORIG[q_rho_loc](i));
      
			if(thermal) {
				fields.uConORIG[q_Eges_loc](i) = TransEth2E(rhoinv,psq,Bsq,fields.uPriORIG[q_Eges_loc](i));
			} else {
				fields.uConORIG[q_Eges_loc](i) = TransT2E(Problem, rhoinv,psq,Bsq,fields.uPriORIG[q_Eges_loc](i));
			}
		}

#if(CRSWITCH_DUAL_ENERGY == CRONOS_ON)
		// Computation of entropy for dual energy corrections
		for (int i = -2; i <= gdata.mx[dir]+1; ++i){
			double Etherm = fields.ptherm(i)/(Problem.gamma - 1.);
			fields.uPri[q_Eadd_loc](i) = Etherm*pow(fields.uCon[q_rho_loc](i),
			                                    1.-Problem.gamma);
			fields.uCon[q_Eadd_loc](i) = fields.uPri[q_Eadd_loc](i);
		}
#endif

	}

	// Total pressure:
	for (int i = -2; i <= gdata.mx[dir]+1; ++i){

		REAL Bsq = (sqr(fields.uPriORIG[q_Bx_loc](i)) +
		            sqr(fields.uPriORIG[q_By_loc](i)) +
		            sqr(fields.uPriORIG[q_Bz_loc](i)));
  
		fields.ptotalORIG(i) = fields.pthermORIG(i) + 0.5*Bsq;
	}

	// store inertial fram velocity and compute co-rotating frame velocity
#if (USE_COROTATION == CRONOS_ON)
	store_uInert(gdata, Pos, fields, dir);
#endif

}




void Transformations::get_Cons_MHD(const Data &gdata, const ProblemType &Problem,
		const EquationOfState  &eos, phys_fields_0D &fields, int ix, int iy, int iz, int face) {

	NumArray<REAL> Pos(3);
	Pos(0) = gdata.getCen_x(ix);
	Pos(1) = gdata.getCen_y(iy);
	Pos(2) = gdata.getCen_z(iz);

	NumArray<REAL> iPos(3);
	iPos(0) = ix;
	iPos(1) = iy;
	iPos(2) = iz;

	int dir = face / 2;

	// Take care of specific face
	switch (face) {
	case 0:
		Pos(0) = gdata.getEdgL_x(ix);
		break;
	case 1:
		Pos(0) = gdata.getEdgL_x(ix+1);
		break;
	case 2:
		Pos(1) = gdata.getEdgL_y(iy);
		break;
	case 3:
		Pos(1) = gdata.getEdgL_y(iy+1);
		break;
	case 4:
		Pos(2) = gdata.getEdgL_z(iz);
		break;
	case 5:
		Pos(2) = gdata.getEdgL_z(iz+1);
		break;
	}

	// Compute thermal pressure:
	if(ENERGETICS == FULL) {
		if(thermal) { // Compute thermal pressure from Etherm
			fields.ptherm = (Problem.gamma-1)*fields.uPri[q_Eges_loc];
		} else { // Compute thermal pressure from Temperature
				fields.ptherm = fields.uPri[q_rho_loc]*fields.uPri[q_Eges_loc];
		}
	} else {
			fields.ptherm = eos.pressure(gdata, Problem, fields.uPri[q_rho_loc],
					Pos(0), Pos(1), Pos(2));
	}

	// Set / compute conservative variables:

	fields.uCon[q_rho_loc] = fields.uPri[q_rho_loc]; // Density

	fields.uCon[q_sx_loc] = fields.uPri[q_sx_loc]; // save v_x
	fields.uCon[q_sy_loc] = fields.uPri[q_sy_loc]; // save v_y
	fields.uCon[q_sz_loc] = fields.uPri[q_sz_loc]; // save v_z
	fields.uCon[q_sx_loc] *= fields.uPri[q_rho_loc]; // v_x -> s_x
	fields.uCon[q_sy_loc] *= fields.uPri[q_rho_loc]; // v_y -> s_y
	fields.uCon[q_sz_loc] *= fields.uPri[q_rho_loc]; // v_z -> s_z

	fields.uCon[q_Bx_loc] = fields.uPri[q_Bx_loc]; // B_x
	fields.uCon[q_By_loc] = fields.uPri[q_By_loc]; // B_y
	fields.uCon[q_Bz_loc] = fields.uPri[q_Bz_loc]; // B_y

	// When working with the angular momentum: do trafo.
#if (USE_ANGULAR_MOMENTUM == TRUE)
	// I will omit this for the time being - let's see if it is
	// important indeed...
	int ishift = 0;
	if(face % 2 > 0) ishift = 1;
	int shift_vec[3] = {0,0,0};
	shift_vec[dir] = ishift;

	fields.uCon[q_sy_loc] *= gdata.h1(ix, iy, iz,
			shift_vec[0],shift_vec[1],shift_vec[2]);

#endif

	// Thermal energy / temperature / overall energy
	if(ENERGETICS == FULL) {
		REAL rhoinv = 1./fields.uPri[q_rho_loc];

		REAL Bsq = (sqr(fields.uPri[q_Bx_loc]) +
				sqr(fields.uPri[q_By_loc]) +
				sqr(fields.uPri[q_Bz_loc]));

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

	REAL Bsq = (sqr(fields.uPri[q_Bx_loc]) +
			sqr(fields.uPri[q_By_loc]) +
			sqr(fields.uPri[q_Bz_loc]));

	fields.ptotal = fields.ptherm + 0.5*Bsq;

	// store inertial fram velocity and compute co-rotating frame velocity
#if (USE_COROTATION == CRONOS_ON)
	store_uInert(gdata, fields, ix, iy, iz);
#endif

}
