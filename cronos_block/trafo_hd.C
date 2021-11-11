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


void Transformations::get_Cons_HD(const Data &gdata, const ProblemType &Problem,
		const EquationOfState  &eos, phys_fields_0D &fields, int face) const {

	int dir = face / 2;

	// Compute thermal pressure:
	if(ENERGETICS == FULL) {
		if(thermal) { // Compute thermal pressure from Etherm
			fields.ptherm = (Problem.gamma-1.)*fields.uPri(q_Eges_loc);
		} else { // Compute thermal pressure from Temperature
			fields.ptherm = fields.uPri(q_rho_loc)*fields.uPri(q_Eges_loc);
		}
	} else {
		NumArray<REAL> Pos(3);
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
			fields.uCon[q_Eges_loc] = TransEth2E_HD(rhoinv,psq,Bsq,
					fields.uPri[q_Eges_loc]);
		} else {
			fields.uCon[q_Eges_loc] = TransT2E_HD(Problem, rhoinv,psq,Bsq,
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
