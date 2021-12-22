#ifndef CRONOS_TRANSFORMATIONS_H
#define CRONOS_TRANSFORMATIONS_H

#include "CException.H"
#include "data.H"
#include "gridfunc.H"
#include "problem.H"
#include "eos.H"
#include "vector.H"
#include "fields_local.H"

class Transformations {
public:
	Transformations(const CronosFluid&, ProblemType &, bool, int iFluid=0);
    void TransMomen2Vel(Data &, gridFunc &, ProblemType &);
    void TransVel2Momen(Data &, gridFunc &, ProblemType &);
    void TransEth2E(const Data &, gridFunc &, ProblemType &) const;
    void TransE2Eth(Data &, gridFunc &, ProblemType &);
	void TransE2Eth(Data &, gridFunc &, ProblemType &, int, bool);
    void TransT2E(const Data &, gridFunc &, ProblemType &) const;
	void TransT2Eth(const Data &, gridFunc &, ProblemType &) const;
    REAL TransEth2T(Data &, gridFunc &, ProblemType &);
    void set_thermal(bool thermal);
private:
    REAL TNorm;
	REAL kBoverMeanMolWeight_num;
#if (USE_COROTATION == CRONOS_ON)
	double omegaZ;
#endif
	int TempErr, thermal;
	int fluidType, iFluid;
	int q_rho, q_sx, q_sy, q_sz, q_Bx, q_By, q_Bz, q_Eges, q_Eadd;
	int q_rho_loc, q_sx_loc, q_sy_loc, q_sz_loc;
	int q_Bx_loc, q_By_loc, q_Bz_loc, q_Eges_loc, q_Eadd_loc;
	bool TPhys;
	bool magFluid;

};


#endif