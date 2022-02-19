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
	void TransPrim2Cons(Data &gdata, gridFunc &gfunc, ProblemType &Problem);
	void TransCons2Prim(Data &gdata, gridFunc &gfunc, ProblemType &Problem);
    void TransMomen2Vel(Data &, gridFunc &, ProblemType &);
    void TransVel2Momen(Data &, gridFunc &, ProblemType &);
    void TransEth2E(const Data &, gridFunc &, ProblemType &) const;
    void TransE2Eth(Data &, gridFunc &, ProblemType &);
	void TransE2Eth(Data &, gridFunc &, ProblemType &, int, bool);
	double TransE2T(Data &, gridFunc &, ProblemType &);
    void TransT2E(const Data &, gridFunc &, ProblemType &) const;
	void TransT2Eth(const Data &, gridFunc &, ProblemType &) const;
    double TransEth2T(Data &, gridFunc &, ProblemType &);
	void get_Cons(const Data &, const ProblemType &, const EquationOfState &,
			phys_fields_0D &, int ix, int iy, int iz, int face);
    void set_thermal(bool thermal);

	double TransEth2E(double, double, double, double) const;
	double TransT2E(const ProblemType &, double, double, double, double) const;

#if (USE_COROTATION == CRONOS_ON)
	void TransCorotToInert(Data &, gridFunc &, ProblemType &);
	void TransInertToCorot(Data &, gridFunc &, ProblemType &);
	double TransCorotToInert_x(Data &, double vCorot, int ix, int iy, int iz);
	double TransCorotToInert_y(Data &, double vCorot, int ix, int iy, int iz);
	double TransCorotToInert_z(Data &, double vCorot, int ix, int iy, int iz);
	double TransInertToCorot_x(Data &, double vInert, int ix, int iy, int iz);
	double TransInertToCorot_y(Data &, double vInert, int ix, int iy, int iz);
	double TransInertToCorot_z(Data &, double vInert, int ix, int iy, int iz);
//	void setOmega(double OmegaZ);
	void src_Corotating(Data &, ProblemType &, NumMatrix<double, 3> nom[]);
	void store_uInert(Data &gdata, phys_fields_0D &fields, int ix, int iy, int iz);
#endif

private:
	double TransEth2E_HD(double, double, double, double);
	double TransEth2E_MHD(double, double, double, double) const;
	double TransT2E_HD(ProblemType &, double, double, double, double);
	double TransT2E_MHD(const ProblemType &, double, double, double, double) const;
    double TNorm;
	double kBoverMeanMolWeight_num;
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