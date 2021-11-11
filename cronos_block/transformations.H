#ifndef CRONOS_TRANSFORMATIONS_H
#define CRONOS_TRANSFORMATIONS_H

#include "CException.H"
#include "data.H"
#include "gridfunc.H"
#include "problem.H"
#include "eos.H"
#include "vector.H"
#include "fields_local.H"

/*********************************************************\
 *
 * CLASS: Transformations:
 *
 * Supplies routines to convert to and fro between 
 * conservative and different kinds of characteristic 
 * variables
 *
\*********************************************************/

class Transformations {
public:
//	Transformations(ProblemType &, bool);
	Transformations(const CronosFluid&, ProblemType &, bool, int iFluid=0);
	void TransPrim2Cons(Data &gdata, gridFunc &gfunc,
	                    ProblemType &Problem); 				//used
	void TransCons2Prim(Data &gdata, gridFunc &gfunc,
	                    ProblemType &Problem);				//used
	void TransE2Eth(Data &, gridFunc &, ProblemType &, int, bool);  //used
	void TransE2Eth(Data &, gridFunc &, ProblemType &); 			//used
	void TransEth2E(const Data &, gridFunc &, ProblemType &) const; //used
	void TransT2Eth(const Data &, gridFunc &, ProblemType &) const; //used
	REAL TransEth2T(Data &, gridFunc &, ProblemType &);				//used
	void TransT2E(const Data &, gridFunc &, ProblemType &) const;	//used
	REAL TransE2T(Data &, gridFunc &, ProblemType &);				//used
#if (USE_ANGULAR_MOMENTUM == TRUE)
	void TransAngMom2Vel(Data &, gridFunc &, ProblemType &); 		//used
	void TransAngMom2Momen(Data &, gridFunc &, ProblemType &);		//used
	void TransVel2AngMom(Data &, gridFunc &, ProblemType &);		//used
	void TransMomen2AngMom(Data &, gridFunc &, ProblemType &);		//
#endif
#if (USE_COROTATION == CRONOS_ON)
	void TransCorotToInert(Data &, gridFunc &, ProblemType &); //used
	void TransInertToCorot(Data &, gridFunc &, ProblemType &); //used
//	REAL TransCorotToInert(Data &, REAL vCorot, int dir, int ix, int iy, int iz);
//	REAL TransInertToCorot(Data &, REAL vInert, int dir, int ix, int iy, int iz);
	REAL TransCorotToInert_x(Data &, REAL vCorot, int ix, int iy, int iz); //NOT USED
	REAL TransCorotToInert_y(Data &, REAL vCorot, int ix, int iy, int iz); //NOT USED
	REAL TransCorotToInert_z(Data &, REAL vCorot, int ix, int iy, int iz); //NOT USED
	REAL TransInertToCorot_x(Data &, REAL vInert, int ix, int iy, int iz); //used
	REAL TransInertToCorot_y(Data &, REAL vInert, int ix, int iy, int iz); //used
	REAL TransInertToCorot_z(Data &, REAL vInert, int ix, int iy, int iz); //used
//	void setOmega(double OmegaZ);
	void src_Corotating(Data &, ProblemType &, NumMatrix<REAL, 3> nom[]); //used
	void store_uInert(Data &gdata, phys_fields_0D &fields, int ix, int iy, int iz); //used
#endif
	void TransMomen2Vel(Data &, gridFunc &, ProblemType &); //used
	void TransVel2Momen(Data &, gridFunc &, ProblemType &); //used
	void EntropyCorrection(Data &, gridFunc &, ProblemType &, int);//used
	void computeEntropyFromE(Data &, gridFunc &, ProblemType &); //used

	inline REAL TransEth2E_HD(REAL, REAL, REAL, REAL);
	inline REAL TransT2E_HD(ProblemType &, REAL, REAL, REAL, REAL);
	void get_Cons_HD(const Data &, const ProblemType &, const EquationOfState &,
				phys_fields_0D &, int ix, int iy, int iz, int face) const;
	//void get_ConsUser(Data &, phys_fields_1D &, int);
	void set_thermal(bool thermal);
//	void reset_Indices(ProblemType &);
	void reset_Indices(CronosFluid &);
		
private:
	inline double get_EMag(const Data &, int ix, int iy, int iz) const;
	REAL max(const REAL &, const REAL &, const REAL &, const REAL &,
	         const REAL &, const REAL &, const REAL &);
	REAL min(const REAL &, const REAL &, const REAL &, const REAL &,
	         const REAL &, const REAL &, const REAL &);
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
