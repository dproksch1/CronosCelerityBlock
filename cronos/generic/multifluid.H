#ifndef CRONOS_MULTIFLUID_H_INC
#define CRONOS_MULTIFLUID_H_INC = 1

#include "matrix.H"
#include "cronos.H"
#include "fluid.H"

class CronosMultifluid {
public:
		CronosMultifluid(int numFluids, int energetics, bool apply_userFields);
		void setup();
		void set_dualEnergy(int aux_Energy);
		void unset_dualEnergy();
		void compute_Variables(int n_add, int n_subs);
		int get_numFluids() const;
		int get_fluidType(int iFluid) const;
		int get_N_OMINT(int iFluid) const;
		int get_N_OMINT() const;
		int get_N_OMEGA() const;
		int get_N_OMINT_ALL() const;
		int get_N_OM() const;
		int get_N_OMINT_USER() const;
		int get_N_OMEGA_USER() const;
		int get_N_OM_USER() const;
		int get_i_magFluid() const;
		int get_q_rho(int iFluid) const;
		int get_q_sx(int iFluid) const;
		int get_q_sy(int iFluid) const;
		int get_q_sz(int iFluid) const;
		int get_q_Bx(int iFluid) const;
		int get_q_By(int iFluid) const;
		int get_q_Bz(int iFluid) const;
		int get_q_Bx() const;
		int get_q_By() const;
		int get_q_Bz() const;
		int get_q_Eges(int iFluid) const;
		int get_q_Eadd(int iFluid) const;
		int get_IndexLocal(int iFieldGlobal) const;
		int get_FluidIndex(int iFieldGlobal) const;
		bool with_magField();
		CronosFluid *fluids;
private:
		int numFluids, N_OMINT_all;
		int n_omInt_total, n_Omega_total, n_om_total;
		int n_omIntUser_total, n_omUser_total;
		NumArray<int> N_OMINTs, N_OMINTsUser;
		NumArray<int> fluidTypes;
		NumArray<int> userFields;
		NumArray<int> index_local, index_local_user;
		NumArray<int> fluidIndex;
//		NumArray<int> q_rho, q_sx, q_sy, q_sz;
//		NumArray<int> q_rho, q_sx, q_sy, q_sz, q_Bx, q_By, q_Bz;
//		NumArray<int> q_Eges, q_Eadd;
//		int q_Bx, q_By, q_Bz;
		int i_magFluid;
		int energetics;
		int aux_energyType;
		bool apply_userFields;
		bool use_dualEnergy;
		bool has_magField;
};

#endif
