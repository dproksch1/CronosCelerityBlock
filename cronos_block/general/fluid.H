#ifndef CRONOS_FLUID_H_INCLUDED
#define CRONOS_FLUID_H_INCLUDED

#include <string>
#include <vector>
#include <matrix.H>

class CronosFluid {
public:
	CronosFluid();
	CronosFluid(int type, int energetics, int n_add, int n_sub,
			int num_user=0, int first_index=0, int first_index_user=0);
	void setup(int type, int energetics, int n_add, int n_sub,
			int num_user=0, int first_index=0, int first_index_user=0);
	void set_dualEnergy(int aux_Energy);
	void unset_dualEnergy();
	int add_Field(std::string type, int first_index);
	void add_Ref(std::string type, int field_index);
	void compute_Variables();
	void compute_Numbers();
	int get_q_rho() const;
	int get_q_sx() const;
	int get_q_sy() const;
	int get_q_sz() const;
	int get_q_Bx() const;
	int get_q_By() const;
	int get_q_Bz() const;
	int get_q_Eges() const;
	int get_q_Eadd() const;
	int get_N_OMINT() const;
	int get_N_OMINT_ALL() const;
	int get_N_OM();
	int get_N_OMEGA();
	int get_N_OMINT_USER();
	int get_N_OM_USER();
	int get_fluid_type() const;
	int get_q_rho_global() const;
	int get_q_sx_global() const;
	int get_q_sy_global() const;
	int get_q_sz_global() const;
	int get_q_Bx_global() const;
	int get_q_By_global() const;
	int get_q_Bz_global() const;
	int get_q_Eges_global() const;
	int get_q_Eadd_global() const;
	int get_IndexGlobal(int iField) const;
	int get_IndexGlobalUser(int iField) const;
	bool has_MagField() const;
	void set_GridTypes();
	void set_IndexGlobal();
	void set_Name(std::string);
	std::string get_Name() const;
	std::string get_fieldName(int iField) const;
private:
	NumArray<int> pos_type;
	NumArray<int> index_global, index_global_user;
//	NumArray<int> is_virtual;
	std::vector<std::string> fieldNames;
	std::vector<bool> is_virtual;
	int numFields, numFieldsUser, numVirtuals, fluid_type;
	int energetics;
	int first_index, first_index_user;
	int num_user;
	int n_omInt, n_Omega, n_om;
	int n_omIntFull, n_omIntAll;
	int n_omIntUser, n_omUser;
	int n_add, n_subs;
	// fluid-local indices
	int q_rho, q_sx, q_sy, q_sz, q_Eges, q_Eadd;
	int q_Bx, q_By, q_Bz;
	// global indices
	int q_rho_global, q_sx_global, q_sy_global, q_sz_global, q_Eges_global, q_Eadd_global;
	int q_Bx_global, q_By_global, q_Bz_global;
	int aux_energyType;
	bool use_dualEnergy, have_MagField;
	std::string name;
};

#endif
