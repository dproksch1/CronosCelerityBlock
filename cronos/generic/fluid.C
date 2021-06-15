#include "fluid.H"
#include "cronos.H"
#include <iostream>
#include <stdlib.h>

using namespace std;

// Class holding properties of cronos fluids


CronosFluid::CronosFluid() {
	fluid_type = -1;
	this->numFields = 0;
	this->numVirtuals = 0;
	this->n_omInt = 0;
	this->n_omIntUser = 0;
	this->n_om = 0;
	have_MagField = false;
	this->use_dualEnergy = false;
	name = "fluid00";
}

CronosFluid::CronosFluid(int type, int energetics,
		int n_add, int n_subs,
		int num_user, int first_index, int first_index_user) {
	//! Constructor of fluid class
	/*!
	 * \param type: type of fluid -> hydro(1), MHD(2) or other(3)
	 * \param energetics: isothermal(1), full energy equation(2)
	 * \param first_index: starting index for field (default: 0)
	 * \param num_user: number of user fields (default: 0)
	 * \param use_dualEnergy: check for dual energy option (default: false)
	 * */

	this->numFields = 0;
	this->numVirtuals = 0;
	this->n_omInt = 0;
	this->n_omIntUser = 0;
	this->n_om = 0;
	have_MagField = false;
	this->use_dualEnergy = false;

	setup(type, energetics, n_add, n_subs, first_index,
			num_user);


}


void CronosFluid::setup(int type, int energetics, int n_add, int n_sub,
		int num_user, int first_index, int first_index_user) {
	//! Setup for fluid class
	/*!
	 * \param type: type of fluid -> hydro(1), MHD(2) or other(3)
	 * \param energetics: isothermal(1), full energy equation(2)
	 * \param first_index: starting index for field (default: 0)
	 * \param num_user: number of user fields (default: 0)
	 * \param use_dualEnergy: check for dual energy option (default: false)
	 * */

	this->fluid_type = type;
	this->energetics = energetics;
	this->n_add = n_add;
	this->n_subs = n_sub;
	this->first_index = first_index;
	this->first_index_user = first_index_user;
	this->num_user = num_user;
	this->numFields = 0;
	this->numVirtuals = 0;
	have_MagField = false;
	name = "fluid00";

	this->n_omInt= 0; // Number of variables to integrate
	this->n_om = 0;

	// Make field identifiers for standard fields
	if(type>0 && type<3) {
		compute_Variables();
	}
}


void CronosFluid::set_dualEnergy(int aux_Energy) {
	this->use_dualEnergy = true;
	this->aux_energyType = aux_Energy;
}

void CronosFluid::unset_dualEnergy() {
	this->use_dualEnergy = false;
}


void CronosFluid::add_Ref(std::string type, int field_index) {
	//! Add reference to other field
	/*! This routine is exclusively for multifluid calculations.
	 * This purpose is to store a reference to a field in another fluid.
	 * /param type: name of field
	 * /param field_index: index of other fluid
	 */
	if(type=="rho") {
		q_rho = field_index;
	} else if (type=="s_x") {
		q_sx = field_index;
	} else if (type=="s_y") {
		q_sy = field_index;
	} else if (type=="s_z") {
		q_sz = field_index;
	} else if (type=="B_x") {
		q_Bx = field_index;
	} else if (type=="B_y") {
		q_By = field_index;
	} else if (type=="B_z") {
		q_Bz = field_index;
	} else if (type=="Eges") {
		q_Eges = field_index;
	} else if (type=="Eadd") {
		q_Eadd = field_index;
	} else {
		std::cerr << " CronosFluid::Error: not applicable for user fields ";
		std::cerr << std::endl;
		exit(-29);
	}

	// Indicate that fields are only references
	is_virtual.push_back(true);
	numVirtuals++;

}


int CronosFluid::add_Field(std::string type, int my_index) {
	//! Add specific field to fluid
	/*! \param type: of field
	 * \param my_index: index where field should be stored
	 * */

	if(numFields==0) {
		first_index = my_index;
	}
	if(numFieldsUser==0) {
		first_index_user = my_index;
	}

	bool is_userField(false);
	if(type=="rho") {
		q_rho = (my_index++);
		fieldNames.push_back("rho");
	} else if (type=="s_x") {
		q_sx = (my_index++);
		fieldNames.push_back("v_x");
	} else if (type=="s_y") {
		q_sy = (my_index++);
		fieldNames.push_back("v_y");
	} else if (type=="s_z") {
		q_sz = (my_index++);
		fieldNames.push_back("v_z");
	} else if (type=="B_x") {
		q_Bx = (my_index++);
		have_MagField = true;
		fieldNames.push_back("B_x");
	} else if (type=="B_y") {
		q_By = (my_index++);
		have_MagField = true;
		fieldNames.push_back("B_y");
	} else if (type=="B_z") {
		q_Bz = (my_index++);
		have_MagField = true;
		fieldNames.push_back("B_z");
	} else if (type=="Eges") {
		q_Eges = (my_index++);
		fieldNames.push_back("Etherm");
	} else if (type=="Eadd") {
		q_Eadd = (my_index++);
#if(AUX_ENERGY == ENTROPY)
		fieldNames.push_back("Entropy");
#else
		fieldNames.push_back("Etherm");
#endif
	} else {
		// Adding a user field
		my_index++;
		is_userField = true;
	}

	// Indicate that fields are stored locally
	is_virtual.push_back(false);

	if(is_userField) {
		n_omIntUser++;
	} else {
		n_omInt++;
	}

	numFields++;
	return my_index;
}

void CronosFluid::compute_Variables() {
	//! Compute all variables for fluid
	if(energetics == ISOTHERMAL) {
		if(fluid_type == CRONOS_MHD) {
			n_omInt = 7;
		} else if (fluid_type == CRONOS_HYDRO) {
			n_omInt = 4;
		}
	} else if (energetics == FULL) {
		if(fluid_type == CRONOS_MHD) {
			n_omInt = 8;
		} else if (fluid_type == CRONOS_HYDRO) {
			n_omInt = 5;
		}
	}
	// If dual energy is used, add one additional field to each fluid
	if(use_dualEnergy) {
		n_omInt += 1;
	}
	// Add user fields
//	n_omInt += num_user;


	// Now set q-value for different fluids
	numFields = 0;
	q_rho = (numFields++);
	q_rho_global = q_rho + first_index;
	fieldNames.push_back("rho");
	q_sx  = (numFields++);
	q_sx_global = q_sx + first_index;
	fieldNames.push_back("v_x");
	q_sy  = (numFields++);
	q_sy_global = q_sy + first_index;
	fieldNames.push_back("v_y");
	q_sz  = (numFields++);
	q_sz_global = q_sz + first_index;
	fieldNames.push_back("v_z");


	if(fluid_type == CRONOS_MHD) {
		q_Bx = (numFields++);
		q_Bx_global = q_Bx + first_index;
		fieldNames.push_back("B_x");
		q_By = (numFields++);
		q_By_global = q_By + first_index;
		fieldNames.push_back("B_y");
		q_Bz = (numFields++);
		q_Bz_global = q_Bz + first_index;
		fieldNames.push_back("B_z");
		have_MagField = true;
	} else {
		q_Bx = -42;
		q_By = -42;
		q_Bz = -42;
	}

	if(energetics == FULL) {
		q_Eges = (numFields++);
		q_Eges_global = q_Eges + first_index;
		fieldNames.push_back("Etherm");
	}

	if(use_dualEnergy) {
		q_Eadd = (numFields++);
		q_Eadd_global = q_Eadd + first_index;
		cout << aux_energyType << " " << ENTROPY << endl;
//#if(AUX_ENERGY == ENTROPY)
		if(aux_energyType == 2) {
			fieldNames.push_back("Entropy");
		} else {
//#else
			fieldNames.push_back("Etherm");
		}
//#endif
	} else {
		q_Eadd = -42;
	}

	// Now allow for user fields
	numFieldsUser = num_user;
	n_omIntUser = numFieldsUser;

	if(n_omInt != numFields) {
		std::cerr << " Error: something is wrong in fluid.C " << std::endl;
		std::cerr << n_omInt << " " << numFields << std::endl;
		exit(3);
	}

	compute_Numbers();
	set_GridTypes();
	set_IndexGlobal();
}

void CronosFluid::compute_Numbers() {
	//! Compute total number of different arrays for fluid

	n_Omega = n_omInt + n_add + n_subs;
	n_omIntAll = n_omInt + n_omIntUser;
//	n_om = 2*n_omInt + n_add + n_subs;
	n_om = n_omInt + n_add + n_subs;
//	cout <<  " n_om " << n_omInt << " " << n_add << " " << n_subs << endl;
	n_omUser = 2*n_omIntUser;

}

void CronosFluid::set_GridTypes() {
	//! Store grid type storage
	/*! Prepare array holding positions on the grid
	 * Options are:
	 * 0 -> cell centered
	 * 1 -> x-face centered
	 * 2 -> y-face centered
	 * 3 -> z-face centered
	 */
	pos_type.resize(numFields+first_index);
	pos_type(q_rho) = 0;
	pos_type(q_sx) = 0;
	pos_type(q_sy) = 0;
	pos_type(q_sz) = 0;
	if(fluid_type == CRONOS_MHD) {
		pos_type(q_Bx) = 1;
		pos_type(q_By) = 2;
		pos_type(q_Bz) = 3;
	}
	if (energetics == FULL) {
		pos_type(q_Eges) = 0;
	}
	if(use_dualEnergy) {
		pos_type(q_Eadd) = 0;
	}
	// Indices of user fields
	for(int q_index=numFields; q_index<n_omInt+first_index; ++q_index) {
		pos_type(q_index) = 0;
	}
}


void CronosFluid::set_IndexGlobal() {
	//! Make list of global indices (only useful for multifluid)
	index_global.resize(numFields);
	for(int iField=0; iField<numFields; ++iField) {
		index_global(iField) = iField+first_index;
	}

	// Do the same for user fields
	index_global_user.resize(numFieldsUser);
	for(int iField=0; iField<numFieldsUser; ++iField) {
		index_global_user(iField) = iField+first_index_user;
	}
}

int CronosFluid::get_IndexGlobal(int iField) const {
	//! Output of global index from local index
	assert(iField>=0 && iField<numFields);
	return index_global(iField);
}

int CronosFluid::get_IndexGlobalUser(int iField) const {
	//! Output of global index from local index
	assert(iField>=0 && iField<numFieldsUser);
	return index_global_user(iField);
}

int CronosFluid::get_q_rho() const {
	//! Return index of density field
	return q_rho;
}

int CronosFluid::get_q_rho_global() const {
	//! Return index of density field
	return q_rho_global;
}

int CronosFluid::get_q_sx() const {
	//! Return index of x-component of velocity
	return q_sx;
}

int CronosFluid::get_q_sy() const {
	//! Return index of y-component of velocity
	return q_sy;
}

int CronosFluid::get_q_sz() const {
	//! Return index of z-component of velocity
	return q_sz;
}

int CronosFluid::get_q_sx_global() const {
	//! Return index of x-component of velocity
	return q_sx_global;
}

int CronosFluid::get_q_sy_global() const {
	//! Return index of y-component of velocity
	return q_sy_global;
}

int CronosFluid::get_q_sz_global() const {
	//! Return index of z-component of velocity
	return q_sz_global;
}

int CronosFluid::get_q_Bx() const {
	//! Return index of x-component of magnetic field
	return q_Bx;
}

int CronosFluid::get_q_By() const {
	//! Return index of y-component of magnetic field
	return q_By;
}

int CronosFluid::get_q_Bz() const {
	//! Return index of z-component of magnetic field
	return q_Bz;
}

int CronosFluid::get_q_Bx_global() const {
	//! Return index of x-component of magnetic field
	return q_Bx_global;
}

int CronosFluid::get_q_By_global() const {
	//! Return index of y-component of magnetic field
	return q_By_global;
}

int CronosFluid::get_q_Bz_global() const {
	//! Return index of z-component of magnetic field
	return q_Bz_global;
}

int CronosFluid::get_q_Eges() const {
	//! Return index of total energy
	return q_Eges;
}

int CronosFluid::get_q_Eges_global() const {
	//! Return index of total energy
	return q_Eges_global;
}

int CronosFluid::get_q_Eadd() const {
	//! Return index of additional energy field
	return q_Eadd;
}

int CronosFluid::get_q_Eadd_global() const {
	//! Return index of additional energy field
	return q_Eadd_global;
}

int CronosFluid::get_N_OMINT() const {
	return n_omInt;
}

int CronosFluid::get_N_OMEGA() {
	return n_Omega;
}

int CronosFluid::get_N_OM() {
	return n_om;
}

int CronosFluid::get_N_OMINT_ALL() const {
	return n_omIntAll;
}

int CronosFluid::get_N_OMINT_USER() {
	return num_user;
}

int CronosFluid::get_N_OM_USER() {
	return n_omUser;
}

int CronosFluid::get_fluid_type() const {
	//! Get type of fluid (hydro or magnetic)
	return fluid_type;
}

bool CronosFluid::has_MagField() const {
	//! Return wheter magnetic field is included in local field
	return have_MagField;
}

void CronosFluid::set_Name(std::string fluidName) {
	//! Set a name for the fluid
	name = fluidName;
}

std::string CronosFluid::get_Name() const {
	//! Get name of fluid
	return name;
}

std::string CronosFluid::get_fieldName(int iField) const {
	//! Get name of field
	return fieldNames[iField];
}

//N_OMINT_USER_all += userFields(iFluid);
//
