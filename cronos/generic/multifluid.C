#include "multifluid.H"
#include <cmath>
#include <iostream>
#include <stdio.h>

using namespace std;

CronosMultifluid::CronosMultifluid(int numFluids, int energetics,
		bool apply_userFields) {
	this->numFluids = numFluids; // default value is single fluid
	this->energetics = energetics;
	fluidTypes.resize(numFluids);
	N_OMINTs.resize(numFluids);
	N_OMINTsUser.resize(numFluids);
	userFields.resize(numFluids);
//	q_rho.resize(numFluids);
//	q_sx.resize(numFluids);
//	q_sy.resize(numFluids);
//	q_sz.resize(numFluids);
//	q_Eges.resize(numFluids);
//	q_Eadd.resize(numFluids);

	// Make array of fluid classes
	fluids = new CronosFluid[numFluids];

	has_magField = false;
	i_magFluid = -1;
//	q_Bx = -1;
//	cout << " done " << endl;
//	q_By = -1;
//	q_Bz = -1;


	N_OMINT_all = 0;
	n_om_total= 0;

	n_omIntUser_total = 0;
	n_omUser_total = 0;

	this->apply_userFields = apply_userFields;
	this->use_dualEnergy = false;
}


void CronosMultifluid::set_dualEnergy(int aux_Energy) {
	this->use_dualEnergy = true;
	this->aux_energyType = aux_Energy;
}

void CronosMultifluid::unset_dualEnergy() {
	this->use_dualEnergy = false;
}




void CronosMultifluid::compute_Variables(int n_add, int n_subs) {
	//! To be applied after(!) user setup -- computes all dependent variables

	int allIndices=0;
	int allIndicesUser=0;
	n_omIntUser_total = 0;
	has_magField = false;

	int lognum = static_cast<int>(log10(numFluids)+1);
	char cnum[255];
	sprintf(cnum,"%i",lognum);
	std::string format = "%";
	format += cnum;
	format += "i";
	char fluidNum[255];

	// Loop over all fluids - setup of individual fluids
	for(int iFluid=0; iFluid<numFluids; ++iFluid) {

		// Setup of fluid
		if(fluidTypes(iFluid) == CRONOS_MHD) {
			cout << " Making mag with " << has_magField << " " << fluidTypes(iFluid) << endl;
			if(!has_magField) {

				if(use_dualEnergy) {
					fluids[iFluid].set_dualEnergy(aux_energyType);
				} else {
					fluids[iFluid].unset_dualEnergy();
				}

				fluids[iFluid].setup(fluidTypes(iFluid),energetics,n_add, n_subs,
						userFields(iFluid), allIndices, allIndicesUser);

				i_magFluid = iFluid;
				has_magField = true;
			} else {
				// If magnetic field has already been stored - use only hydro field in
				if(use_dualEnergy) {
					fluids[iFluid].set_dualEnergy(aux_energyType);
				} else {
					fluids[iFluid].unset_dualEnergy();
				}

				// ref to magnetic field
				fluids[iFluid].setup(CRONOS_HYDRO,energetics,n_add, n_subs,
						userFields(iFluid), allIndices, allIndicesUser);

				fluids[iFluid].add_Ref("B_x",get_q_Bx());
				fluids[iFluid].add_Ref("B_y",get_q_By());
				fluids[iFluid].add_Ref("B_z",get_q_Bz());
			}
		} else {
			// Add hydro fields
			if(use_dualEnergy) {
				fluids[iFluid].set_dualEnergy(aux_energyType);
			} else {
				fluids[iFluid].unset_dualEnergy();
			}

			fluids[iFluid].setup(fluidTypes(iFluid),energetics,n_add, n_subs,
					userFields(iFluid), allIndices, allIndicesUser);

		}

		// Set name of fluid:
		sprintf(fluidNum, format.c_str(), iFluid);
		string fluidName = "fluid";
		fluidName += fluidNum;
		fluids[iFluid].set_Name(fluidName);

		N_OMINTs(iFluid) = fluids[iFluid].get_N_OMINT();
		allIndices += N_OMINTs(iFluid);

		// add up user fields
		N_OMINTsUser(iFluid) = fluids[iFluid].get_N_OMINT_USER();
		allIndicesUser += N_OMINTsUser(iFluid);
		n_omIntUser_total += userFields(iFluid);
	}


	N_OMINT_all = allIndices + allIndicesUser;
	// Total number of fields to be integrated
	n_omInt_total = allIndices;
	n_Omega_total = n_omInt_total + n_add + n_subs;
	n_om_total = 2*n_omInt_total + n_add + n_subs;

	n_omUser_total = 2*n_omIntUser_total;

	// Store relation of local index to global index
	index_local.resize(N_OMINT_all);
	fluidIndex.resize(N_OMINT_all);
	int iFieldGlobal=0;
	for(int iFluid=0; iFluid<numFluids; ++iFluid) {
		// Loop over all fields of fluid
		for(int iField=0; iField<N_OMINTs(iFluid); ++iField) {
			index_local(iFieldGlobal) = iField;
			fluidIndex(iFieldGlobal) = iFluid;
			iFieldGlobal++;
		}
	}

}

int CronosMultifluid::get_IndexLocal(int iFieldGlobal) const {
	assert(iFieldGlobal>=0 && iFieldGlobal<N_OMINT_all);
	return index_local(iFieldGlobal);
}

int CronosMultifluid::get_FluidIndex(int iFieldGlobal) const {
	assert(iFieldGlobal>=0 && iFieldGlobal<N_OMINT_all);
	return fluidIndex(iFieldGlobal);
}

int CronosMultifluid::get_fluidType(int iFluid) const {
	return fluidTypes(iFluid);
}

int CronosMultifluid::get_q_rho(int iFluid) const {
//	return q_rho(iFluid);
	return fluids[iFluid].get_q_rho_global();
}

int CronosMultifluid::get_q_sx(int iFluid) const {
	return fluids[iFluid].get_q_sx_global();
}

int CronosMultifluid::get_q_sy(int iFluid) const {
	return fluids[iFluid].get_q_sy_global();
}

int CronosMultifluid::get_q_sz(int iFluid) const {
	return fluids[iFluid].get_q_sz_global();
}

int CronosMultifluid::get_q_Bx(int iFluid) const {
	if( fluidTypes(iFluid) == CRONOS_MHD ) {
		return fluids[i_magFluid].get_q_Bx_global();
	} else {
		return fluids[iFluid].get_q_Bx_global();
	}
}

int CronosMultifluid::get_q_By(int iFluid) const {
	if( fluidTypes(iFluid) == CRONOS_MHD ) {
		return fluids[i_magFluid].get_q_By_global();
	} else {
		return fluids[iFluid].get_q_By_global();
	}
//	return q_By(iFluid);
}

int CronosMultifluid::get_q_Bz(int iFluid) const {
	if( fluidTypes(iFluid) == CRONOS_MHD ) {
		return fluids[i_magFluid].get_q_Bz_global();
	} else {
		return fluids[iFluid].get_q_Bz_global();
	}
//	return q_Bz(iFluid);
}


int CronosMultifluid::get_q_Bx() const {
	if(has_magField) {
		return fluids[i_magFluid].get_q_Bx_global();
	} else {
		return -42;
	}
}

int CronosMultifluid::get_q_By() const {
	if(has_magField) {
		return fluids[i_magFluid].get_q_By_global();
	} else {
		return -42;
	}
}

int CronosMultifluid::get_q_Bz() const {
	if(has_magField) {
		return fluids[i_magFluid].get_q_Bz_global();
	} else {
		return -42;
	}
}

bool CronosMultifluid::with_magField() {
	return has_magField;
}


int CronosMultifluid::get_q_Eges(int iFluid) const {
	return fluids[iFluid].get_q_Eges_global();
}

int CronosMultifluid::get_q_Eadd(int iFluid) const {
	return fluids[iFluid].get_q_Eadd_global();
}

int CronosMultifluid::get_i_magFluid() const {
	if(has_magField) {
		return i_magFluid;
	} else {
		return -1;
	}
}

int CronosMultifluid::get_numFluids() const {
	//! Return number of fluids
	return numFluids;
}

int CronosMultifluid::get_N_OMINT() const {
//	return N_OMINT_all;
	return n_omInt_total;
}

int CronosMultifluid::get_N_OMEGA() const {
//	return n_Omega;
	return n_Omega_total;
}

int CronosMultifluid::get_N_OM() const {
//	return N_OM_all;
	return n_om_total;
}

int CronosMultifluid::get_N_OMINT_USER() const {
	return n_omIntUser_total;
}

int CronosMultifluid::get_N_OMEGA_USER() const {
	return n_omIntUser_total;
}

int CronosMultifluid::get_N_OM_USER() const {
	return n_omUser_total;
}

int CronosMultifluid::get_N_OMINT_ALL() const {
//	return N_OMINT_all + N_OMINT_USER_all;
	return n_omInt_total + n_omIntUser_total;
}

int CronosMultifluid::get_N_OMINT(int iFluid) const {
	assert(iFluid >= 0 && iFluid < numFluids);
	return N_OMINTs(iFluid);
}
