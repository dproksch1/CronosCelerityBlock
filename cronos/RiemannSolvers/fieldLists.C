#include "fieldLists.H"

//fieldLists::fieldLists(int dir) {
//	//! Constructur especially for user fields
//	assert(dir >= 0 && dir < 3);
//
//}

#if(FLUID_TYPE == CRONOS_MULTIFLUID)
//fieldLists::fieldLists(const CronosMultifluid &fluids, int dir) {
fieldLists::fieldLists(const CronosFluid &fluid, int dir) {
	assert(dir >= 0 && dir < 3);


	// Loop over all fluids
	int n_omInt = fluid.get_N_OMINT();

	if(fluid.get_fluid_type() == CRONOS_HYDRO) {

		for(int q=0; q<n_omInt; ++q) {
			GenType[q] = "Centered";
		}

	} else if(fluid.get_fluid_type() == CRONOS_MHD) {

		for(int q=0; q<n_omInt; ++q) {
			GenType[q] = "Centered";
		}
		// Check whether magnetic field is included in fluid
		if(fluid.has_MagField()) {
			// Find mag field components
			int q_Bx = fluid.get_q_Bx();
			int q_By = fluid.get_q_By();
			int q_Bz = fluid.get_q_Bz();

			// Depending on the direction set some fields to be different
			if(dir == 0) {
				GenType[q_Bx] = "Parallel";
				GenType[q_By] = "Perpendicular";
				GenType[q_Bz] = "Perpendicular";
			} else if (dir == 1) {
				GenType[q_Bx] = "Perpendicular";
				GenType[q_By] = "Parallel";
				GenType[q_Bz] = "Perpendicular";
			} else {
				GenType[q_Bx] = "Perpendicular";
				GenType[q_By] = "Perpendicular";
				GenType[q_Bz] = "Parallel";
			}
		}
	}

}



#else



fieldLists::fieldLists(const CronosFluid &fluid, int dir) {
	assert(dir >= 0 && dir < 3);

	int n_omInt = fluid.get_N_OMINT();
#if (FLUID_TYPE == CRONOS_HYDRO)

		for(int q=0; q<n_omInt; ++q) {
			GenType[q] = "Centered";
		}

#elif (FLUID_TYPE == CRONOS_MHD)
		// Find mag field components
		int q_Bx = fluid.get_q_Bx();
		int q_By = fluid.get_q_By();
		int q_Bz = fluid.get_q_Bz();

		for(int q=0; q<N_OMINT; ++q) {
			GenType[q] = "Centered";
		}

		// Depending on the direction set some fields to be different
		if(dir == 0) {
			GenType[q_Bx] = "Parallel";
			GenType[q_By] = "Perpendicular";
			GenType[q_Bz] = "Perpendicular";
		} else if (dir == 1) {
			GenType[q_Bx] = "Perpendicular";
			GenType[q_By] = "Parallel";
			GenType[q_Bz] = "Perpendicular";
		} else {
			GenType[q_Bx] = "Perpendicular";
			GenType[q_By] = "Perpendicular";
			GenType[q_Bz] = "Parallel";
		}
#endif
}
#endif


std::string fieldLists::getGenType(int q) {
	/* Return generic type of field 
	   Types are

	   Centered
	   Parallel
	   Perpendicular
	 */
	assert(q>=0 && q<N_OMINT);

	return GenType[q];
	
}


bool fieldLists::isGenType(int q, const std::string &type) {
	if(GenType[q] == type) {
		return true;
	} else {
		return false;
	}
}


bool fieldLists::isGenTypeCentered(int q) {
	if(GenType[q] == "Centered") {
		return true;
	} else {
		return false;
	}
}


bool fieldLists::isGenTypeParallel(int q) {
	if(GenType[q] == "Parallel") {
		return true;
	} else {
		return false;
	}
}


bool fieldLists::isGenTypePerpendicular(int q) {
	if(GenType[q] == "Perpendicular") {
		return true;
	} else {
		return false;
	}
}


//#if(FLUID_TYPE == CRONOS_MULTIFLUID)
////fieldLists::fieldLists(const CronosMultifluid &fluids, int dir) {
//fieldLists::fieldLists(const CronosFluid &fluid, int dir) {
//	assert(dir >= 0 && dir < 3);
//
//
//	// Loop over all fluids
//	int n_omInt = fluid.get_N_OMINT();
//
//		if(fluids.get_fluidType(iFluid) == CRONOS_HYDRO) {
//
//			for(int q=0; q<n_omInt; ++q) {
//				GenType[q] = "Centered";
//			}
//
//		} else if(fluids.get_fluidType(iFluid) == CRONOS_MHD) {
//
//			for(int q=0; q<n_omInt; ++q) {
//				GenType[q] = "Centered";
//			}
//			// Check wether magnetic field is included in fluid
//			if(iFluid == fluids.get_i_magFluid()) {
//				// Find mag field components
//				int q_Bx = fluids.get_q_Bx();
//				int q_By = fluids.get_q_By();
//				int q_Bz = fluids.get_q_Bz();
//
//				// Depending on the direction set some fields to be different
//				if(dir == 0) {
//					GenType[q_Bx] = "Parallel";
//					GenType[q_By] = "Perpendicular";
//					GenType[q_Bz] = "Perpendicular";
//				} else if (dir == 1) {
//					GenType[q_Bx] = "Perpendicular";
//					GenType[q_By] = "Parallel";
//					GenType[q_Bz] = "Perpendicular";
//				} else {
//					GenType[q_Bx] = "Perpendicular";
//					GenType[q_By] = "Perpendicular";
//					GenType[q_Bz] = "Parallel";
//				}
//			}
//		}
//
//	}
//
//
//}
//
