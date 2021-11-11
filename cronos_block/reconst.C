#include "reconst.H"
#include "reconst_const.H"
#include "reconst_2nd.H"
#include "reconst_WENO.H"
#include "CException.H"
#include "queue.H"

#include <stdlib.h>
#include <iostream>

using namespace std;


SingleReconstruction(const Data& gdata, int dir, int substep) :
        SingleReconstruction(gdata, dir, substep)
{
    int limiterID = -1;

    if (value_exists("Limiter_" + ToString(substep)))
    {
        limiterID = static_cast<int>(value("Limiter_" + ToString(substep)));
    }
    else
    {
        limiterID = static_cast<int>(value("Limiter"));
    }

    //Limiter = new limiter(limiterID);
    Limiter = limiter();

    if (gdata.rank == 0) {
        cout << "  Using Reconst: " + Limiter.get_Name() + " - " << "GENERIC" << " - dir " << dir << " - substep " << substep << endl;
    }

}

SingleReconstruction(const Data& gdata, const CronosFluid& fluid, int dir, int qReconst, int substep) :
		SingleReconstruction(gdata, fluid, dir, qReconst, substep)
{
    int limiterID = -1;
    string fieldName = fluid.get_fieldName(qReconst);


    if (value_exists(fieldName + "_Limiter_" + ToString(substep)) && substep > -1)
    {
        limiterID = static_cast<int>(value(fieldName + "_Limiter_" + ToString(substep)));
    }
    else if (value_exists(fieldName + "_Limiter"))
    {
        limiterID = static_cast<int>(value(fieldName + "_Limiter"));
    }
    else if (value_exists("Limiter_" + ToString(substep)) && substep > -1)
    {
        limiterID = static_cast<int>(value("Limiter_" + ToString(substep)));
    }
    else
    {
        limiterID = static_cast<int>(value("Limiter"));
    }

    //Limiter = new limiter(limiterID);
    //Limiter = limiter<minmod>(-1);
    Limiter = limiter();


    if (gdata.rank == 0) {
        cout << "  Using Reconst: " + Limiter.get_Name() + " - " << fluid.get_fieldName(qReconst) << " - dir " << dir << " - substep " << substep << endl;
    }
}

SingleReconstruction::~SingleReconstruction () {
}

/**
 * Compute limited version of derivatives for block-structured version of code
 * */
void SingleReconstruction::prepareDerivs(const Data& gdata, int ix, int iy, int iz) {
    //! Compute derivate from at given position

    getDerivs(gdata, ix, iy, iz);

    deriv_x = Limiter.compute(dudxp_q, dudx0_q, dudxm_q);
    deriv_y = Limiter.compute(dudyp_q, dudy0_q, dudym_q);
    deriv_z = Limiter.compute(dudzp_q, dudz0_q, dudzm_q);

    //deriv_x = Limiter.compute(dud_q[DudDir::_x]);
    //deriv_y = Limiter.compute(dud_q[DudDir::_y]);
    //deriv_z = Limiter.compute(dud_q[DudDir::_z]);
}

void SingleReconstruction::get_Vals_EW(const Data& gdata, phys_fields_0D& xFieldsW,
    phys_fields_0D &xFieldsE, int ix, int iy, int iz) {
        {

            // Take into account shifted collocation points
            REAL shift(0.);
#if (GEOM != CARTESIAN)
#if (SHIFTED_COLLOCATION == TRUE)
            shift = sources->shift_Geom_WE(gdata, i);
#endif
#endif

            int q = qReconst;

#if (NON_LINEAR_GRID == CRONOS_ON)
            REAL delx = gdata.getCen_dx(0, ix);
            xFieldsW.uPri(q) = gdata.om[q](ix, iy, iz) - (0.5 + shift) * deriv_x * delx;
            xFieldsE.uPri(q) = xFieldsW.uPri(q) + deriv_x * delx;
#else
            xFieldsW.uPri(q) = gdata.om[q](ix, iy, iz) - (0.5 + shift) * deriv_x;
            xFieldsE.uPri(q) = xFieldsW.uPri(q) + deriv_x;
#endif

        }
}

void SingleReconstruction::get_Vals_SN(const Data& gdata, phys_fields_0D& xFieldsS,
    phys_fields_0D &xFieldsN, int ix, int iy, int iz) {
        {

#if (NON_LINEAR_GRID == CRONOS_ON)
            REAL dely = gdata.getCen_dx(1, iy);	// iterate over all indices
#endif

            int q = qReconst;

#if (NON_LINEAR_GRID == CRONOS_ON)
            xFieldsS.uPri(q) = gdata.om[q](ix, iy, iz) - 0.5 * deriv_y * dely;
            xFieldsN.uPri(q) = xFieldsS.uPri(q) + deriv_y * dely;
#else
            xFieldsS.uPri(q) = gdata.om[q](ix, iy, iz) - 0.5 * deriv_y;
            xFieldsN.uPri(q) = xFieldsS.uPri(q) + deriv_y;
#endif
        }
}

void SingleReconstruction::get_Vals_BT(const Data& gdata, phys_fields_0D& xFieldsB, 
    phys_fields_0D &xFieldsT, int ix, int iy, int iz) {
        {

#if (NON_LINEAR_GRID == CRONOS_ON)
            REAL delz = gdata.getCen_dx(2, iz);	// iterate over all indices
#endif

            int q = qReconst;

#if (NON_LINEAR_GRID == CRONOS_ON)
            xFieldsB.uPri(q) = gdata.om[q](ix, iy, iz) - 0.5 * deriv_z * delz;
            xFieldsT.uPri(q) = xFieldsB.uPri(q) + deriv_z * delz;
#else
            xFieldsB.uPri(q) = gdata.om[q](ix, iy, iz) - 0.5 * deriv_z;
            xFieldsT.uPri(q) = xFieldsB.uPri(q) + deriv_z;
#endif

        }
}

Reconstruction::Reconstruction(const Data &gdata, int dir, const CronosFluid &fluid, int substep) {

	assert(dir>=0  && dir<=2);
	this->dir = dir;

	this->substep = substep;


#if (FLUID_TYPE == CRONOS_HYDRO)

	for(int q=0; q<N_OMINT; ++q) {
		ListNormal.push_back(q);
	}

//#elif (FLUID_TYPE == CRONOS_MHD)
#else

#if (FLUID_TYPE == CRONOS_MHD)
	int q_Bx = fluid.get_q_Bx();
	int q_By = fluid.get_q_By();
	int q_Bz = fluid.get_q_Bz();
	int n_omInt = fluid.get_N_OMINT();
//	int q_Bx = gdata.fluid.get_q_Bx();
//	int q_By = gdata.fluid.get_q_By();
//	int q_Bz = gdata.fluid.get_q_Bz();
//	int n_omInt = gdata.fluid.get_N_OMINT();
	bool has_magField = true;
#elif (FLUID_TYPE == CRONOS_MULTIFLUID)
	int q_Bx = gdata.fluids->get_q_Bx();
	int q_By = gdata.fluids->get_q_By();
	int q_Bz = gdata.fluids->get_q_Bz();
	int n_omInt = fluid.get_N_OMINT();
	bool has_magField = fluid.has_MagField();
#endif

	for(int q=0; q<n_omInt; ++q) {
		if(!has_magField || (q!=q_Bx && q!=q_By && q!=q_Bz)) {
			ListNormal.push_back(q);
		}
	}

	if(has_magField) {
		if(dir == 0) {
			ListParallel.push_back(q_Bx);
			ListPerp.push_back(q_By);
			ListPerp.push_back(q_Bz);
		} else if (dir == 1) {
			ListPerp.push_back(q_Bx);
			ListParallel.push_back(q_By);
			ListPerp.push_back(q_Bz);
		} else {
			ListPerp.push_back(q_Bx);
			ListPerp.push_back(q_By);
			ListParallel.push_back(q_Bz);
		}


	}

#endif
	
	set_singleReconstructions(gdata, fluid);
}

Reconstruction::Reconstruction(const Data &gdata, int dir, int num, int substep) {

	assert(dir >= 0 && dir <= 2);
	this->dir = dir;

	this->substep = substep;

	for(int q=0; q<num; ++q) {
		ListNormal.push_back(q);
	}
	
	set_singleReconstructions(gdata);
}


Reconstruction::~Reconstruction() {
}

void Reconstruction::compute(const Data& gdata, std::vector<phys_fields_0D> &allFields,
		int ix, int iy, int iz, Direction dir) {
//void Reconstruction::compute(const Data& gdata, phys_fields_0D allFields[],
//		int ix, int iy, int iz) {
	//! Reconstruction for block-structured code

	// Todo: Andere Generalisierung für punktweise rekonstruktion einführen?
	for(int q=0; q<ListNormal.size(); ++q) {
		
		// First compute all derivatives
		ListReconstructionNormal[q]->prepareDerivs(gdata, ix, iy, iz);

		if (dir == -1 || dir == DirX) {
			// Get Values in east and west direction:
			ListReconstructionNormal[q]->get_Vals_EW(gdata, allFields[0], allFields[1], ix, iy, iz);
		}

		if (dir == -1 || dir == DirY) {
			// Get Values in South and North directions:
			ListReconstructionNormal[q]->get_Vals_SN(gdata, allFields[2], allFields[3], ix, iy, iz);
		}

		if (dir == -1 || dir == DirZ) {
			// Get Values in Bottom and Top directions:
			ListReconstructionNormal[q]->get_Vals_BT(gdata, allFields[4], allFields[5], ix, iy, iz);
		}
	}

	return;

}