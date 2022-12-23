#include "reconst_block.H"
#include "CException.H"

#include <stdlib.h>
#include <iostream>

#if(STANDALONE_USAGE==FALSE)
#include "reconst.H"
#endif

using namespace std;


SingleReconstruction_Block::SingleReconstruction_Block(const Data& gdata, int dir, int substep)
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

SingleReconstruction_Block::SingleReconstruction_Block(const Data& gdata, const CronosFluid& fluid, int dir, int qReconst, int substep)
{
	derivPerp.resize(Index::set(-2),Index::set(gdata.mx[dir]+1));
	derivM.resize(Index::set(-2), Index::set(gdata.mx[dir] + 1));

	assert(dir>=0  && dir<=2);
	this->dir = dir;

	this->qReconst = qReconst;
	this->substep = substep;

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

    Limiter = limiter();


    if (gdata.rank == 0) {
        cout << "  Using Reconst: " + Limiter.get_Name() + " - " << fluid.get_fieldName(qReconst) << " - dir " << dir << " - substep " << substep << endl;
    }
}

SingleReconstruction_Block::~SingleReconstruction_Block () {
}

void SingleReconstruction_Block::getDerivs(const Data &gdata, int ix, int iy, int iz) {
	//! Compute derivate from at given position

	// Start with getting the properties of the cell
#if (NON_LINEAR_GRID == CRONOS_ON)

	double delxm = 0.5*(gdata.getCen_dx(0, ix-1) + gdata.getCen_dx(0, ix  ));
	double delxp = 0.5*(gdata.getCen_dx(0, ix  ) + gdata.getCen_dx(0, ix+1));
	double delx0 = delxm + delxp;

	double delym = 0.5*(gdata.getCen_dx(1, iy-1) + gdata.getCen_dx(1, iy  ));
	double delyp = 0.5*(gdata.getCen_dx(1, iy  ) + gdata.getCen_dx(1, iy+1));
	double dely0 = delym + delyp;

	double delzm = 0.5*(gdata.getCen_dx(2, iz-1) + gdata.getCen_dx(2, iz  ));
	double delzp = 0.5*(gdata.getCen_dx(2, iz  ) + gdata.getCen_dx(2, iz+1));
	double delz0 = delzm + delzp;
#endif

	int q = qReconst;
#if (NON_LINEAR_GRID == CRONOS_ON)

	dudxm_q = (gdata.om[q](ix  ,iy,iz) - gdata.om[q](ix-1,iy,iz))/delxm;
	dudx0_q = (-sqr(delxp)*gdata.om[q](ix-1,iy,iz) +
			(sqr(delxp) - sqr(delxm))*gdata.om[q](ix,iy,iz) +
			sqr(delxm)*gdata.om[q](ix+1,iy,iz))/(delxp*delxm*(delx0));
	dudxp_q = (gdata.om[q](ix+1,iy,iz) - gdata.om[q](ix  ,iy,iz))/delxp;

	dudym_q = (gdata.om[q](ix,iy  ,iz) - gdata.om[q](ix,iy-1,iz))/delym;
	dudy0_q = (-sqr(delyp)*gdata.om[q](iy,iy-1,iz) +
			(sqr(delyp) - sqr(delym))*gdata.om[q](ix,iy,iz) +
			sqr(delym)*gdata.om[q](iy,iy+1,iz))/(delyp*delym*(dely0));
	dudyp_q = (gdata.om[q](ix,iy+1,iz) - gdata.om[q](ix,iy  ,iz))/delyp;

	dudzm_q = (gdata.om[q](ix,iy,iz  ) - gdata.om[q](ix,iy,iz-1))/delzm;
	dudz0_q = (-sqr(delzp)*gdata.om[q](ix,iy,iz-1) +
			(sqr(delzp) - sqr(delzm))*gdata.om[q](ix,iy,iz) +
			sqr(delzm)*gdata.om[q](ix,iy,iz+1))/(delzp*delzm*(delz0));
	dudzp_q = (gdata.om[q](ix,iy,iz+1) - gdata.om[q](ix,iy,iz  ))/delzp;
#else

	dudxm_q =  gdata.om[q](ix  ,iy,iz) - gdata.om[q](ix-1,iy,iz);
//	double dudx0 = (gdata.om[q](ix+1,iy,iz) - gdata.om[q](ix-1,iy,iz))*0.5;
	dudxp_q =  gdata.om[q](ix+1,iy,iz) - gdata.om[q](ix  ,iy,iz);
	dudx0_q = 0.5*(dudxm_q+dudxp_q);

	dudym_q =  gdata.om[q](ix,iy  ,iz) - gdata.om[q](ix,iy-1,iz);
//	double dudy0 = (gdata.om[q](ix,iy+1,iz) - gdata.om[q](ix,iy-1,iz))*0.5;
	dudyp_q =  gdata.om[q](ix,iy+1,iz) - gdata.om[q](ix,iy  ,iz);
	dudy0_q = 0.5*(dudym_q + dudyp_q);

	dudzm_q =  gdata.om[q](ix,iy,iz  ) - gdata.om[q](ix,iy,iz-1);
//	double dudz0 = (gdata.om[q](ix,iy,iz+1) - gdata.om[q](ix,iy,iz-1))*0.5;
	dudzp_q =  gdata.om[q](ix,iy,iz+1) - gdata.om[q](ix,iy,iz  );
	dudz0_q = 0.5*(dudzm_q + dudzp_q);
#endif

}

/**
 * Compute limited version of derivatives for block-structured version of code
 * */
void SingleReconstruction_Block::prepareDerivs(const Data& gdata, int ix, int iy, int iz) {
    //! Compute derivate from at given position

    getDerivs(gdata, ix, iy, iz);

    deriv_x = Limiter.compute(dudxp_q, dudx0_q, dudxm_q);
    deriv_y = Limiter.compute(dudyp_q, dudy0_q, dudym_q);
    deriv_z = Limiter.compute(dudzp_q, dudz0_q, dudzm_q);

    //deriv_x = Limiter.compute(dud_q[DudDir::_x]);
    //deriv_y = Limiter.compute(dud_q[DudDir::_y]);
    //deriv_z = Limiter.compute(dud_q[DudDir::_z]);
}

void SingleReconstruction_Block::get_Vals_EW(const Data& gdata, phys_fields_0D& xFieldsW,
    phys_fields_0D &xFieldsE, int ix, int iy, int iz) {
        {

            // Take into account shifted collocation points
            double shift(0.);
#if (GEOM != CARTESIAN)
#if (SHIFTED_COLLOCATION == TRUE)
            shift = sources->shift_Geom_WE(gdata, i);
#endif
#endif

            int q = qReconst;

#if (NON_LINEAR_GRID == CRONOS_ON)
            double delx = gdata.getCen_dx(0, ix);
            xFieldsW.uPri(q) = gdata.om[q](ix, iy, iz) - (0.5 + shift) * deriv_x * delx;
            xFieldsE.uPri(q) = xFieldsW.uPri(q) + deriv_x * delx;
#else
            xFieldsW.uPri(q) = gdata.om[q](ix, iy, iz) - (0.5 + shift) * deriv_x;
            xFieldsE.uPri(q) = xFieldsW.uPri(q) + deriv_x;
#endif

        }
}

void SingleReconstruction_Block::get_Vals_SN(const Data& gdata, phys_fields_0D& xFieldsS,
    phys_fields_0D &xFieldsN, int ix, int iy, int iz) {
        {

#if (NON_LINEAR_GRID == CRONOS_ON)
            double dely = gdata.getCen_dx(1, iy);	// iterate over all indices
#endif

            int q = qReconst;

#if (NON_LINEAR_GRID == CRONOS_ON)
            xFieldsS.uPri(q) = gdata.om[q](ix, iy, iz) - 0.5 * deriv_y * dely;
            xFieldsN.uPri(q) = xFieldsS.uPri(q) + deriv_y * dely;
#else
            xFieldsS.uPri(q) = gdata.om[q](ix, iy, iz) - 0.5 * deriv_y;
            xFieldsN.uPri(q) = xFieldsS.uPri(q) + deriv_y;		

			/*if (q == 1) {
				cout << "uPri[FaceSouth][" << q << "] = " << xFieldsS.uPri(q) << endl;
				cout << "gdata: " << ix << "," << iy << ","<< iz << " " << gdata.om[q](ix, iy, iz) << endl;
				cout << "derivs: " << deriv_x << ", " << deriv_y << ","<< deriv_z << endl;
				cout << "uPri[FaceNorth][" << q << "] = " << xFieldsN.uPri(q) << endl;
			}*/
#endif
        }
}

void SingleReconstruction_Block::get_Vals_BT(const Data& gdata, phys_fields_0D& xFieldsB, 
    phys_fields_0D &xFieldsT, int ix, int iy, int iz) {
        {

#if (NON_LINEAR_GRID == CRONOS_ON)
            double delz = gdata.getCen_dx(2, iz);	// iterate over all indices
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

Reconstruction_Block::Reconstruction_Block(const Data &gdata, int dir, const CronosFluid &fluid, int substep) {

	assert(dir>=0  && dir<=2);
	this->dir = dir;

	this->substep = substep;


#if (FLUID_TYPE == CRONOS_HYDRO)

	for(int q=0; q<N_OMINT; ++q) {
		ListNormal.push_back(q);
	}

#endif
	
	set_singleReconstructions(gdata, fluid);
}

Reconstruction_Block::Reconstruction_Block(const Data &gdata, int dir, int num, int substep) {

	assert(dir >= 0 && dir <= 2);
	this->dir = dir;

	this->substep = substep;

	for(int q=0; q<num; ++q) {
		ListNormal.push_back(q);
	}
	
	set_singleReconstructions(gdata);
}


Reconstruction_Block::~Reconstruction_Block() {
}

void Reconstruction_Block::compute(const Data& gdata, std::vector<phys_fields_0D> &allFields,
		int ix, int iy, int iz, Direction dir) {

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

template <typename ... Ts>
SingleReconstruction_Block * get_reconst_block(string ch_reconst, const Ts & ... ts)
{
	SingleReconstruction_Block * reconst;

	/*if(ch_reconst == "WENO") {
		reconst = new SingleReconstruction_WENO(ts...);
	}
	else if(ch_reconst == "constant"){
		reconst = new SingleReconstruction_constant(ts...);
	}
	else*/ if(ch_reconst == "slope_limiter") {
		reconst = new SingleReconstruction_Block(ts...);
	}
	else {
	    throw CException("No such reconstruction available: " + ch_reconst);
	}

	return reconst;
}

#if(STANDALONE_USAGE==TRUE)
string read_reconst(string fieldName, const int substep) {
	string reconstName = "";

	if (value_exists(fieldName + "_Reconstruction_" + ToString(substep)) && substep > -1 && fieldName != "")
	{
		reconstName = svalue(fieldName + "_Reconstruction_" + ToString(substep));
	}
	else if(value_exists(fieldName + "_Reconstruction") && fieldName != "")
	{
		reconstName = svalue(fieldName + "_Reconstruction");
	}
	else if (value_exists("Reconstruction_" + ToString(substep)) && substep > -1)
	{
		reconstName = svalue("Reconstruction_" + ToString(substep));
	}
	else if(value_exists("Reconstruction"))
	{
		reconstName = svalue("Reconstruction");
	}
	else
	{
		reconstName = "slope_limiter";
	}

	return reconstName;
}
#endif

void Reconstruction_Block::set_singleReconstructions(const Data & gdata) {

	string ch_reconst = read_reconst("", substep);

	for(iter = ListNormal.begin(); iter != ListNormal.end(); ++iter) {
		ListReconstructionNormal.push_back(
				get_reconst_block(ch_reconst, gdata, dir, substep)
		);
	}

	for(iter = ListParallel.begin(); iter != ListParallel.end(); ++iter) {
		ListReconstructionPar.push_back(
				get_reconst_block(ch_reconst, gdata, dir, substep)
		);
	}

	for(iter = ListPerp.begin(); iter != ListPerp.end(); ++iter) {
		ListReconstructionPerp.push_back(
				get_reconst_block(ch_reconst, gdata, dir, substep)
		);
	}
}

void Reconstruction_Block::set_singleReconstructions(const Data & gdata, const CronosFluid &fluid) {

	string ch_reconst = "";

	for(iter = ListNormal.begin(); iter != ListNormal.end(); ++iter) {
		ch_reconst = read_reconst(fluid.get_fieldName(*iter), substep);
		ListReconstructionNormal.push_back(
				get_reconst_block(ch_reconst, gdata, fluid, dir, *iter, substep)
		);
	}

	for(iter = ListParallel.begin(); iter != ListParallel.end(); ++iter) {
		ch_reconst = read_reconst(fluid.get_fieldName(*iter), substep);
		ListReconstructionPar.push_back(
				get_reconst_block(ch_reconst, gdata, fluid, dir, *iter, substep)
		);
	}

	for(iter = ListPerp.begin(); iter != ListPerp.end(); ++iter) {
		ch_reconst = read_reconst(fluid.get_fieldName(*iter), substep);
		ListReconstructionPerp.push_back(
				get_reconst_block(ch_reconst, gdata, fluid, dir, *iter, substep)
		);
	}
}