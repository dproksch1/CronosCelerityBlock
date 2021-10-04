#include "reconst.H"
#include "reconst_const.H"
#include "reconst_2nd.H"
#include "reconst_WENO.H"
#include "CException.H"
#include "queue.H"

#include <stdlib.h>
#include <iostream>

using namespace std;


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



template <typename ... Ts>
SingleReconstruction * get_reconst(string ch_reconst, const Ts & ... ts)
{
	SingleReconstruction * reconst;

	/*if(ch_reconst == "WENO") {
		reconst = new SingleReconstruction_WENO(ts...);
	}
	else if(ch_reconst == "constant"){
		reconst = new SingleReconstruction_constant(ts...);
	}
	else*/ if(ch_reconst == "slope_limiter") {
		reconst = new SingleReconstruction_2nd(ts...);
	}
	else {
	    throw CException("No such reconstruction available: " + ch_reconst);
	}

	return reconst;
}

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

void Reconstruction::set_singleReconstructions(const Data & gdata) {

	string ch_reconst = read_reconst("", substep);

	for(iter = ListNormal.begin(); iter != ListNormal.end(); ++iter) {
		ListReconstructionNormal.push_back(
				get_reconst(ch_reconst, gdata, dir, substep)
		);
	}

	for(iter = ListParallel.begin(); iter != ListParallel.end(); ++iter) {
		ListReconstructionPar.push_back(
				get_reconst(ch_reconst, gdata, dir, substep)
		);
	}

	for(iter = ListPerp.begin(); iter != ListPerp.end(); ++iter) {
		ListReconstructionPerp.push_back(
				get_reconst(ch_reconst, gdata, dir, substep)
		);
	}
}

void Reconstruction::set_singleReconstructions(const Data & gdata, const CronosFluid &fluid) {

	string ch_reconst = "";

	for(iter = ListNormal.begin(); iter != ListNormal.end(); ++iter) {
		ch_reconst = read_reconst(fluid.get_fieldName(*iter), substep);
		ListReconstructionNormal.push_back(
				get_reconst(ch_reconst, gdata, fluid, dir, *iter, substep)
		);
	}

	for(iter = ListParallel.begin(); iter != ListParallel.end(); ++iter) {
		ch_reconst = read_reconst(fluid.get_fieldName(*iter), substep);
		ListReconstructionPar.push_back(
				get_reconst(ch_reconst, gdata, fluid, dir, *iter, substep)
		);
	}

	for(iter = ListPerp.begin(); iter != ListPerp.end(); ++iter) {
		ch_reconst = read_reconst(fluid.get_fieldName(*iter), substep);
		ListReconstructionPerp.push_back(
				get_reconst(ch_reconst, gdata, fluid, dir, *iter, substep)
		);
	}
}


void Reconstruction::compute(Queue& queue, const Data &gdata,
                             fields_1D &fields,
                             phys_fields_1D &physValL,
                             phys_fields_1D &physValR) {
	// Doing normal reconstruction
	int q=0;
	for(iter = ListNormal.begin(); iter != ListNormal.end(); ++iter) {
		ListReconstructionNormal[q]->computeNormal(
					queue,
					gdata,
					fields.omLocORIG[ListNormal[q]],
					physValL.uPriORIG[ListNormal[q]],
		            physValR.uPriORIG[ListNormal[q]],
		            fields.derivORIG[ListNormal[q]],
					physValL.uPriSYCL[ListNormal[q]],
		            physValR.uPriSYCL[ListNormal[q]],
					fields.omLocSYCL[ListNormal[q]],
					fields.derivSYCL[ListNormal[q]]);
		q++;
	}
  
	// Doing parallel reconstruction for special fields
	q=0;
	for(iter=ListParallel.begin(); iter != ListParallel.end(); ++iter) {
		ListReconstructionPar[q]->computePar(
					gdata,
					fields.omLocORIG[ListParallel[q]],
					physValL.uPriORIG[ListParallel[q]],
					physValR.uPriORIG[ListParallel[q]]);
		q++;
	}
  
	// Doing perpendicular reconstruction for special fields
	q=0;
	for(iter = ListPerp.begin(); iter != ListPerp.end(); ++iter) {
		ListReconstructionPerp[q]->computePerp(
					queue,
					gdata,
					fields.omLocPORIG[ListPerp[q]],
		            fields.omLocMORIG[ListPerp[q]],
		            physValL.uPriORIG[ListPerp[q]],
		            physValR.uPriORIG[ListPerp[q]],
		            fields.derivORIG[ListPerp[q]],
					fields.omLocPSYCL[ListPerp[q]],
					fields.omLocMSYCL[ListPerp[q]],
					fields.derivSYCL[ListPerp[q]]);
		q++;
	}
}

void Reconstruction::compute(const Data& gdata, std::vector<phys_fields_0D> &allFields,
		int ix, int iy, int iz, Direction dir) {
//void Reconstruction::compute(const Data& gdata, phys_fields_0D allFields[],
//		int ix, int iy, int iz) {
	//! Reconstruction for block-structured code

	// Todo: Andere Generalisierung für punktweise rekonstruktion einführen?
	for(int q=0; q<ListNormal.size(); ++q) {
		// First compute all derivatives
		ListReconstructionNormal[q]->perpareDerivs(gdata, ix, iy, iz);

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

void Reconstruction::compute(const Data &gdata,
                             NumMatrix<REAL,1> input[],
                             NumMatrix<REAL,1> lhs[],
                             NumMatrix<REAL,1> rhs[],
                             NumMatrix<REAL,1> deriv[]) {
	assert(false && "Not implemented");
	//// Doing normal reconstruction
	//int q=0;
	//for(iter = ListNormal.begin(); iter != ListNormal.end(); ++iter) {
	//	ListReconstructionNormal[q]->computeNormal(
	//				gdata,
	//	            input[ListNormal[q]],
	//				lhs[ListNormal[q]],
	//				rhs[ListNormal[q]],
	//	            deriv[ListNormal[q]]);
	//	q++;
	//}

	//// Doing parallel reconstruction for special fields
	//q=0;
	//for(iter=ListParallel.begin(); iter != ListParallel.end(); ++iter) {
	//}

	//// Doing perpendicular reconstruction for special fields
	//q=0;
	//for(iter = ListPerp.begin(); iter != ListPerp.end(); ++iter) {
	//}


}

//SingleReconstruction::SingleReconstruction() { }

SingleReconstruction::SingleReconstruction(const Data &gdata, const CronosFluid &fluid, int dir, int qReconst, int substep) {

	derivPerp.resize(Index::set(-2),Index::set(gdata.mx[dir]+1));
	derivM.resize(Index::set(-2), Index::set(gdata.mx[dir] + 1));

	//dudx = Buffer<REAL, 1>(Range<1>(DudIndex::Max_Index));
	//dud_q = std::vector<Buffer<REAL, 1>>(DudDir::Max_Dir, Range<1>(DudIndex::Max_Index));

	assert(dir>=0  && dir<=2);
	this->dir = dir;

	this->qReconst = qReconst;
	this->substep = substep;
}


SingleReconstruction::SingleReconstruction(const Data& gdata, int dir, int substep) : SingleReconstruction(gdata, {}, dir, -1, substep) {

	//derivPerp.resize(Index::set(-2),Index::set(gdata.mx[dir]+1));
	//derivM.resize(Index::set(-2),Index::set(gdata.mx[dir]+1));

	//dudx = Buffer<REAL, 1>(Range<1>(DudIndex::NUM_ELEMENTS));
	//dud_q = Buffer<REAL, 2>(Range<2>(DudDir::NUM_ELEMENTS, DudIndex::NUM_ELEMENTS));

	//assert(dir>=0  && dir<=2);
	//this->dir = dir;

	//this->qReconst = -1;
	//this->substep = substep;
}


SingleReconstruction::~SingleReconstruction() {
}

void SingleReconstruction::getDerivs(const Data &gdata, int ix, int iy, int iz) {
	//! Compute derivate from at given position

	// Start with getting the properties of the cell
#if (NON_LINEAR_GRID == CRONOS_ON)
	REAL delxm = 0.5*(gdata.getCen_dx(0, ix-1) + gdata.getCen_dx(0, ix  ));
	REAL delxp = 0.5*(gdata.getCen_dx(0, ix  ) + gdata.getCen_dx(0, ix+1));
	REAL delx0 = delxm + delxp;

	REAL delym = 0.5*(gdata.getCen_dx(1, iy-1) + gdata.getCen_dx(1, iy  ));
	REAL delyp = 0.5*(gdata.getCen_dx(1, iy  ) + gdata.getCen_dx(1, iy+1));
	REAL dely0 = delym + delyp;

	REAL delzm = 0.5*(gdata.getCen_dx(2, iz-1) + gdata.getCen_dx(2, iz  ));
	REAL delzp = 0.5*(gdata.getCen_dx(2, iz  ) + gdata.getCen_dx(2, iz+1));
	REAL delz0 = delzm + delzp;
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
	//auto dudx_acc = dud_q[_x].get_access<cl::sycl::access::mode::read_write>();
	//auto dudy_acc = dud_q[_y].get_access<cl::sycl::access::mode::read_write>();
	//auto dudz_acc = dud_q[_z].get_access<cl::sycl::access::mode::read_write>();
	//auto& dudxm_q = dudx_acc[SingleReconstruction::DudIndex::_m];
	//auto& dudxp_q = dudx_acc[SingleReconstruction::DudIndex::_p];
	//auto& dudx0_q = dudx_acc[SingleReconstruction::DudIndex::_0];
	//auto& dudym_q = dudy_acc[SingleReconstruction::DudIndex::_m];
	//auto& dudyp_q = dudy_acc[SingleReconstruction::DudIndex::_p];
	//auto& dudy0_q = dudy_acc[SingleReconstruction::DudIndex::_0];
	//auto& dudzm_q = dudz_acc[SingleReconstruction::DudIndex::_m];
	//auto& dudzp_q = dudz_acc[SingleReconstruction::DudIndex::_p];
	//auto& dudz0_q = dudz_acc[SingleReconstruction::DudIndex::_0];
	dudxm_q =  gdata.om[q](ix  ,iy,iz) - gdata.om[q](ix-1,iy,iz);
//	REAL dudx0 = (gdata.om[q](ix+1,iy,iz) - gdata.om[q](ix-1,iy,iz))*0.5;
	dudxp_q =  gdata.om[q](ix+1,iy,iz) - gdata.om[q](ix  ,iy,iz);
	dudx0_q = 0.5*(dudxm_q+dudxp_q);

	dudym_q =  gdata.om[q](ix,iy  ,iz) - gdata.om[q](ix,iy-1,iz);
//	REAL dudy0 = (gdata.om[q](ix,iy+1,iz) - gdata.om[q](ix,iy-1,iz))*0.5;
	dudyp_q =  gdata.om[q](ix,iy+1,iz) - gdata.om[q](ix,iy  ,iz);
	dudy0_q = 0.5*(dudym_q + dudyp_q);

	dudzm_q =  gdata.om[q](ix,iy,iz  ) - gdata.om[q](ix,iy,iz-1);
//	REAL dudz0 = (gdata.om[q](ix,iy,iz+1) - gdata.om[q](ix,iy,iz-1))*0.5;
	dudzp_q =  gdata.om[q](ix,iy,iz+1) - gdata.om[q](ix,iy,iz  );
	dudz0_q = 0.5*(dudzm_q + dudzp_q);
#endif

}

void SingleReconstruction::getDeriv(Queue& queue, const Data &gdata, const NumMatrix<REAL,1> &inputORIG, int iPos, NumMatrix<REAL,1>& dudxORIG, Buffer<REAL, 1> inputSYCL, Buffer<REAL, 1> dudxSYCL) {
#if (NON_LINEAR_GRID == CRONOS_ON)
		// First version -- without second order correction
		REAL delxm = 0.5*(gdata.getCen_dx(dir, iPos-1) +
		                  gdata.getCen_dx(dir, iPos  ));
		REAL delxp = 0.5*(gdata.getCen_dx(dir, iPos  ) +
		                  gdata.getCen_dx(dir, iPos+1));
		REAL delx0 = delxm + delxp;

		dudxp = (input(iPos+1) - input(iPos  ))/delxp;
//		dudx0 = (input(iPos+1) - input(iPos-1))/delx0;

		dudx0 = (-sqr(delxp)*input(iPos-1) + (sqr(delxp) - sqr(delxm))*input(iPos)
				+ sqr(delxm)*input(iPos+1))/(delxp*delxm*(delx0));

		dudxm = (input(iPos  ) - input(iPos-1))/delxm;
#else
		dudxORIG(_p) = inputORIG(iPos + 1) - inputORIG(iPos);
		dudxORIG(_0) = (inputORIG(iPos + 1) - inputORIG(iPos - 1)) * 0.5;
		dudxORIG(_m) = inputORIG(iPos) - inputORIG(iPos - 1);

		//int iPosSYCL = iPos + 3;
		//auto dudx_acc = dudxSYCL.get_access<cl::sycl::access::mode::discard_write>();
		//auto input_acc = inputSYCL.get_access<cl::sycl::access::mode::read>();
		//dudx_acc[_p] = input_acc[iPosSYCL + 1] - input_acc[iPosSYCL];
		//dudx_acc[_0] = (input_acc[iPosSYCL + 1] - input_acc[iPosSYCL - 1]) * 0.5;
		//dudx_acc[_m] = input_acc[iPosSYCL] - input_acc[iPosSYCL - 1];
#endif

}


Reconstruction2D::Reconstruction2D(const Data &gdata, const int &dir) {

	assert(dir >= 0 && dir < DIM);
	if(dir == 0) {
		this->dir0 = 1;
		this->dir1 = 2;
	} else if (dir == 1) {
		this->dir0 = 2;
		this->dir1 = 0;
	} else {
		this->dir0 = 0;
		this->dir1 = 1;
	}
}


void Reconstruction2D::compute(const Data &gdata,
                               fields_2D &fields,
                               phys_fields_2D &physValLL,
                               phys_fields_2D &physValLR,
                               phys_fields_2D &physValRL,
                               phys_fields_2D &physValRR) {
	// Doing normal reconstruction
	if(dir0 == 1 && dir1 == 2) {
		computeNormal(gdata, fields.v_y, fields.dvydx[0], fields.dvydx[1],
		              physValLL.v_y, physValLR.v_y, physValRL.v_y, physValRR.v_y);
		computeNormal(gdata, fields.v_z, fields.dvzdx[0], fields.dvzdx[1],
		              physValLL.v_z, physValLR.v_z, physValRL.v_z, physValRR.v_z);
    
		computePerp(gdata, fields.B_y, fields.dBydx, 2,
		            physValLL.B_y, physValLR.B_y, physValRL.B_y, physValRR.B_y);
		computePerp(gdata, fields.B_z, fields.dBzdx, 1,
		            physValLL.B_z, physValLR.B_z, physValRL.B_z, physValRR.B_z);
	} else if (dir0 == 2 && dir1 == 0) {
		computeNormal(gdata, fields.v_z, fields.dvzdx[0], fields.dvzdx[1],
		              physValLL.v_z, physValLR.v_z, physValRL.v_z, physValRR.v_z);
		computeNormal(gdata, fields.v_x, fields.dvxdx[0], fields.dvxdx[1],
		              physValLL.v_x, physValLR.v_x, physValRL.v_x, physValRR.v_x);
    
		computePerp(gdata, fields.B_z, fields.dBzdx, 0,
		            physValLL.B_z, physValLR.B_z, physValRL.B_z, physValRR.B_z);
		computePerp(gdata, fields.B_x, fields.dBxdx, 2,
		            physValLL.B_x, physValLR.B_x, physValRL.B_x, physValRR.B_x);
	} else if(dir0 == 0 && dir1 == 1) {
		computeNormal(gdata, fields.v_x, fields.dvxdx[0], fields.dvxdx[1],
		              physValLL.v_x, physValLR.v_x, physValRL.v_x, physValRR.v_x);
		computeNormal(gdata, fields.v_y, fields.dvydx[0], fields.dvydx[1],
		              physValLL.v_y, physValLR.v_y, physValRL.v_y, physValRR.v_y);
    
		computePerp(gdata, fields.B_x, fields.dBxdx, 1,
		            physValLL.B_x, physValLR.B_x, physValRL.B_x, physValRR.B_x);
		computePerp(gdata, fields.B_y, fields.dBydx, 0,
		            physValLL.B_y, physValLR.B_y, physValRL.B_y, physValRR.B_y);

	}
}






void Reconstruction2D::getDeriv(const Data &gdata, 
		const NumMatrix<REAL,2> &input,
		int dir, int i, int j) {
	if(dir == dir0) {

#if (NON_LINEAR_GRID == CRONOS_ON)
		// First version -- without second order correction
		REAL delxm = 0.5*(gdata.getCen_dx(dir, i-1) +
				gdata.getCen_dx(dir, i  ));
		REAL delxp = 0.5*(gdata.getCen_dx(dir, i+1) +
				gdata.getCen_dx(dir, i  ));
		REAL delx0 = delxm + delxp;
#endif


#if (NON_LINEAR_GRID == CRONOS_ON)
		dudxp = (input(i+1,j) - input(i  ,j))/delxp;
		dudx0 = (input(i+1,j) - input(i-1,j))/delx0;
		dudxm = (input(i  ,j) - input(i-1,j))/delxm;
#else
		dudxp =  input(i+1,j) - input(i  ,j);
		dudx0 = (input(i+1,j) - input(i-1,j))/2.;
		dudxm =  input(i  ,j) - input(i-1,j);
#endif

	} else if (dir == dir1) {

#if (NON_LINEAR_GRID == CRONOS_ON)
		REAL delxm = 0.5*(gdata.getCen_dx(dir, j-1) +
				gdata.getCen_dx(dir, j  ));
		REAL delxp = 0.5*(gdata.getCen_dx(dir, j+1) +
				gdata.getCen_dx(dir, j  ));
		REAL delx0 = delxm + delxp;

		dudyp = (input(i,j+1) - input(i,j  ))/delxp;
		dudy0 = (input(i,j+1) - input(i,j-1))/delx0;
		dudym = (input(i,j  ) - input(i,j-1))/delxm;
#else
		dudyp =  input(i,j+1) - input(i,j  );
		dudy0 = (input(i,j+1) - input(i,j-1))/2.;
		dudym =  input(i,j  ) - input(i,j-1);
#endif

	}
}
