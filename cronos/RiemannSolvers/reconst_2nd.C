#include "reconst_2nd.H"
#include <iostream>
#include <stdlib.h>

using namespace std;

//SingleReconstruction_2nd::SingleReconstruction_2nd(const Data &gdata, const CronosFluid &fluid, int dir, int qReconst, int substep):
//	SingleReconstruction(gdata, fluid, dir, qReconst, substep)
//{
//	int limiterID = -1;
//	string fieldName = fluid.get_fieldName(qReconst);
//
//
//	if (value_exists(fieldName + "_Limiter_" + ToString(substep)) && substep > -1)
//	{
//		limiterID = static_cast<int>(value(fieldName + "_Limiter_" + ToString(substep)));
//	}
//	else if(value_exists(fieldName + "_Limiter"))
//	{
//		limiterID = static_cast<int>(value(fieldName + "_Limiter"));
//	}
//	else if (value_exists("Limiter_" + ToString(substep)) && substep > -1)
//	{
//		limiterID = static_cast<int>(value("Limiter_" + ToString(substep)));
//	} else
//	{
//		limiterID = static_cast<int>(value("Limiter"));
//	}
//
//	Limiter = new limiter(limiterID);
//
//
//	if (gdata.rank == 0) {
//		cout << "  Using Reconst: " + Limiter->get_Name() + " - " << fluid.get_fieldName(qReconst) << " - dir " << dir << " - substep " << substep << endl;
//	}
//}


//SingleReconstruction_2nd::SingleReconstruction_2nd(const Data & gdata, int dir, int substep):
//	SingleReconstruction(gdata, dir, substep)
//{
//	int limiterID = -1;
//
//	if (value_exists("Limiter_" + ToString(substep)))
//	{
//		limiterID = static_cast<int>(value("Limiter_" + ToString(substep)));
//	} else
//	{
//		limiterID = static_cast<int>(value("Limiter"));
//	}
//
//	Limiter = new limiter(limiterID);
//
//	if (gdata.rank == 0) {
//		cout << "  Using Reconst: " + Limiter->get_Name() + " - " << "GENERIC" << " - dir " << dir << " - substep " << substep << endl;
//	}
//
//}

//SingleReconstruction_2nd::~SingleReconstruction_2nd() {
//	delete Limiter;
//}





//void SingleReconstruction_2nd::get_Vals_EW(const Data &gdata, phys_fields_0D &xFieldsW,
//		phys_fields_0D &xFieldsE, int ix, int iy, int iz)
//{
//
//	// Take into account shifted collocation points
//	REAL shift(0.);
//#if (GEOM != CARTESIAN)
//#if (SHIFTED_COLLOCATION == TRUE)
//	shift = sources->shift_Geom_WE(gdata, i);
//#endif
//#endif
//
//	int q = qReconst;
//
//#if (NON_LINEAR_GRID == CRONOS_ON)
//	REAL delx = gdata.getCen_dx(0, ix );
//	xFieldsW.uPri(q) = gdata.om[q](ix,iy,iz) - (0.5+shift)*deriv_x*delx;
//	xFieldsE.uPri(q) = xFieldsW.uPri(q) + deriv_x*delx;
//#else
//	xFieldsW.uPri(q) = gdata.om[q](ix,iy,iz) - (0.5+shift)*deriv_x;
//	xFieldsE.uPri(q) = xFieldsW.uPri(q) + deriv_x;
//#endif
//
//}
//
//void SingleReconstruction_2nd::get_Vals_SN(const Data &gdata, phys_fields_0D &xFieldsS,
//		phys_fields_0D &xFieldsN, int ix, int iy, int iz)
//{
//
//#if (NON_LINEAR_GRID == CRONOS_ON)
//	REAL dely = gdata.getCen_dx(1, iy );	// iterate over all indices
//#endif
//
//	int q = qReconst;
//
//#if (NON_LINEAR_GRID == CRONOS_ON)
//	xFieldsS.uPri(q) = gdata.om[q](ix,iy,iz) - 0.5*deriv_y*dely;
//	xFieldsN.uPri(q) = xFieldsS.uPri(q) + deriv_y*dely;
//#else
//	xFieldsS.uPri(q) = gdata.om[q](ix,iy,iz) - 0.5*deriv_y;
//	xFieldsN.uPri(q) = xFieldsS.uPri(q) + deriv_y;
//#endif
//}
//
//void SingleReconstruction_2nd::get_Vals_BT(const Data &gdata, phys_fields_0D &xFieldsB,
//		phys_fields_0D &xFieldsT, int ix, int iy, int iz)
//{
//
//#if (NON_LINEAR_GRID == CRONOS_ON)
//	REAL delz = gdata.getCen_dx(2, iz );	// iterate over all indices
//#endif
//
//	int q = qReconst;
//
//#if (NON_LINEAR_GRID == CRONOS_ON)
//	xFieldsB.uPri(q) = gdata.om[q](ix,iy,iz) - 0.5*deriv_z*delz;
//	xFieldsT.uPri(q) = xFieldsB.uPri(q) + deriv_z*delz;
//#else
//	xFieldsB.uPri(q) = gdata.om[q](ix,iy,iz) - 0.5*deriv_z;
//	xFieldsT.uPri(q) = xFieldsB.uPri(q) + deriv_z;
//#endif
//
//}



//void SingleReconstruction_2nd::prepareDerivs(const Data &gdata, int ix, int iy, int iz) {
//	//! Compute derivate from at given position
//
//	getDerivs(gdata, ix, iy, iz);
//
//	//deriv_x = Limiter->compute(dudxp_q, dudx0_q, dudxm_q);
//	//deriv_x = Limiter->compute(dudyp_q, dudy0_q, dudym_q);
//	//deriv_x = Limiter->compute(dudzp_q, dudz0_q, dudzm_q);
//
//	deriv_x = Limiter->compute(dud_q[DudDir::_x]);
//	deriv_y = Limiter->compute(dud_q[DudDir::_y]);
//	deriv_z = Limiter->compute(dud_q[DudDir::_z]);
//}
//
//void SingleReconstruction_2nd::getDeriv_limit(Queue& queue, const Data &gdata,
//		const NumMatrix<REAL,1> &inputORIG, NumMatrix<REAL,1> &derivORIG, Buffer<REAL, 1>& inputSYCL, Buffer<REAL, 1>& derivSYCL) {
//	queue.submit([&](Handler& cgh) {
//		auto dudx_acc = dudx.get_access<cl::sycl::access::mode::discard_write>();
//		auto derivSYCL_acc = derivSYCL.get_access<cl::sycl::access::mode::discard_write>();
//		cgh.single_task<class myKernel2>([=]() {
//			for (int i = -2; i <= gdata.mx[dir] + 1; ++i) {
//				//getDeriv(queue, gdata, inputORIG, i, inputSYCL);
//				auto dudxp = dudx_acc[_p];
//				auto dudx0 = dudx_acc[_p];
//				auto dudxm = dudx_acc[_p];
//				dudxp = inputORIG(i + 1) - inputORIG(i);
//				dudx0 = (inputORIG(i + 1) - inputORIG(i - 1)) * 0.5;
//				dudxm = inputORIG(i) - inputORIG(i - 1);
//
//				//derivORIG(i) = Limiter->compute(dudx);
//				derivSYCL_acc[i] = Limiter->compute(dudxp, dudx0, dudxm);
//			}
//		});
//	});
//}
//
//
//void SingleReconstruction_2nd::computeNormal(Queue& queue, const Data &gdata,
//		NumMatrix<REAL,1> &inputORIG, NumMatrix<REAL,1> &lhs, NumMatrix<REAL,1> &rhs, NumMatrix<REAL,1> &derivORIG, Buffer<REAL,1> &inputSYCL, Buffer<REAL, 1>& derivSYCL)
//{
//	//std::cout << "SingleReconstruction_2nd" << std::endl << std::flush;
//	// Compute limited derivatives
//	getDeriv_limit(queue, gdata, inputORIG, derivORIG, inputSYCL, derivSYCL);
//
//	for (int i = -2; i <= gdata.mx[dir]+1; ++i){
//
//		// Take into account shifted collocation points
//		REAL shift(0.);
//#if (GEOM != CARTESIAN)
//#if (SHIFTED_COLLOCATION == TRUE)
//		if(dir == 0) {
//			shift = sources->shift_Geom_WE(gdata, i);
//		}
//#endif
//#endif
//
//#if (NON_LINEAR_GRID == CRONOS_ON)
//		REAL delx = gdata.getCen_dx(dir, i  );
//		lhs(i) = input(i) - (0.5+shift)*deriv(i)*delx;
//		rhs(i) = lhs(i) + deriv(i)*delx;
//#else
//		lhs(i) = inputORIG(i) - (0.5+shift)*derivORIG(i);
//		rhs(i) = lhs(i) + derivORIG(i);
//#endif
//
//	}
//}

//void SingleReconstruction_2nd::computePar(const Data &gdata,
//		NumMatrix<REAL,1> &inputPar, NumMatrix<REAL,1> &lhs, NumMatrix<REAL,1> &rhs)
//{
//	for (int i = -2; i <= gdata.mx[dir]+1; ++i){
//		lhs(i) = inputPar(i-1);
//		rhs(i) = inputPar(i);
//	}
//}


//void SingleReconstruction_2nd::computePerp(Queue& queue, const Data &gdata,
//		NumMatrix<REAL,1> &inputPerpPORIG, NumMatrix<REAL,1> &inputPerpMORIG,
//		NumMatrix<REAL,1> &lhs, NumMatrix<REAL,1> &rhs, NumMatrix<REAL,1> &derivPORIG, Buffer<REAL, 1>& inputPerpPSYCL, Buffer<REAL, 1>& inputPerpMSYCL, Buffer<REAL, 1>& derivPSYCL)
//{
//	getDeriv_limit(queue, gdata, inputPerpPORIG, derivPORIG, inputPerpPSYCL, derivPSYCL);
//	getDeriv_limit(queue, gdata, inputPerpMORIG, derivM, inputPerpMSYCL, derivPSYCL);
//	for (int i = -2; i <= gdata.mx[dir]+1; ++i){
//
//#if (NON_LINEAR_GRID == CRONOS_ON)
//		REAL delx = gdata.getCen_dx(dir, i  );
//		derivPerp(i) = (derivPORIG(i) + derivM(i))*0.5;
//		lhs(i) = (inputPerpPORIG(i) + inputPerpMORIG(i) - derivPerp(i)*delx)*0.5;
//		rhs(i) =  lhs(i) + derivPerp(i)*delx;
//#else
//		derivPerp(i) = (derivPORIG(i) + derivM(i))*0.5;
//		lhs(i) = (inputPerpPORIG(i) + inputPerpMORIG(i) - derivPerp(i))*0.5;
//		rhs(i) =  lhs(i) + derivPerp(i);
//#endif
//
//	}  
//}




//Reconstruction2D_2nd::Reconstruction2D_2nd(const Data &gdata, const int &dir)  : Reconstruction2D(gdata, dir){
//	int substep = 0;
//	int limiterID = -1;
//
//	if (limiterID == -1 && substep > -1)
//	{
//		if (value_exists("Limiter_" + ToString(substep)))
//		{
//			limiterID = static_cast<int>(value("Limiter_" + ToString(substep)));
//		} else
//		{
//			limiterID = static_cast<int>(value("Limiter"));
//		}
//	}
//
//	cout << "2D reconst" << "\t on substep " << substep << " using limiter: " << limiterID << endl;
//
//	Limiter = new limiter(limiterID);
//}



//void Reconstruction2D_2nd::computeNormal(const Data &gdata,
//                                     NumMatrix<REAL,2> &input,
//                                     NumMatrix<REAL,2> &deriv_dir0,
//                                     NumMatrix<REAL,2> &deriv_dir1,
//                                     NumMatrix<REAL,2> &recLL,
//                                     NumMatrix<REAL,2> &recLR,
//                                     NumMatrix<REAL,2> &recRL,
//                                     NumMatrix<REAL,2> &recRR)
//{
//
//#ifdef CRONOS_SAVEMEM
//	/* Will only compute derivatives if not saved from earlier
//	   Otherwise they are taken from input.
//	*/
//	getDerivs(gdata, input, deriv_dir0, deriv_dir1);
//#endif
//	for (int i = -1; i <= gdata.mx[dir0]+1; ++i){
//		for (int j = -1; j <= gdata.mx[dir1]+1; ++j){
//
//#if (NON_LINEAR_GRID == CRONOS_ON)
//			REAL delx0 = gdata.getCen_dx(dir0, i  );
//			REAL delx1 = gdata.getCen_dx(dir1, j  );
//
//			recLL(i,j) = input(i,j) - 0.5*(deriv_dir0(i,j)*delx0 +
//			                               deriv_dir1(i,j)*delx1);
//			recLR(i,j) = recLL(i,j) + deriv_dir1(i,j)*delx1;
//			recRL(i,j) = recLL(i,j) + deriv_dir0(i,j)*delx0;
//			recRR(i,j) = recRL(i,j) + deriv_dir1(i,j)*delx1;
//#else
//			recLL(i,j) = input(i,j) - 0.5*deriv_dir0(i,j) - 0.5*deriv_dir1(i,j);
//			recLR(i,j) = recLL(i,j) + deriv_dir1(i,j);
//			recRL(i,j) = recLL(i,j) + deriv_dir0(i,j);
//			recRR(i,j) = recRL(i,j) + deriv_dir1(i,j);
//#endif
//
//		}
//	}
//}


//void Reconstruction2D_2nd::computePerp(const Data &gdata,
//                                   NumMatrix<REAL,2> &input,
//                                   NumMatrix<REAL,2> &deriv_inp,
//                                   const int &dir,
//                                   NumMatrix<REAL,2> &recLL,
//                                   NumMatrix<REAL,2> &recLR,
//                                   NumMatrix<REAL,2> &recRL,
//                                   NumMatrix<REAL,2> &recRR)
//{
//	assert(dir == dir0 || dir == dir1);
//
//#ifdef CRONOS_SAVEMEM
//	getDeriv_limit(gdata, input, deriv_inp, dir);
//#endif
//	if(dir == dir0) {
//		for (int i = -1; i <= gdata.mx[dir0]+1; ++i){
//			for (int j = -1; j <= gdata.mx[dir1]+1; ++j){
//
//#if (NON_LINEAR_GRID == CRONOS_ON)
//				REAL delxm = gdata.getCen_dx(dir, i);
//				recLL(i,j) = input(i,j-1) - 0.5*deriv_inp(i,j-1)*delxm;
//				recRL(i,j) = recLL(i,j) + deriv_inp(i,j-1)*delxm;
//				REAL delxp = delxm; // For orthogonal grids
//				recLR(i,j) = input(i,  j) - 0.5*deriv_inp(i,j  )*delxp;
//				recRR(i,j) = recLR(i,j) + deriv_inp(i,j  )*delxp;
//#else
//				recLL(i,j) = input(i,j-1) - 0.5*deriv_inp(i,j-1);
//				recRL(i,j) = recLL(i,j) + deriv_inp(i,j-1);
//				recLR(i,j) = input(i,  j) - 0.5*deriv_inp(i,j  );
//				recRR(i,j) = recLR(i,j) + deriv_inp(i,j  );
//#endif
//
//			}
//		}
//	} else if (dir == dir1) {
//		for (int i = -1; i <= gdata.mx[dir0]+1; ++i){
//			for (int j = -1; j <= gdata.mx[dir1]+1; ++j){
//#if (NON_LINEAR_GRID == CRONOS_ON)
//				REAL delxm = gdata.getCen_dx(dir, j);
//				recLL(i,j) = input(i-1,j) - 0.5*deriv_inp(i-1,j)*delxm;
//				recLR(i,j) = recLL(i,j) + deriv_inp(i-1,j)*delxm;
//				REAL delxp = delxm;
//				recRL(i,j) = input(i  ,j) - 0.5*deriv_inp(i  ,j)*delxp;
//				recRR(i,j) = recRL(i,j) + deriv_inp(i  ,j)*delxp;
//#else
//				recLL(i,j) = input(i-1,j) - 0.5*deriv_inp(i-1,j);
//				recLR(i,j) = recLL(i,j) + deriv_inp(i-1,j);
//				recRL(i,j) = input(i  ,j) - 0.5*deriv_inp(i  ,j);
//				recRR(i,j) = recRL(i,j) + deriv_inp(i  ,j);
//#endif
//			}
//		}
//	} 
//  
//}



//void Reconstruction2D_2nd::getDerivs(const Data &gdata,
//                                 const NumMatrix<REAL,2> &input,
//                                 NumMatrix<REAL,2> &deriv_dir0,
//                                 NumMatrix<REAL,2> &deriv_dir1) {
//
//	for (int i = -2; i <= gdata.mx[dir0]+1; ++i){
//		for (int j = -2; j <= gdata.mx[dir1]+1; ++j){
//
//			// Compute derivatives and store in dudxp, dudx0, dudxm
//			getDeriv(gdata, input, dir0, i, j);
//
//			deriv_dir0(i,j) =  Limiter->compute(dudxp,dudx0,dudxm);
//		}
//	}
//
//
//	for (int i = -2; i <= gdata.mx[dir0]+1; ++i){
//		for (int j = -2; j <= gdata.mx[dir1]+1; ++j){
//
//			// Compute derivatives and store in dudyp, dudy0, dudym
//			getDeriv(gdata, input, dir1, i, j);
//
//			deriv_dir1(i,j) =  Limiter->compute(dudyp,dudy0,dudym);
//		}
//	}
//
//}


//void Reconstruction2D_2nd::getDeriv_limit(const Data &gdata,
//                                const NumMatrix<REAL,2> &input,
//                                NumMatrix<REAL,2> &deriv,
//                                int dir) {
//	if(dir == dir0) {
//		for (int i = -2; i <= gdata.mx[dir0]+1; ++i){
//			for (int j = -2; j <= gdata.mx[dir1]+1; ++j){
//
//				// Compute derivatives and store in dudxp, dudx0, dudxm
//				getDeriv(gdata, input, dir0, i, j);
//
//				deriv(i,j) =  Limiter->compute(dudxp,dudx0,dudxm);
//			}
//		}
//	} else if (dir == dir1) {
//		for (int i = -2; i <= gdata.mx[dir0]+1; ++i){
//			for (int j = -2; j <= gdata.mx[dir1]+1; ++j){
//
//				// Compute derivatives and store in dudyp, dudy0, dudym
//				getDeriv(gdata, input, dir1, i, j);
//
//				deriv(i,j) =  Limiter->compute(dudyp,dudy0,dudym);
//			}
//		}
//	}
//}


