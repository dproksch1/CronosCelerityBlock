#include "reconst_const.H"
#include <iostream>
#include <stdlib.h>

using namespace std;

SingleReconstruction_constant::SingleReconstruction_constant(const Data &gdata, const CronosFluid &fluid, int dir, int qReconst, int substep):
	SingleReconstruction(gdata, fluid, dir, qReconst, substep)
{
	if (gdata.rank == 0) {
		cout << "  Using Reconst: " << "constant" << " - " << fluid.get_fieldName(qReconst) << " - dir " << dir << " - substep " << substep << endl;
	}
}


SingleReconstruction_constant::SingleReconstruction_constant(const Data & gdata, int dir, int substep):
	SingleReconstruction(gdata, dir, substep)
{
	if (gdata.rank == 0) {
		cout << "  Using Reconst: " << "constant" << " - " << "GENERIC" << " - dir " << dir << " - substep " << substep << endl;
	}
}

SingleReconstruction_constant::~SingleReconstruction_constant() {
}


void SingleReconstruction_constant::get_Vals_EW(const Data &gdata, phys_fields_0D &xFieldsW,
		phys_fields_0D &xFieldsE, int ix, int iy, int iz)
{
	xFieldsW.uPri(qReconst) = gdata.om[qReconst](ix,iy,iz);
	xFieldsE.uPri(qReconst) = xFieldsW.uPri(qReconst);
}

void SingleReconstruction_constant::get_Vals_SN(const Data &gdata, phys_fields_0D &xFieldsS,
		phys_fields_0D &xFieldsN, int ix, int iy, int iz)
{
	xFieldsS.uPri(qReconst) = gdata.om[qReconst](ix,iy,iz);
	xFieldsN.uPri(qReconst) = xFieldsS.uPri(qReconst);
}

void SingleReconstruction_constant::get_Vals_BT(const Data &gdata, phys_fields_0D &xFieldsB,
		phys_fields_0D &xFieldsT, int ix, int iy, int iz)
{
	xFieldsB.uPri(qReconst) = gdata.om[qReconst](ix,iy,iz);
	xFieldsT.uPri(qReconst) = xFieldsB.uPri(qReconst);
}


void SingleReconstruction_constant::computeNormal(const Data &gdata,
		NumMatrix<REAL,1> &inputORIG, NumMatrix<REAL,1> &lhs, NumMatrix<REAL,1> &rhs, NumMatrix<REAL,1> &deriv, Buffer<REAL, 1>& inputSYCL)
{
	//std::cout << "SingleReconstruction_constant" << std::endl << std::flush;
	for (int i = -2; i <= gdata.mx[dir]+1; ++i){
		lhs(i) = inputORIG(i);
		rhs(i) = inputORIG(i);

	}
}

void SingleReconstruction_constant::computePar(const Data &gdata,
		NumMatrix<REAL,1> &inputPar, NumMatrix<REAL,1> &lhs, NumMatrix<REAL,1> &rhs)
{
	for (int i = -2; i <= gdata.mx[dir]+1; ++i){
		lhs(i) = inputPar(i-1);
		rhs(i) = inputPar(i);
	}
}


void SingleReconstruction_constant::computePerp(const Data &gdata,
		NumMatrix<REAL,1> &inputPerpPORIG, NumMatrix<REAL,1> &inputPerpMORIG,
		NumMatrix<REAL,1> &lhs, NumMatrix<REAL,1> &rhs, NumMatrix<REAL,1> &derivP, Buffer<REAL, 1>& inputPerpPSYCL, Buffer<REAL, 1>& inputPerpMSYCL)
{
	for (int i = -2; i <= gdata.mx[dir]+1; ++i){
		lhs(i) = (inputPerpPORIG(i) + inputPerpMORIG(i))*0.5;
		rhs(i) =  lhs(i);
	}  
}


Reconstruction2D_constant::Reconstruction2D_constant(const Data &gdata, const int &dir)  : Reconstruction2D(gdata, dir){
}



void Reconstruction2D_constant::computeNormal(const Data &gdata,
                                     NumMatrix<REAL,2> &input,
                                     NumMatrix<REAL,2> &deriv_dir0,
                                     NumMatrix<REAL,2> &deriv_dir1,
                                     NumMatrix<REAL,2> &recLL,
                                     NumMatrix<REAL,2> &recLR,
                                     NumMatrix<REAL,2> &recRL,
                                     NumMatrix<REAL,2> &recRR)
{

	for (int i = -1; i <= gdata.mx[dir0]+1; ++i){
		for (int j = -1; j <= gdata.mx[dir1]+1; ++j){
			recLL(i,j) = input(i,j);
			recLR(i,j) = recLL(i,j);
			recRL(i,j) = recLL(i,j);
			recRR(i,j) = recRL(i,j);
		}
	}
}


void Reconstruction2D_constant::computePerp(const Data &gdata,
                                   NumMatrix<REAL,2> &input,
                                   NumMatrix<REAL,2> &deriv_inp,
                                   const int &dir,
                                   NumMatrix<REAL,2> &recLL,
                                   NumMatrix<REAL,2> &recLR,
                                   NumMatrix<REAL,2> &recRL,
                                   NumMatrix<REAL,2> &recRR)
{
	assert(dir == dir0 || dir == dir1);

	if(dir == dir0) {
		for (int i = -1; i <= gdata.mx[dir0]+1; ++i){
			for (int j = -1; j <= gdata.mx[dir1]+1; ++j){
				recLL(i,j) = input(i,j-1);
				recRL(i,j) = recLL(i,j);
				recLR(i,j) = input(i,  j);
				recRR(i,j) = recLR(i,j);
			}
		}
	} else if (dir == dir1) {
		for (int i = -1; i <= gdata.mx[dir0]+1; ++i){
			for (int j = -1; j <= gdata.mx[dir1]+1; ++j){
				recLL(i,j) = input(i-1,j);
				recLR(i,j) = recLL(i,j);
				recRL(i,j) = input(i  ,j) ;
				recRR(i,j) = recRL(i,j);
			}
		}
	} 
  
}



