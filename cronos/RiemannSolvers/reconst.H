
#ifndef CRONOS_RECONST_H
#define CRONOS_RECONST_H = 1

#include "matrix.H"
#include "data.H"
#include "limiter.H"
#include "fields_local.H"
#include <vector>

namespace gpu {

	enum {

		FaceWest = 0,
		FaceEast,
		FaceSouth,
		FaceNorth,
		FaceBottom,
		FaceTop,

		FaceMax,

	};

	enum {

		TypeUPri = 0,
		TypeUCon,

		TypeMax,

	};

	//class Reconstruction {

		struct Dud {
			REAL _m;
			REAL _0;
			REAL _p;

		};

	//public:

		inline REAL minmodCompute(REAL Theta, REAL& up, REAL& uc, REAL& um) {
			up *= Theta;
			um *= Theta;
			if (up * um > 0.) {
				if (up > 0) {
					return fmin(up, (fmin(uc, um)));
				} else if (up < 0) {
					return fmax(up, (fmax(uc, um)));
				} else {
					return 0.0;
					//cerr << " Error " << endl;
					//exit(1);
				}
			} else {
				return 0.;
			}
		}

		template<typename AccessorType>
		auto getDerivs(AccessorType& omAcc, size_t ix, size_t iy, size_t iz) {
			//size_t q = qReconst;
			Dud x;
			Dud y;
			Dud z;

			x._m = omAcc[ix][iy][iz] - omAcc[ix - 1][iy][iz];
			//	REAL dudx0 = (omAcc(ix+1,iy,iz) - omAcc(ix-1,iy,iz))*0.5;
			x._p = omAcc[ix + 1][iy][iz] - omAcc[ix][iy][iz];
			x._0 = 0.5 * (x._m + x._p);

			y._m = omAcc[ix][iy][iz] - omAcc[ix][iy - 1][iz];
			//	REAL dudy0 = (omAcc(ix,iy+1,iz) - omAcc(ix,iy-1,iz))*0.5;
			y._p = omAcc[ix][iy + 1][iz] - omAcc[ix][iy][iz];
			y._0 = 0.5 * (y._m + y._p);

			z._m = omAcc[ix][iy][iz] - omAcc[ix][iy][iz - 1];
			//	REAL dudz0 = (omAcc(ix,iy,iz+1) - omAcc(ix,iy,iz-1))*0.5;
			z._p = omAcc[ix][iy][iz + 1] - omAcc[ix][iy][iz];
			z._0 = 0.5 * (z._m + z._p);

			return std::tuple{x, y, z};
		}

		template<typename AccessorType>
		auto prepareDerivs(AccessorType& omAcc, size_t ix, size_t iy, size_t iz) {

			//! Compute derivate from at given position
			auto [dudx, dudy, dudz] = getDerivs(omAcc, ix, iy, iz);

			// TODO: read from project parameters
			REAL Theta = 1;

			REAL deriv_x = minmodCompute(Theta, dudx._p, dudx._0, dudx._m);
			REAL deriv_y = minmodCompute(Theta, dudy._p, dudy._0, dudy._m);
			REAL deriv_z = minmodCompute(Theta, dudz._p, dudz._0, dudz._m);

			return std::tuple{ deriv_x, deriv_y, deriv_z };
		}


		template<typename AccessorType, typename FieldAccessorType>
		void get_Vals_EW(AccessorType& omAcc, FieldAccessorType& fieldAcc, size_t ix, size_t iy, size_t iz, REAL deriv_x) {
			REAL shift(0.);
			fieldAcc[FaceWest][TypeUPri] = omAcc[ix][iy][iz] - (0.5 + shift) * deriv_x;
			//xFieldsW.uPri(q) = omAcc[ix][iy][iz] - (0.5 + shift) * deriv_x;
			fieldAcc[FaceEast][TypeUPri] = fieldAcc[FaceWest][TypeUPri] + deriv_x;
			//xFieldsE.uPri(q) = xFieldsW.uPri(q) + deriv_x;
		}

		template<typename AccessorType, typename FieldAccessorType>
		void get_Vals_SN(AccessorType& omAcc, FieldAccessorType& fieldAcc, size_t ix, size_t iy, size_t iz, REAL deriv_y) {
			REAL shift(0.);
			fieldAcc[FaceSouth][TypeUPri] = omAcc[ix][iy][iz] - 0.5 * deriv_y;
			//xFieldsS.uPri(q) = gdata.om[q](ix, iy, iz) - 0.5 * deriv_y;
			fieldAcc[FaceNorth][TypeUPri] = fieldAcc[FaceSouth][TypeUPri] + deriv_y;
			//xFieldsN.uPri(q) = xFieldsS.uPri(q) + deriv_y;

		}

		template<typename AccessorType, typename FieldAccessorType>
		void get_Vals_BT(AccessorType& omAcc, FieldAccessorType& fieldAcc, size_t ix, size_t iy, size_t iz, REAL deriv_z) {
			REAL shift(0.);
			fieldAcc[FaceBottom][TypeUPri] = omAcc[ix][iy][iz] - 0.5 * deriv_z;
			//xFieldsB.uPri(q) = gdata.om[q](ix, iy, iz) - 0.5 * deriv_z;
			fieldAcc[FaceTop][TypeUPri] = fieldAcc[FaceBottom][TypeUPri] + deriv_z;
			//xFieldsT.uPri(q) = xFieldsB.uPri(q) + deriv_z;

		}

		//std::vector<SingleReconstruction> ListReconstructionNormal;
		//std::vector<size_t> ListNormal;

		template<typename AccessorType, typename FieldAccessorType>
		void compute(AccessorType& omAcc, FieldAccessorType& fieldAcc, size_t ix, size_t iy, size_t iz) {
			// Todo: Andere Generalisierung für punktweise rekonstruktion einführen?
			//for (size_t q = 0; q < ListNormal.size(); ++q) {
				// First compute all derivatives
				//ListReconstructionNormal[q].prepareDerivs(omAcc, ix, iy, iz);
				auto [deriv_x, deriv_y, deriv_z] = prepareDerivs(omAcc, ix, iy, iz);

				// Get Values in east and west direction:
				//ListReconstructionNormal[q].get_Vals_EW(omAcc, *allFields[0], *allFields[1], ix, iy, iz);
				get_Vals_EW(omAcc, fieldAcc, ix, iy, iz, deriv_x);

				// Get Values in South and North directions:
				//ListReconstructionNormal[q].get_Vals_SN(omAcc, *allFields[2], *allFields[3], ix, iy, iz);
				get_Vals_SN(omAcc, fieldAcc, ix, iy, iz, deriv_y);

				// Get Values in Bottom and Top directions:
				//ListReconstructionNormal[q].get_Vals_BT(omAcc, *allFields[4], *allFields[5], ix, iy, iz);
				get_Vals_BT(omAcc, fieldAcc, ix, iy, iz, deriv_z);
			//}
		}

	//};

}

class SingleReconstruction{
public:

	enum DudIndex {
		_p,
		_0, 
		_m,
		Max_Index,
	};

	enum DudDir {
		_x,
		_y,
		_z,
		Max_Dir,
	};

	SingleReconstruction() {}
	SingleReconstruction(const Data &, int dir, int substep);
	SingleReconstruction(const Data &, const CronosFluid &fluid, int dir, int qReconst, int substep);
	virtual ~SingleReconstruction();

	/**
	 * Compute reconstruction of cell-centered variable
	 * */
	virtual void computeNormal(Queue& /*queue*/, const Data& /*gdata*/, NumMatrix<REAL, 1>& /*inputORIG*/, NumMatrix<REAL, 1>& /*lhs*/,
		NumMatrix<REAL, 1>& /*rhs*/, NumMatrix<REAL, 1>& /*derivORIG*/, Buffer<REAL, 1> /*lhsSYCL*/, Buffer<REAL, 1> /*rhsSYCL*/, Buffer<REAL, 1> /*inputSYCL*/, Buffer<REAL, 1> /*derivSYCL*/) = 0;

	/**
	 * Compute reconstruction of face-centered variable in direction of respective cell face
	 * */
	virtual void computePar(const Data &, NumMatrix<REAL,1> &, NumMatrix<REAL,1> &,
	                NumMatrix<REAL,1> &) = 0;

	/**
	 * Compute reconstruction of face-centered variable in perpendicular direction
	 * */
	virtual void computePerp(Queue& /*queue*/, const Data& /*gdata*/, NumMatrix<REAL, 1>& /*inputPerpPORIG*/, NumMatrix<REAL, 1>& /*inputPerpMORIG*/,
							 NumMatrix<REAL, 1>& /*lhs*/, NumMatrix<REAL, 1>& /*rhs*/, NumMatrix<REAL, 1>& /*derivP*/,
							 Buffer<REAL, 1> /*inputPerpPSYCL*/, Buffer<REAL, 1> /*inputPerpMSYCL*/, Buffer<REAL, 1> /*derivSYCL*/) = 0;

	virtual void prepareDerivs(const Data& gdata, int ix, int iy, int iz) {};

	virtual void get_Vals_EW(const Data &, phys_fields_0D &xFieldsW,
			phys_fields_0D &xFieldsE, int, int, int) = 0;
	virtual void get_Vals_SN(const Data &, phys_fields_0D &xFieldsS,
			phys_fields_0D &xFieldsN, int, int, int) = 0;
	virtual void get_Vals_BT(const Data &, phys_fields_0D &xFieldsB,
			phys_fields_0D &xFieldsT, int, int, int) = 0;
protected:
	/**
	 * Compute left-handed, centred, and right-handed derivatives in given direction
	 * */
	void getDeriv(Queue& /*queue*/, const Data& /*gdata*/, const NumMatrix<REAL, 1>& /*inputORIG*/, int /*iPos*/, NumMatrix<REAL, 1>& /*dudxORIG*/, Buffer<REAL, 1> /*inputSYCL*/, Buffer<REAL, 1> /*dudxSYCL*/);

	/**
	 * Compute left-handed, centred, and right-handed derivatives in all spatial directions
	 * --> block-structured version
	 * */
	void getDerivs(const Data &, int ix, int iy, int iz);
	int dir;
	int qReconst;
	int substep;

	NumMatrix<REAL, 1> derivPerp, derivM;
	REAL deriv_x, deriv_y, deriv_z;
	REAL dudxp_q, dudx0_q, dudxm_q;
	REAL dudyp_q, dudy0_q, dudym_q;
	REAL dudzp_q, dudz0_q, dudzm_q;
	//std::vector<Buffer<REAL, 1>> dud_q = std::vector<Buffer<REAL, 1>>(DudDir::Max_Dir, Range<1>(DudIndex::Max_Index));

	//REAL dudxp, dudx0, dudxm;
	//Buffer<REAL, 1> dudx;
};

class Reconstruction{
public:
	Reconstruction(){}
	Reconstruction(const Data &, int dir, const CronosFluid &fluid, int substep=-1);
	Reconstruction(const Data &, int dir, int num, int substep=-1);
	virtual ~Reconstruction();

	/**
	 * Compute reconstruction in case of 1D arrays
	 * */
	void compute(Queue& /*queue*/, const Data &, fields_1D &, phys_fields_1D &, phys_fields_1D &);

	/**
	 * Compute reconstruction for block-structured code
	 * */
	void compute(const Data &, std::vector<phys_fields_0D> &allFields, int, int, int, Direction dir = (Direction)-1);

	/**
	 * Compute reconstruction in case of 1D arrays
	 * alternative function parameters
	 * */
	void compute(const Data &, NumMatrix<REAL,1> [], NumMatrix<REAL,1> [],
	             NumMatrix<REAL,1> [], NumMatrix<REAL,1> []);

protected:
	void set_singleReconstructions(const Data &);
	void set_singleReconstructions(const Data &, const CronosFluid &fluid);

	int dir;
	int substep;
	std::vector<int>::iterator iter;
	std::vector<int> ListNormal;
	std::vector<int> ListParallel;
	std::vector<int> ListPerp;

	std::vector<SingleReconstruction*> ListReconstructionNormal;
	std::vector<SingleReconstruction*> ListReconstructionPar;
	std::vector<SingleReconstruction*> ListReconstructionPerp;

};




class Reconstruction2D{
public:
	Reconstruction2D(const Data &, const int & dir);
	void compute(const Data &, fields_2D &,
	             phys_fields_2D &, phys_fields_2D &,
	             phys_fields_2D &, phys_fields_2D &);
	virtual void computeNormal(const Data &,
	                   NumMatrix<REAL,2> &, NumMatrix<REAL,2> &,
	                   NumMatrix<REAL,2> &, NumMatrix<REAL,2> &,
	                   NumMatrix<REAL,2> &, NumMatrix<REAL,2> &,
	                   NumMatrix<REAL,2> &) = 0;
	virtual void computePerp(const Data &, NumMatrix<REAL,2> &, NumMatrix<REAL,2> &,
	                 const int &,
	                 NumMatrix<REAL,2> &, NumMatrix<REAL,2> &,
	                 NumMatrix<REAL,2> &, NumMatrix<REAL,2> &) = 0;
protected:
	void getDeriv(const Data &, const NumMatrix<REAL,2> &,
	              int dir, int i, int j);
	REAL dudxp, dudx0, dudxm;
	REAL dudyp, dudy0, dudym;
	int dir0, dir1;
};
#endif
