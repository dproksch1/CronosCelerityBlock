#ifndef CRONOS_FIELDS_1D
#define CRONOS_FIELDS_1D = 1

#include "data.H"
#include "matrix.H"
#include "constants.H"
#include "queue.H"
#include <vector>

using OmLocType = std::vector<Buffer<REAL, 1>>;

class fields_1D {
public:

//	fields_1D(const Data &, int);
	fields_1D(const Data &, int, const CronosFluid &fluid);
//	fields_1D(const Data &, int );
	//fields_1D() {}
	fields_1D(const Data &, int , int, int iFluid=0);
	~fields_1D();
	bool isgeneric();
	int get_num();
	//	NumMatrix<REAL, 1> omLoc[N_OMINT];
	NumMatrix<REAL, 1> *omLocORIG;
	OmLocType omLocSYCL;
	NumMatrix<REAL, 1>* derivORIG;
	OmLocType derivSYCL;
//	NumMatrix<REAL, 1> omLocP[N_OMINT];
//	NumMatrix<REAL, 1> omLocM[N_OMINT];
	NumMatrix<REAL, 1> *omLocPORIG;
	NumMatrix<REAL, 1> *omLocMORIG;
	OmLocType omLocPSYCL;
	OmLocType omLocMSYCL;
	
//#if(FLUID_TYPE == CRONOS_MULTIFLUID)
//	NumMatrix<REAL, 1> *ptotal;
//	NumMatrix<REAL, 1> *v_ch_p;
//	NumMatrix<REAL, 1> *v_ch_m;
//#else
	NumMatrix<REAL, 1> ptotalORIG;
	Buffer<REAL, 1> ptotalSYCL;
	NumMatrix<REAL, 1> v_ch_pORIG;
	NumMatrix<REAL, 1> v_ch_mORIG;
	Buffer<REAL, 1> v_ch_pSYCL;
	Buffer<REAL, 1> v_ch_mSYCL;

	NumMatrix<int, 1> carbuncle_flag;
//#endif

	//	NumMatrix<REAL, 1> flux[N_OMINT];
	NumMatrix<REAL, 1>* fluxORIG;
	OmLocType fluxSYCL;
private:
	//	void SetupFields(const Data &);
	int num, dir;
	bool generic;
	bool has_perp;
};


class fields_2D {
public:
	fields_2D(const Data &, const int &);
	~fields_2D();
	NumMatrix<REAL, 2> v_x, v_y, v_z, B_x, B_y, B_z;
	NumMatrix<REAL, 2> dBxdx, dBydx, dBzdx;
	NumMatrix<REAL, 2> dvxdx[2], dvydx[2], dvzdx[2];
	
	NumMatrix<REAL, 2> v_ch_p[2];
	NumMatrix<REAL, 2> v_ch_m[2];
	
	NumMatrix<REAL, 2> v_cor_p[2];
	NumMatrix<REAL, 2> v_cor_m[2];
	
	NumMatrix<REAL, 2> emf;

	NumMatrix<REAL, 2> *fluxWE;
	NumMatrix<REAL, 2> *fluxSN;
	NumMatrix<REAL, 2> *fluxBT;
//	NumMatrix<REAL, 2> fluxWE[7];
//	NumMatrix<REAL, 2> fluxSN[7];
//	NumMatrix<REAL, 2> fluxBT[7];
	NumMatrix<REAL, 2> Ex, Ey, Ez;

	// For special form of the Stone CT
	NumMatrix<REAL, 2> dExdyjm1_4, dExdyjm3_4, dExdzkm1_4, dExdzkm3_4;
	NumMatrix<REAL, 2> dEydxim1_4, dEydxim3_4, dEydzkm1_4, dEydzkm3_4;
	NumMatrix<REAL, 2> dEzdxim1_4, dEzdxim3_4, dEzdyjm1_4, dEzdyjm3_4;
private:
	int dir0, dir1;
};


class phys_fields_2D {
public:
	phys_fields_2D(const Data &, const int &);
	NumMatrix<REAL, 2> v_x, v_y, v_z, B_x, B_y, B_z;
	NumMatrix<REAL, 2> emf;
private:
	int dir0, dir1;
};


//class fields_0D {
//public:
//	fields_0D(const Data &){}
//	REAL omLoc[N_OMINT];
//	REAL omLocP[N_OMINT];
//	REAL omLocM[N_OMINT];
//	REAL deriv[N_OMINT];
//	REAL derivP[N_OMINT];
//	REAL derivM[N_OMINT];
//
//	REAL v_ch_p[N_OMINT];
//	REAL v_ch_m[N_OMINT];
//
//	REAL flux[N_OMINT];
//};

class phys_fields_1D {
public:
	phys_fields_1D(const Data &, const int &);
	phys_fields_1D(const Data &, int, int);
	//phys_fields_1D() {}
	~phys_fields_1D();
	NumMatrix<REAL, 1> pthermORIG;
	Buffer<REAL, 1> pthermSYCL;
	NumMatrix<REAL, 1> ptotalORIG;
	Buffer<REAL, 1> ptotalSYCL;
	//	NumMatrix<REAL, 1> uCon[N_OMINT];
	//	NumMatrix<REAL, 1> *uLoc;
	NumMatrix<REAL, 1> *uPriORIG;
	OmLocType uPriSYCL;
	NumMatrix<REAL, 1>* uConORIG;
	OmLocType uConSYCL;
	NumMatrix<REAL, 1> *fluxORIG;
	OmLocType fluxSYCL;
#if (USE_COROTATION == CRONOS_ON)
	NumMatrix<REAL, 1> uInertial[3];
#endif

};


class phys_fields_0D {
public:
	phys_fields_0D(const Data &, int, const CronosFluid &fluid);
	phys_fields_0D() {}
	~phys_fields_0D() {};
	bool isgeneric();
	int get_num() const;

	NumArray<REAL> uPri;
	NumArray<REAL> uCon;
	NumArray<REAL> flux_phys;
#if (USE_COROTATION == CRONOS_ON)
	NumArray<REAL> uInertial;
#endif

	double ptotal, ptherm;
	int carbuncle_flag;
private:
	int face, num;
};


class num_fields_0D {
public:
	num_fields_0D(const Data &, const CronosFluid &fluid);
	num_fields_0D() {};
	~num_fields_0D() {};
	REAL v_ch_p;
	REAL v_ch_m;
	NumArray<REAL> flux_num;
	double ptotal_num;
	int carbuncle_flag;
private:
	int num;
};

#endif
