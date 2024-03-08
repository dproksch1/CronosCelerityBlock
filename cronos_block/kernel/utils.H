#ifndef CRONOS_UTILS_H
#define CRONOS_UTILS_H = 1

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
		TypeFluxPhys,

		TypeMax,
	};

	enum {

		fluidConst_q_rho = 0,
		fluidConst_q_sx,
		fluidConst_q_sy,
		fluidConst_q_sz,
		fluidConst_q_Eges,
		fluidConst_q_Eadd,
		fluidConst_q_Bx,
		fluidConst_q_By,
		fluidConst_q_Bz,

		fluidConstMax,
	};

	enum {

		NumV_ch_p = 0,
		NumV_ch_m,

		NumV_Max,
	};

}

typedef double nomMat[N_OMINT];
struct nom_t {
	nomMat mat;
};

typedef double uPriMat[gpu::FaceMax][N_OMINT];
struct uPri_t {
	uPriMat mat;
};


#endif