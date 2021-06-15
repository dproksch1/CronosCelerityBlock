#include "fields_local.H"
#include "fieldLists.H"
#include "queue.H"
#include <iostream>

using namespace std;

fields_1D::fields_1D(const Data &gdata, int dir, const CronosFluid &fluid)
	: ptotalSYCL (Buffer<REAL, 1>(Range<1>(gdata.mx[dir] + 4))),
	  v_ch_pSYCL (Buffer<REAL, 1>(Range<1>(gdata.mx[dir] + 4))),
	  v_ch_mSYCL (Buffer<REAL, 1>(Range<1>(gdata.mx[dir] + 4)))
{


	this->dir = dir;
	generic = true; // Code native stuff

#if (FLUID_TYPE == CRONOS_MULTIFLUID)
	fieldLists List(fluid, dir);
	this->num = fluid.get_N_OMINT();
//	fieldLists List(gdata.fluids, dir);
//	this->num = gdata.fluids->get_N_OMINT_ALL();
//	fieldLists List(gdata.fluids->fluids, gdata.fluids->get_fluidType(iFluid), dir);
#else
	fieldLists List(fluid, dir);
//	fieldLists List(gdata.fluid, dir);
	this->num = N_OMINT;
#endif

	omLocORIG = new NumMatrix<REAL,1> [num];
	auto gdataSize = gdata.mx[dir] + 2 + 3 + 1;
	omLocSYCL = std::vector(num, Buffer<REAL, 1>(Range<1>(gdataSize)));
	for (auto& e : omLocSYCL) {
		clearBufferOnHost(e);
	}

	// , Buffer<REAL, 1>(Range<1>(gdataSize)));
	derivORIG = new NumMatrix<REAL, 1>[num];
	derivSYCL = std::vector(num, Buffer<REAL, 1>(Range<1>(gdataSize)));
	for (auto& e : derivSYCL) {
		clearBufferOnHost(e);
		//auto e_acc = e.get_access<cl::sycl::access::mode::discard_write>();
		//for (int i = 0; i < gdataSize; ++i) {
		//	e_acc[i] = 0.0;
		//}
	}
	fluxORIG = new NumMatrix<REAL,1> [num];
	omLocPORIG = new NumMatrix<REAL,1> [num];
	omLocMORIG = new NumMatrix<REAL,1> [num];
	omLocPSYCL = std::vector(num, Buffer<REAL, 1>(Range<1>(gdataSize)));
	for (auto& e : omLocPSYCL) {
		clearBufferOnHost(e);
		//auto e_acc = e.get_access<cl::sycl::access::mode::discard_write>();
		//for (int i = 0; i < gdataSize; ++i) {
		//	e_acc[i] = 0.0;
		//}
	}
	omLocMSYCL = std::vector(num, Buffer<REAL, 1>(Range<1>(gdataSize)));
	for (auto& e : omLocMSYCL) {
		clearBufferOnHost(e);
		//auto e_acc = e.get_access<cl::sycl::access::mode::discard_write>();
		//for (int i = 0; i < gdataSize; ++i) {
		//	e_acc[i] = 0.0;
		//}
	}
	has_perp = true;

	for(int q=0; q<num; ++q) {
		// Identify fields that require different reconstruction
		if(List.isGenTypePerpendicular(q)) {
			omLocPORIG[q].resize(Index::set(-3),Index::set(gdata.mx[dir]+2));
			omLocPORIG[q].clear();
			omLocMORIG[q].resize(Index::set(-3),Index::set(gdata.mx[dir]+2));
			omLocMORIG[q].clear();
		} else {
			omLocORIG[q].resize(Index::set(-3),Index::set(gdata.mx[dir]+2));
			omLocORIG[q].clear();
		}
//		cout << " Loop: " << q << " " << num << " " << dir << " " << omLoc[q].getHigh(0) << endl << endl;
	}

	for(int q=0; q<num; ++q) {
		derivORIG[q].resize(Index::set(-3),Index::set(gdata.mx[dir]+2));
		derivORIG[q].clear();
	}

//#if(FLUID_TYPE == CRONOS_MULTIFLUID)
//	// In multifluid case make one vChar and one ptotal for each fluid
//	int numFluids = gdata.fluids->get_numFluids();
//
//	ptotal = new NumMatrix<REAL,1> [numFluids];
//	for(int iFluid; iFluid<numFluids; ++iFluid) {
//		ptotal[iFluid].resize(Index::set(-2),Index::set(gdata.mx[dir]+1));
//	}
//
//	v_ch_p = new NumMatrix<REAL,1> [numFluids];
//	v_ch_m = new NumMatrix<REAL,1> [numFluids];
//	for(int iFluid; iFluid<numFluids; ++iFluid) {
//		v_ch_p[iFluid].resize(Index::set(-2),Index::set(gdata.mx[dir]+1));
//		v_ch_m[iFluid].resize(Index::set(-2),Index::set(gdata.mx[dir]+1));
//	}
//#else
	ptotalORIG.resize(Index::set(-2),Index::set(gdata.mx[dir]+1));
	ptotalORIG.clear();

	int size = gdata.mx[dir] + 1 + 2 + 1;

	clearBufferOnHost(ptotalSYCL);

	v_ch_pORIG.resize(Index::set(-2),Index::set(gdata.mx[dir]+1));
	v_ch_pORIG.clear();

	clearBufferOnHost(v_ch_pSYCL);
	//auto v_ch_p_acc = v_ch_pSYCL.get_access<cl::sycl::access::mode::discard_write>();
	//for (int i = 0; i < size; ++i) {
	//	v_ch_p_acc[i] = 0.0;
	//}

	v_ch_mORIG.resize(Index::set(-2),Index::set(gdata.mx[dir]+1));
	v_ch_mORIG.clear();

	clearBufferOnHost(v_ch_mSYCL);
	//auto v_ch_m_acc = v_ch_mSYCL.get_access<cl::sycl::access::mode::discard_write>();
	//for (int i = 0; i < size; ++i) {
	//	v_ch_m_acc[i] = 0.0;
	//}

//#endif

	if(gdata.use_carbuncleFlag) {
		carbuncle_flag.resize(Index::set(-2),Index::set(gdata.mx[dir]+2));
		carbuncle_flag.clear();
	}

	for(int q=0; q<num; ++q) {
		fluxORIG[q].resize(Index::set(-2),Index::set(gdata.mx[dir]+1));
		fluxORIG[q].clear();
	}

	fluxSYCL = std::vector(num, Buffer<REAL, 1>(Range<1>(size)));
	for (auto& e : fluxSYCL) {
		clearBufferOnHost(e);
		//auto e_acc = e.get_access<cl::sycl::access::mode::discard_write>();
		//for (int i = 0; i < size; ++i) {
		//	e_acc[i] = 0.0;
		//}
	}
}


fields_1D::fields_1D(const Data &gdata, int dir, int num, int iFluid)
	: ptotalSYCL (Buffer<REAL, 1>(Range<1>(gdata.mx[dir] + 4))),
	  v_ch_pSYCL (Buffer<REAL, 1>(Range<1>(gdata.mx[dir] + 4))),
	  v_ch_mSYCL (Buffer<REAL, 1>(Range<1>(gdata.mx[dir] + 4)))
{

	this->num = num;
	this->dir = dir;
	generic = false; // User dependent stuff

	omLocORIG = new NumMatrix<REAL,1> [num];
	auto gdataSize = gdata.mx[dir] + 2 + 3;
	omLocSYCL = std::vector(num, Buffer<REAL, 1>(Range<1>(gdataSize)));
	derivORIG = new NumMatrix<REAL,1> [num];
	fluxORIG = new NumMatrix<REAL,1> [num];
	has_perp = false;

#if (FLUID_TYPE == CRONOS_MULTIFLUID)
	fieldLists List(gdata.fluids->fluids[iFluid], dir);
	//	fieldLists List(gdata.fluids->fluids, gdata.fluids->get_fluidType(iFluid), dir);
#else
	fieldLists List(gdata.fluid, dir);
#endif

	for(int q=0; q<num; ++q) {
		omLocORIG[q].resize(Index::set(-3),Index::set(gdata.mx[dir]+2));
		omLocORIG[q].clear();

		derivORIG[q].resize(Index::set(-3),Index::set(gdata.mx[dir]+2));
		derivORIG[q].clear();

		fluxORIG[q].resize(Index::set(-2),Index::set(gdata.mx[dir]+1));
		fluxORIG[q].clear();
	}

//#if(FLUID_TYPE == CRONOS_MULTIFLUID)
//	// In multifluid case make one vChar and one ptotal for each fluid
//	int numFluids = gdata.fluids->get_numFluids();
//
//	ptotal = new NumMatrix<REAL,1> [numFluids];
//	for(int iFluid; iFluid<numFluids; ++iFluid) {
//		ptotal[iFluid].resize(Index::set(-2),Index::set(gdata.mx[dir]+1));
//	}
//
//	v_ch_p = new NumMatrix<REAL,1> [numFluids];
//	v_ch_m = new NumMatrix<REAL,1> [numFluids];
//	for(int iFluid; iFluid<numFluids; ++iFluid) {
//		v_ch_p[iFluid].resize(Index::set(-2),Index::set(gdata.mx[dir]+1));
//		v_ch_m[iFluid].resize(Index::set(-2),Index::set(gdata.mx[dir]+1));
//	}
//#else
	v_ch_pORIG.resize(Index::set(-2),Index::set(gdata.mx[dir]+1));
	v_ch_pORIG.clear();

	//const int size = gdata.mx[dir] + 1 + 2 + 1;
	clearBufferOnHost(v_ch_pSYCL);
	//auto v_ch_p_acc = v_ch_pSYCL.get_access<cl::sycl::access::mode::discard_write>();
	//for (int i = 0; i < size; ++i) {
	//	v_ch_p_acc[i] = 0.0;
	//}

	v_ch_mORIG.resize(Index::set(-2),Index::set(gdata.mx[dir]+1));
	v_ch_mORIG.clear();

	clearBufferOnHost(v_ch_mSYCL);
	//auto v_ch_m_acc = v_ch_mSYCL.get_access<cl::sycl::access::mode::discard_write>();
	//for (int i = 0; i < size; ++i) {
	//	v_ch_m_acc[i] = 0.0;
	//}

	ptotalORIG.resize(Index::set(-2),Index::set(gdata.mx[dir]+1));
	ptotalORIG.clear();
	clearBufferOnHost(ptotalSYCL);
//#endif
}


fields_1D::~fields_1D() {
	delete [] omLocORIG;
	delete [] derivORIG;
	delete [] fluxORIG;
	if(has_perp) {
		delete [] omLocPORIG;
		delete [] omLocMORIG;
	}
//#if(FLUID_TYPE == CRONOS_MULTIFLUID)
//	delete [] v_ch_m;
//	delete [] v_ch_p;
//	delete [] ptotal;
//#endif
}


bool fields_1D::isgeneric() {
	return generic;
}

int fields_1D::get_num() {
	return num;
}


fields_2D::fields_2D(const Data &gdata, const int &dir) {

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

#if(FLUID_TYPE == CRONOS_MULTIFLUID)
	int iMagFluid = gdata.fluids->get_i_magFluid();
	int num = gdata.fluids->get_N_OMINT(iMagFluid);
	int q_Bx = gdata.fluids->get_q_Bx();
	int q_By = gdata.fluids->get_q_By();
	int q_Bz = gdata.fluids->get_q_Bz();
#else
	int num = N_OMINT;
	int q_Bx = gdata.fluid.get_q_Bx();
	int q_By = gdata.fluid.get_q_By();
	int q_Bz = gdata.fluid.get_q_Bz();
#endif

	fluxWE = new NumMatrix<REAL,2> [num];
	fluxSN = new NumMatrix<REAL,2> [num];
	fluxBT = new NumMatrix<REAL,2> [num];

#if (CT_TYPE == STONE)
	if(dir == 0) {
		Ex.resize(Index::set(-1, -1),
		          Index::set(gdata.mx[dir0]+1, gdata.mx[dir1]+1));
		fluxSN[q_Bz].resize(Index::set(-2, -2),
		                 Index::set(gdata.mx[dir0]+1, gdata.mx[dir1]+1));
		fluxBT[q_By].resize(Index::set(-2, -2),
		                 Index::set(gdata.mx[dir0]+1, gdata.mx[dir1]+1));

#if (STONE_TYPE == STONE_CENTRE)
		fluxSN[0].resize(Index::set(-2, -2),
		                 Index::set(gdata.mx[dir0]+1, gdata.mx[dir1]+1));
		fluxBT[0].resize(Index::set(-2, -2),
		                 Index::set(gdata.mx[dir0]+1, gdata.mx[dir1]+1));

		dExdyjm1_4.resize(Index::set(-2, -2),
		                  Index::set(gdata.mx[dir0]+1, gdata.mx[dir1]+1));
		dExdyjm3_4.resize(Index::set(-2, -2),
		                  Index::set(gdata.mx[dir0]+1, gdata.mx[dir1]+1));
		dExdzkm1_4.resize(Index::set(-2, -2),
		                  Index::set(gdata.mx[dir0]+1, gdata.mx[dir1]+1));
		dExdzkm3_4.resize(Index::set(-2, -2),
		                  Index::set(gdata.mx[dir0]+1, gdata.mx[dir1]+1));
#endif

	} else if(dir == 1) {
		Ey.resize(Index::set(-1, -1),
		          Index::set(gdata.mx[dir0]+1, gdata.mx[dir1]+1));
		fluxWE[q_Bz].resize(Index::set(-2, -2),
		                 Index::set(gdata.mx[dir0]+1, gdata.mx[dir1]+1));
		fluxBT[q_Bx].resize(Index::set(-2, -2),
		                 Index::set(gdata.mx[dir0]+1, gdata.mx[dir1]+1));

#if (STONE_TYPE == STONE_CENTRE)
		fluxWE[0].resize(Index::set(-2, -2),
		                 Index::set(gdata.mx[dir0]+1, gdata.mx[dir1]+1));
		fluxBT[0].resize(Index::set(-2, -2),
		                 Index::set(gdata.mx[dir0]+1, gdata.mx[dir1]+1));

		dEydxim1_4.resize(Index::set(-2, -2),
		                  Index::set(gdata.mx[dir0]+1, gdata.mx[dir1]+1));
		dEydxim3_4.resize(Index::set(-2, -2),
		                  Index::set(gdata.mx[dir0]+1, gdata.mx[dir1]+1));
		dEydzkm1_4.resize(Index::set(-2, -2),
		                  Index::set(gdata.mx[dir0]+1, gdata.mx[dir1]+1));
		dEydzkm3_4.resize(Index::set(-2, -2),
		                  Index::set(gdata.mx[dir0]+1, gdata.mx[dir1]+1));
#endif
	} else {
		Ez.resize(Index::set(-1, -1),
		          Index::set(gdata.mx[dir0]+1, gdata.mx[dir1]+1));
		fluxWE[q_By].resize(Index::set(-2, -2),
		                 Index::set(gdata.mx[dir0]+1, gdata.mx[dir1]+1));
		fluxSN[q_Bx].resize(Index::set(-2, -2),
		                 Index::set(gdata.mx[dir0]+1, gdata.mx[dir1]+1));

#if (STONE_TYPE == STONE_CENTRE)
		fluxWE[0].resize(Index::set(-2, -2),
		                 Index::set(gdata.mx[dir0]+1, gdata.mx[dir1]+1));
		fluxSN[0].resize(Index::set(-2, -2),
		                 Index::set(gdata.mx[dir0]+1, gdata.mx[dir1]+1));

		dEzdxim1_4.resize(Index::set(-2, -2),
		                  Index::set(gdata.mx[dir0]+1, gdata.mx[dir1]+1));
		dEzdxim3_4.resize(Index::set(-2, -2),
		                  Index::set(gdata.mx[dir0]+1, gdata.mx[dir1]+1));
		dEzdyjm1_4.resize(Index::set(-2, -2),
		                  Index::set(gdata.mx[dir0]+1, gdata.mx[dir1]+1));
		dEzdyjm3_4.resize(Index::set(-2, -2),
		                  Index::set(gdata.mx[dir0]+1, gdata.mx[dir1]+1));
#endif

	}
#endif

	if(dir0 == 0 || dir1 == 0) {
		v_x.resize(Index::set(-2, -2),
		           Index::set(gdata.mx[dir0]+1, gdata.mx[dir1]+1));
		B_x.resize(Index::set(-2, -2),
		           Index::set(gdata.mx[dir0]+1, gdata.mx[dir1]+1));
		dBxdx.resize(Index::set(-2, -2),
		             Index::set(gdata.mx[dir0]+1, gdata.mx[dir1]+1));

		for(int num=0; num<2; ++num) {
			dvxdx[num].resize(Index::set(-2, -2),
			                  Index::set(gdata.mx[dir0]+1, gdata.mx[dir1]+1));
		}
	}
   
	if(dir0 == 1 || dir1 == 1) {
		v_y.resize(Index::set(-2, -2),
		           Index::set(gdata.mx[dir0]+1, gdata.mx[dir1]+1));
		B_y.resize(Index::set(-2, -2),
		           Index::set(gdata.mx[dir0]+1, gdata.mx[dir1]+1));
		dBydx.resize(Index::set(-2, -2),
		             Index::set(gdata.mx[dir0]+1, gdata.mx[dir1]+1));

		for(int num=0; num<2; ++num) {
			dvydx[num].resize(Index::set(-2, -2),
			                  Index::set(gdata.mx[dir0]+1, gdata.mx[dir1]+1));
		}
	}

	if(dir0 == 2 || dir1 == 2) {
		v_z.resize(Index::set(-2, -2),
		           Index::set(gdata.mx[dir0]+1, gdata.mx[dir1]+1));
		B_z.resize(Index::set(-2, -2),
		           Index::set(gdata.mx[dir0]+1, gdata.mx[dir1]+1));
		dBzdx.resize(Index::set(-2, -2),
		             Index::set(gdata.mx[dir0]+1, gdata.mx[dir1]+1));

		for(int num=0; num<2; ++num) {
			dvzdx[num].resize(Index::set(-2, -2),
			                  Index::set(gdata.mx[dir0]+1, gdata.mx[dir1]+1));
		}
	}

#if (CT_TYPE == CONSISTENT)
	for(int num=0; num<2; ++num) {

		v_ch_p[num].resize(Index::set(-2, -2),
		                   Index::set(gdata.mx[dir0]+1, gdata.mx[dir1]+1));
		v_ch_m[num].resize(Index::set(-2, -2),
		                   Index::set(gdata.mx[dir0]+1, gdata.mx[dir1]+1));
		v_cor_p[num].resize(Index::set(-2, -2),
		                    Index::set(gdata.mx[dir0]+1, gdata.mx[dir1]+1));
		v_cor_m[num].resize(Index::set(-2, -2),
		                    Index::set(gdata.mx[dir0]+1, gdata.mx[dir1]+1));

		v_ch_p[num].set_constVal(0.);
		v_ch_m[num].set_constVal(0.);
		v_cor_p[num].set_constVal(0.);
		v_cor_m[num].set_constVal(0.);

	}
#endif
	emf.resize(Index::set(-2, -2),
	           Index::set(gdata.mx[dir0]+1, gdata.mx[dir1]+1));

}


fields_2D::~fields_2D() {
	//! Destructor
	delete [] fluxWE;
	delete [] fluxSN;
	delete [] fluxBT;
}

phys_fields_1D::phys_fields_1D(const Data &gdata, const int &dir) 
	: pthermSYCL (Buffer<REAL, 1>(Range<1>(gdata.mx[dir] + 4))),
	  ptotalSYCL (Buffer<REAL, 1>(Range<1>(gdata.mx[dir] + 4)))
{

#if (FLUID_TYPE == CRONOS_MULTIFLUID)
	int n_omInt = gdata.fluids->get_N_OMINT_ALL();
#else
	int n_omInt = gdata.fluid.get_N_OMINT();
#endif

	uPriORIG = new NumMatrix<REAL,1> [n_omInt];
	auto size = gdata.mx[dir] + 1 + 2 + 1;
	uPriSYCL = std::vector(n_omInt, Buffer<REAL, 1>(Range<1>(size)));
	for (auto& e : uPriSYCL) {
		clearBufferOnHost(e);
		//auto e_acc = e.get_access<cl::sycl::access::mode::discard_write>();
		//for (int i = 0; i < size; ++i) {
		//	e_acc[i] = 0.0;
		//}
	}
	uConORIG = new NumMatrix<REAL, 1>[n_omInt];
	uConSYCL = std::vector(n_omInt, Buffer<REAL, 1>(Range<1>(size)));
	for (auto& e : uConSYCL) {
		clearBufferOnHost(e);
		//auto e_acc = e.get_access<cl::sycl::access::mode::discard_write>();
		//for (int i = 0; i < size; ++i) {
		//	e_acc[i] = 0.0;
		//}
	}
	fluxORIG = new NumMatrix<REAL,1> [n_omInt];
	fluxSYCL = std::vector(n_omInt, Buffer<REAL, 1>(Range<1>(size)));
	for (auto& e : fluxSYCL) {
		clearBufferOnHost(e);
		//auto e_acc = e.get_access<cl::sycl::access::mode::discard_write>();
		//for (int i = 0; i < size; ++i) {
		//	e_acc[i] = 0.0;
		//}
	}

	pthermORIG.resize(Index::set(-2),Index::set(gdata.mx[dir]+1));
	pthermSYCL = Buffer<REAL, 1>(Range<1>(size));

	ptotalORIG.resize(Index::set(-2), Index::set(gdata.mx[dir] + 1));
	ptotalSYCL = Buffer<REAL, 1>(Range<1>(size));
	
#if (USE_COROTATION == CRONOS_ON)
	// Prepare fields holding velocity in inertial frame
	uInertial[0].resize(Index::set(-2),Index::set(gdata.mx[dir]+1));
	uInertial[1].resize(Index::set(-2),Index::set(gdata.mx[dir]+1));
	uInertial[2].resize(Index::set(-2),Index::set(gdata.mx[dir]+1));
#endif

	for(int q=0; q<n_omInt; ++q) {
		uPriORIG[q].resize(Index::set(-2),Index::set(gdata.mx[dir]+1));
		uPriORIG[q].clear();
		uConORIG[q].resize(Index::set(-2),Index::set(gdata.mx[dir]+1));
		uConORIG[q].clear();

		fluxORIG[q].resize(Index::set(-2),Index::set(gdata.mx[dir]+1));
		fluxORIG[q].clear();
	}
}


phys_fields_1D::phys_fields_1D(const Data &gdata, int dir, int num) 
	: pthermSYCL (Buffer<REAL, 1>(Range<1>(gdata.mx[dir] + 4))),
	  ptotalSYCL (Buffer<REAL, 1>(Range<1>(gdata.mx[dir] + 4)))
{

	uPriORIG = new NumMatrix<REAL,1> [num];
	uConORIG = new NumMatrix<REAL,1> [num];
	fluxORIG = new NumMatrix<REAL,1> [num];

	for(int q=0; q<num; ++q) {
		uPriORIG[q].resize(Index::set(-2),Index::set(gdata.mx[dir]+1));
		uPriORIG[q].clear();
		uConORIG[q].resize(Index::set(-2),Index::set(gdata.mx[dir]+1));
		uConORIG[q].clear();

		fluxORIG[q].resize(Index::set(-2),Index::set(gdata.mx[dir]+1));
		fluxORIG[q].clear();
	}

}

phys_fields_1D::~phys_fields_1D() {
	delete [] uPriORIG;
	delete [] uConORIG;
	delete [] fluxORIG;
}

phys_fields_2D::phys_fields_2D(const Data &gdata, const int &dir) {

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

	if(dir == 1 || dir == 2) {
		v_x.resize(Index::set(-1, -1),
		           Index::set(gdata.mx[dir0]+1, gdata.mx[dir1]+1));
	}
	if(dir == 0 || dir == 2) {
		v_y.resize(Index::set(-1, -1),
		           Index::set(gdata.mx[dir0]+1, gdata.mx[dir1]+1));
	}
	if(dir == 0 || dir == 1) {
		v_z.resize(Index::set(-1, -1),
		           Index::set(gdata.mx[dir0]+1, gdata.mx[dir1]+1));
	}

	if(dir == 1 || dir == 2) {
		B_x.resize(Index::set(-1, -1),
		           Index::set(gdata.mx[dir0]+1, gdata.mx[dir1]+1));
	}
	if(dir == 0 || dir == 2) {
		B_y.resize(Index::set(-1, -1),
		           Index::set(gdata.mx[dir0]+1, gdata.mx[dir1]+1));
	}
	if(dir == 0 || dir == 1) {
		B_z.resize(Index::set(-1, -1),
		           Index::set(gdata.mx[dir0]+1, gdata.mx[dir1]+1));
	}

	emf.resize(Index::set(-1, -1),
	           Index::set(gdata.mx[dir0]+1, gdata.mx[dir1]+1));
  
}


phys_fields_0D::phys_fields_0D(const Data &gdata, int face, const CronosFluid &fluid) {

#if (FLUID_TYPE == CRONOS_MULTIFLUID)
//	fieldLists List(fluid, dir);
	this->num = fluid.get_N_OMINT();
#else
	this->num = N_OMINT;
#endif

	this->face = face;

	uPri.resize(num);
	uCon.resize(num);
	flux_phys.resize(num);
#if (USE_COROTATION == CRONOS_ON)
	uInertial.resize(num);
#endif

	ptotal = 0.;
	ptherm = 0.;
	carbuncle_flag = 0;

}

int phys_fields_0D::get_num() const {
	return num;
}

num_fields_0D::num_fields_0D(const Data &, const CronosFluid &fluid) {
#if (FLUID_TYPE == CRONOS_MULTIFLUID)
//	fieldLists List(fluid, dir);
	this->num = fluid.get_N_OMINT();
#else
	this->num = N_OMINT;
#endif

	v_ch_p = 1.;
	v_ch_m = 1.;
	flux_num.resize(num);

}
