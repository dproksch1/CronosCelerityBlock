#include "ConstrainedTransport.H"
#include <stdlib.h>
#include "minmod.H"

ConstrainedTransport::ConstrainedTransport(const Data &gdata, int dir,
                                           bool IntegrateA, int i_magFluid) {
	assert(dir >= 0 && dir < DIM);

	this->dir = dir;

	this->IntegrateA = IntegrateA;

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

	// Set field indices
#if(FLUID_TYPE == CRONOS_MULTIFLUID)
	cout << " getting ind " << endl;
	q_sx = gdata.fluids->get_q_sx(i_magFluid);
	q_sy = gdata.fluids->get_q_sy(i_magFluid);
	q_sz = gdata.fluids->get_q_sz(i_magFluid);
	q_Bx = gdata.fluids->get_q_Bx();
	q_By = gdata.fluids->get_q_By();
	q_Bz = gdata.fluids->get_q_Bz();
#else
	q_sx = gdata.fluid.get_q_sx();
	q_sy = gdata.fluid.get_q_sy();
	q_sz = gdata.fluid.get_q_sz();
	q_Bx = gdata.fluid.get_q_Bx();
	q_By = gdata.fluid.get_q_By();
	q_Bz = gdata.fluid.get_q_Bz();
#endif

	physValEmfLL = new phys_fields_2D(gdata, dir);
	physValEmfLR = new phys_fields_2D(gdata, dir);
	physValEmfRL = new phys_fields_2D(gdata, dir);
	physValEmfRR = new phys_fields_2D(gdata, dir);

	// Currently, only second order reconstruction available
	ReconstEmf = new Reconstruction2D_2nd(gdata, dir);

	fieldsEmf = new fields_2D(gdata, dir);

	DissMHD = new DissipationMHD;

	// if(value_exists("Include_CoordinateAxis")) {
	// 	// Check if grid on axis:
	// 	if((gdata.get_EdgeGridding() && 
	// 	    (gdata.global_xb[0] < 0.1*gdata.dx[0])) ||
	// 	   (!gdata.get_EdgeGridding() &&
	// 	    (gdata.global_xb[0] < 0.6*gdata.dx[0]))) {
	// 		with_axis = true;
	// 	} else {
	// 		cerr << " ERROR: Grid does not end on-axis. " << endl;
	// 		cerr << "   switch off /Include_CoordinateAxis/ for code to be allowed to run " << endl;
	// 		exit(3);
	// 	}
	// } else {
	// 	with_axis = false;
	// }
  
}

ConstrainedTransport::~ConstrainedTransport() {
	delete physValEmfLL;
	delete physValEmfLR;
	delete physValEmfRL;
	delete physValEmfRR;

	delete ReconstEmf;

	delete fieldsEmf;

	delete DissMHD;
}


void ConstrainedTransport::get_LocalArrays(const Data &gdata, int layer)
{
#if (FLUID_TYPE != CRONOS_HYDRO)

	if(dir == 0) {
		for (int k = -2; k <= gdata.mx[2]+1; ++k){
			for (int j = -2; j <= gdata.mx[1]+1; ++j){
				fieldsEmf->v_y(j,k) = gdata.om[q_sy](layer,j,k);
				fieldsEmf->v_z(j,k) = gdata.om[q_sz](layer,j,k);
				
				fieldsEmf->B_y(j,k) = gdata.om[q_By](layer,j,k);
				fieldsEmf->B_z(j,k) = gdata.om[q_Bz](layer,j,k);
			}
		}
	} else if (dir == 1) {
		for (int k = -2; k <= gdata.mx[2]+1; ++k){
			for (int i = -2; i <= gdata.mx[0]+1; ++i){
				fieldsEmf->v_x(k,i) = gdata.om[q_sx](i,layer,k);
				fieldsEmf->v_z(k,i) = gdata.om[q_sz](i,layer,k);
					
				fieldsEmf->B_x(k,i) = gdata.om[q_Bx](i,layer,k);
				fieldsEmf->B_z(k,i) = gdata.om[q_Bz](i,layer,k);
			}
		}
	} else {
		for (int j = -2; j <= gdata.mx[1]+1; ++j){
			for (int i = -2; i <= gdata.mx[0]+1; ++i){
				fieldsEmf->v_x(i,j) = gdata.om[q_sx](i,j,layer);
				fieldsEmf->v_y(i,j) = gdata.om[q_sy](i,j,layer);
          
				fieldsEmf->B_x(i,j) = gdata.om[q_Bx](i,j,layer);
				fieldsEmf->B_y(i,j) = gdata.om[q_By](i,j,layer);
			}
		}
	}
#endif
}

#if (USE_COROTATION == CRONOS_ON)
void ConstrainedTransport::TransInertToCorot(Data &gdata, Transformations &Trafo, int layer) {
	//! Transform to corotating frame

	if(dir == 0) {
		for (int k = -2; k <= gdata.mx[2]+1; ++k){
			for (int j = -2; j <= gdata.mx[1]+1; ++j){
				fieldsEmf->v_y(j,k) = Trafo.TransInertToCorot_y(gdata, fieldsEmf->v_y(j,k), layer, j,k);
				fieldsEmf->v_z(j,k) = Trafo.TransInertToCorot_z(gdata, fieldsEmf->v_z(j,k), layer, j,k);
			}
		}
	} else if (dir == 1) {
		for (int k = -2; k <= gdata.mx[2]+1; ++k){
			for (int i = -2; i <= gdata.mx[0]+1; ++i){
				fieldsEmf->v_x(k,i) = Trafo.TransInertToCorot_x(gdata, fieldsEmf->v_x(k,i), i,layer,k);
				fieldsEmf->v_z(k,i) = Trafo.TransInertToCorot_z(gdata, fieldsEmf->v_z(k,i), i,layer,k);

			}
		}
	} else {
		for (int j = -2; j <= gdata.mx[1]+1; ++j){
			for (int i = -2; i <= gdata.mx[0]+1; ++i){
				fieldsEmf->v_x(i,j) = Trafo.TransInertToCorot_x(gdata,fieldsEmf->v_x(i,j),i,j,layer);
				fieldsEmf->v_y(i,j) = Trafo.TransInertToCorot_y(gdata,fieldsEmf->v_y(i,j),i,j,layer);
			}
		}
	}

}
#endif


//#if (USE_COROTATION == CRONOS_ON)
//void ConstrainedTransport::get_NumEMFCT(Data &gdata, ProblemType &Problem,
//                                        Saves &Save,  Transformations &Trafo,
//                                        NumMatrix<REAL, 3> nom [N_OMINT])
//#else
void ConstrainedTransport::get_NumEMFCT(Data &gdata, ProblemType &Problem,
                                        Saves &Save,
                                        NumMatrix<REAL, 3> nom[N_OMINT])
//#endif
{
	
	for (int i = -1; i <= gdata.mx[dir]; ++i){
		ipos.set(dir,i);
		

		// Convert 3D data to 2D slabs
		get_LocalArrays(gdata, i);


		// Retrieve saved data
		Save.retrieve(gdata, *fieldsEmf, dir, i);


		// Doing the reconstruction - if necessary:
		Reconstruct(gdata);

//#if (USE_COROTATION == CRONOS_ON)
//		TransInertToCorot(gdata, Trafo, i);
//#endif

		// Get physical emfs - if necessary:
		get_PhysEmfs(gdata, Problem, *physValEmfLL, ipos);
		get_PhysEmfs(gdata, Problem, *physValEmfLR, ipos);
		get_PhysEmfs(gdata, Problem, *physValEmfRL, ipos);
		get_PhysEmfs(gdata, Problem, *physValEmfRR, ipos);


		// Get characteristic velocities - if necessary:
		get_vChar2D(gdata, Problem, ipos);


		// Get numerical form of emf (necessary routine):
		get_NumEmf2D(gdata, *physValEmfLL, *physValEmfLR, *physValEmfRL,
		             *physValEmfRR, *fieldsEmf);

		bc_EMF(gdata, Problem, *fieldsEmf, dir);

#ifdef PHYSDISS
		// Compute changes by dissipation
		DissMHD->modifyEmf(gdata, Problem, *fieldsEmf, dir, i);
#endif

		// Assign emf:
		get_ChangesEmf(gdata, nom, i);

		// Apply boundary conditions for nom (if necessary)
		bc_noms(gdata, Problem, dir);

	}
}




void ConstrainedTransport::get_ChangesEmf(Data &gdata,
                                          NumMatrix<REAL, 3> nom [N_OMINT],
                                          int layer) {

	if(IntegrateA) {

		if(layer >= 0 && layer <= gdata.mx[dir]) {

			if(dir == 0) {

				for (int k = 0; k <= gdata.mx[2]+1; ++k){
					for (int j = 0; j <= gdata.mx[1]+1; ++j){
						
						nom[q_Bx](layer,j,k) = fieldsEmf->emf(j,k);
						
					}
				}
				
			} else if (dir == 1) {

				// Beware dimensions 1 and 3 are swapped(!)
				for (int k = 0; k <= gdata.mx[2]+1; ++k){
					for (int i = 0; i <= gdata.mx[0]+1; ++i){
						
						nom[q_By](i,layer,k) = fieldsEmf->emf(k,i);
						
					}
				}

			} else {
			
				for (int j = 0; j <= gdata.mx[1]+1; ++j){
					for (int i = 0; i <= gdata.mx[0]+1; ++i){
						
						nom[q_Bz](i,j,layer) = fieldsEmf->emf(i,j);

					}
				}

			}

		}

	} else {

		if(layer >= 0 && layer <= gdata.mx[dir]) {
			
			if(dir == 0) {

#if (NON_LINEAR_GRID == CRONOS_OFF)
				REAL idy = gdata.idx[1];
				REAL idz = gdata.idx[2];
#endif

				// Emf = Ex in this case:
				double i(1.*layer);
				for (int k = 0; k <= gdata.mx[2]; ++k){
#if (NON_LINEAR_GRID == CRONOS_ON)
					REAL idz = gdata.getCen_idx(2,k);
#endif
					for (int j = -1; j <= gdata.mx[1]; ++j){


#if (GEOM != CARTESIAN)

#if (NON_LINEAR_GRID == CRONOS_OFF)
						REAL idA = 1./(gdata.h0(i,j+0.5,k)*gdata.h2(i,j+0.5,k));

						nom[q_By](layer,j,k) += (gdata.h0(i,j+0.5,k+0.5)*fieldsEmf->emf(j+1,k+1) -
								gdata.h0(i,j+0.5,k-0.5)*fieldsEmf->emf(j+1,k  ))*idz*idA;

#else
						REAL idA = 1./(gdata.h0(i,j,k,0,1,0)*
						               gdata.h2(i,j,k,0,1,0));
						
						nom[q_By](layer,j,k) += (gdata.h0(i,j,k,0,1, 1)*fieldsEmf->emf(j+1,k+1) -
								gdata.h0(i,j,k,0,1,-1)*fieldsEmf->emf(j+1,k  ))*idz*idA;
		
#endif

#else
						nom[q_By](layer,j,k) +=  (fieldsEmf->emf(j+1,k+1) -
						                       fieldsEmf->emf(j+1,k  ))*idz;

#endif
						
					}
				}

				
				for (int k = -1; k <= gdata.mx[2]; ++k){
					for (int j = 0; j <= gdata.mx[1]; ++j){
#if (NON_LINEAR_GRID == CRONOS_ON)
						REAL idy = gdata.getCen_idx(1,j);
#endif
					
#if (GEOM != CARTESIAN)

#if (NON_LINEAR_GRID == CRONOS_OFF)

						REAL idA(1./(gdata.h0(i,j,k+0.5)*gdata.h1(i,j,k+0.5)));

						nom[q_Bz](layer,j,k) += -(fieldsEmf->emf(j+1,k+1)*gdata.h0(i,j+0.5,k+0.5) -
						                       fieldsEmf->emf(j  ,k+1)*gdata.h0(i,j-0.5,k+0.5))*idy*idA;

#else

						REAL idA(1./(gdata.h0(i,j,k,0,0,1)*
						             gdata.h1(i,j,k,0,0,1)));

						nom[q_Bz](layer,j,k) += -(fieldsEmf->emf(j+1,k+1)*gdata.h0(i,j,k,0, 1,1) -
						                       fieldsEmf->emf(j  ,k+1)*gdata.h0(i,j,k,0,-1,1))*idy*idA;

#endif

#else

						nom[q_Bz](layer,j,k) += -(fieldsEmf->emf(j+1,k+1) -
						                       fieldsEmf->emf(j  ,k+1))*idy;

#endif

					}
				}

			} else if (dir == 1) {
				
#if (NON_LINEAR_GRID == CRONOS_OFF)
				REAL idx = gdata.idx[0];
				REAL idz = gdata.idx[2];
#endif				

				// Emf = Ey in this case:

				// Beware dimensions 1 and 3 are swapped(!)
				double j(1.*layer);
				for (int k = 0; k <= gdata.mx[2]; ++k){
#if (NON_LINEAR_GRID == CRONOS_ON)
					REAL idz(gdata.getCen_idx(2,k));
#endif
					for (int i = -1; i <= gdata.mx[0]; ++i){
						
#if (GEOM != CARTESIAN)

#if (NON_LINEAR_GRID == CRONOS_OFF)
						REAL idA(1./(gdata.h1(i+0.5,j,k)*gdata.h2(i+0.5,j,k)));
						
						nom[q_Bx](i,layer,k) += -(fieldsEmf->emf(k+1,i+1)*gdata.h1(i+0.5,j,k+0.5) -
						                       fieldsEmf->emf(k  ,i+1)*gdata.h1(i+0.5,j,k-0.5))*idz*idA;

#else
						REAL idA(1./(gdata.h1(i,j,k,1,0,0)*
						             gdata.h2(i,j,k,1,0,0)));

						nom[q_Bx](i,layer,k) += -(gdata.h1(i,j,k,1,0, 1)*fieldsEmf->emf(k+1,i+1) -
						                       gdata.h1(i,j,k,1,0,-1)*fieldsEmf->emf(k  ,i+1))*idz*idA;

#endif

#else

						nom[q_Bx](i,layer,k) += -(fieldsEmf->emf(k+1,i+1) -
						                       fieldsEmf->emf(k  ,i+1))*idz;

#endif

					}
				}

				// Beware dimensions 1 and 3 are swapped(!)
				for (int k = -1; k <= gdata.mx[2]; ++k){
					for (int i = 0; i <= gdata.mx[0]; ++i){
#if (NON_LINEAR_GRID == CRONOS_ON)
						REAL idx(gdata.getCen_idx(0,i));
#endif
						
#if (GEOM != CARTESIAN)

#if (NON_LINEAR_GRID == CRONOS_OFF)

						REAL idA(1./(gdata.h0(i,j,k+0.5)*gdata.h1(i,j,k+0.5)));

						nom[q_Bz](i,layer,k) +=  (fieldsEmf->emf(k+1,i+1)*gdata.h1(i+0.5,j,k+0.5) -
						                       fieldsEmf->emf(k+1,i  )*gdata.h1(i-0.5,j,k+0.5))*idx*idA;

#else
						REAL idA(1./(gdata.h0(i,j,k,0,0,1)*
						             gdata.h1(i,j,k,0,0,1)));

						nom[q_Bz](i,layer,k) += (gdata.h1(i,j,k, 1,0,1)*fieldsEmf->emf(k+1,i+1) -
						                      gdata.h1(i,j,k,-1,0,1)*fieldsEmf->emf(k+1,i  ))*idx*idA;

#endif

#else

						nom[q_Bz](i,layer,k) +=  (fieldsEmf->emf(k+1,i+1) -
						                       fieldsEmf->emf(k+1,i  ))*idx;

#endif
						
					}
				}

			} else {

#if (NON_LINEAR_GRID == CRONOS_OFF)
				REAL idx = gdata.idx[0];
				REAL idy = gdata.idx[1];
#endif

				// Emf = Ez in this case

				double k(1.*layer);
				for (int j = 0; j <= gdata.mx[1]; ++j){
#if (NON_LINEAR_GRID == CRONOS_ON)
					REAL idy(gdata.getCen_idx(1,j));
#endif
					for (int i = -1; i <= gdata.mx[0]; ++i){
						
#if (GEOM != CARTESIAN)

#if (NON_LINEAR_GRID == CRONOS_OFF)
						REAL idA(1./(gdata.h1(i+0.5,j,k)*gdata.h2(i+0.5,j,k)));

						nom[q_Bx](i,j,layer) +=  (fieldsEmf->emf(i+1,j+1)*gdata.h2(i+0.5,j+0.5,k) -
						                       fieldsEmf->emf(i+1,j  )*gdata.h2(i+0.5,j-0.5,k))*idy*idA;

#else
						REAL idA(1./(gdata.h1(i,j,k,1,0,0)*
						             gdata.h2(i,j,k,1,0,0)));

						nom[q_Bx](i,j,layer) += (gdata.h2(i,j,k,1, 1,0)*fieldsEmf->emf(i+1,j+1) -
						                      gdata.h2(i,j,k,1,-1,0)*fieldsEmf->emf(i+1,j  ))*idy*idA;

#endif

#else

						nom[q_Bx](i,j,layer) +=  (fieldsEmf->emf(i+1,j+1) -
						                       fieldsEmf->emf(i+1,j  ))*idy;

#endif
					}
				}

				for (int j = -1; j <= gdata.mx[1]; ++j){
					for (int i = 0; i <= gdata.mx[0]; ++i){
#if (NON_LINEAR_GRID == CRONOS_ON)
						REAL idx(gdata.getCen_idx(0,i));
#endif

#if (GEOM != CARTESIAN)

#if (NON_LINEAR_GRID == CRONOS_OFF)

						REAL idA(1./(gdata.h0(i,j+0.5,k)*gdata.h2(i,j+0.5,k)));

						nom[q_By](i,j,layer) += -(fieldsEmf->emf(i+1,j+1)*gdata.h2(i+0.5,j+0.5,k) -
						                       fieldsEmf->emf(i  ,j+1)*gdata.h2(i-0.5,j+0.5,k))*idx*idA;

#else
						REAL idA(1./(gdata.h0(i,j,k,0,1,0)*
						             gdata.h2(i,j,k,0,1,0)));
						
						nom[q_By](i,j,layer) += -(gdata.h2(i,j,k, 1,1,0)*fieldsEmf->emf(i+1,j+1) -
						                       gdata.h2(i,j,k,-1,1,0)*fieldsEmf->emf(i  ,j+1))*idx*idA;

#endif

#else

						nom[q_By](i,j,layer) += -(fieldsEmf->emf(i+1,j+1) -
						                       fieldsEmf->emf(i  ,j+1))*idx;

#endif

					}
				}

			}
		}

	}
}
