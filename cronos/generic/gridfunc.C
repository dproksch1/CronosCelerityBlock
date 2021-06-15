#include "gridgen.H"
#include "data.H"
#include "problem.H"
#include "util.H"
#include <iomanip>
#include <matrix.H>


/*-----------------------------------------------------------------
  Data connected functions
  -----------------------------------------------------------------*/


gridFunc::gridFunc(Data &gdata)
{


#if (FLUID_TYPE == CRONOS_MULTIFLUID)
	n_om = gdata.fluids->get_N_OM();
	n_omInt = gdata.fluids->get_N_OMINT();
	n_omUser = gdata.fluids->get_N_OM_USER();
	n_omIntUser = gdata.fluids->get_N_OMINT_USER();
	n_omIntAll = gdata.fluids->get_N_OMINT_ALL();

	q_Bx = gdata.fluids->get_q_Bx();
	q_By = gdata.fluids->get_q_By();
	q_Bz = gdata.fluids->get_q_Bz();
#else
#if (FLUID_TYPE == CRONOS_MHD)
	q_Bx = gdata.fluid.get_q_Bx();
	q_By = gdata.fluid.get_q_By();
	q_Bz = gdata.fluid.get_q_Bz();
#else
	// Use a value that is not used (-1 used by user fields -> so, I am using -2)
	q_Bx = -2;
	q_By = -2;
	q_Bz = -2;
#endif

	n_om = N_OM;
	n_omInt = N_OMINT;
#if (OMS_USER == TRUE)
	n_omUser = N_OM_USER;
	n_omIntUser = N_OMINT_USER;
#else
	n_omIntUser = 0;
#endif
	n_omIntAll = N_OMINT_ALL;
#endif

#ifndef parallel
	bc_Type[0] = int(value((char*)"bc_x_bot"));
	bc_Type[1] = int(value((char*)"bc_x_top"));
	bc_Type[2] = int(value((char*)"bc_y_bot"));
	bc_Type[3] = int(value((char*)"bc_y_top"));
	bc_Type[4] = int(value((char*)"bc_z_bot"));
	bc_Type[5] = int(value((char*)"bc_z_top"));

	// if necessary: compute positions for partners on grid axis (so
	// far only for cylindrical grids
#if (GEOM == CYLINDRICAL)
	if(gdata.get_singularity_treatment(0) > 0) {
		compute_AxisPartners(gdata);
		bc_Type[0] = 6;
	}
#elif (GEOM == SPHERICAL)
	if(gdata.get_singularity_treatment(2) > 0) {
		compute_AxisPartners(gdata, 0); // Partners at bottom
		bc_Type[2] = 6;
	}
	if(gdata.get_singularity_treatment(3) > 0) {
		compute_AxisPartners(gdata, 1); // Partners at top
		bc_Type[3] = 6;
	}

#endif

#else
	bc_Type[0] = int(value((char*)"bc_x_bot"));
	bc_Type[1] = int(value((char*)"bc_x_top"));
	bc_Type[2] = int(value((char*)"bc_y_bot"));
	bc_Type[3] = int(value((char*)"bc_y_top"));
	bc_Type[4] = int(value((char*)"bc_z_bot"));
	bc_Type[5] = int(value((char*)"bc_z_top"));

	for(int bound=0; bound<6; ++bound) {
		if(bc_Type[bound] == 6) {
			cerr << " boundary type 6 not to be assigned by user " << endl;
			cerr << " EXITING " << endl;
#ifdef parallel
			MPI_Finalize();
#endif
			exit(-29);
		}
	}



// 	bool atAxis = false;
// #if (GEOM == CYLINDRICAL)
// 	REAL eps(1.e-5*(gdata.global_xe[0]-gdata.global_xb[0]));
// 	bc_Type[0] = -1;
// 	if(gdata.xb[0] >= gdata.global_xb[0]-eps &&
// 	   gdata.xb[0] <= gdata.global_xb[0]+eps) {
// 		atAxis = true;
// 	}
// #elif (GEOM == SPHERICAL)
// 	REAL eps(1.e-5*(gdata.global_xe[1]-gdata.global_xb[1]));
// 	if((gdata.xb[1] >= gdata.global_xb[1]-eps &&
// 	    gdata.xb[1] <= gdata.global_xb[1]+eps) ||
// 	   (gdata.xe[1] >= gdata.global_xe[1]-eps &&
// 	    gdata.xe[1] <= gdata.global_xe[1]+eps)) {
// 		atAxis = true;
// 	}
// #endif


	// if(atAxis) {
#if (GEOM == CYLINDRICAL)
	REAL eps(1.e-5*(gdata.global_xe[0]-gdata.global_xb[0]));
		bc_Type[0] = int(value((char*)"bc_x_bot"));
		if(gdata.get_singularity_treatment(0) > 0) {
			compute_AxisPartners(gdata);
			//			if(gdata.axis_treatment(gdata.phi_dir) == 1) {
			// if(gdata.axis_treatment(0) == 1) {
			if(gdata.get_singularity_treatment(0)==1) {
				bc_Type[0] = 6;
			}
		}
#elif (GEOM == SPHERICAL)
		REAL eps(1.e-5*(gdata.global_xe[1]-gdata.global_xb[1]));
		// if(gdata.get_singularity_treatment(2) > 0 ||
		// 		gdata.get_singularity_treatment(3) > 0) {
		// 	compute_AxisPartners(gdata,0);
		// 	compute_AxisPartners(gdata,1);
		// }
		if(gdata.get_singularity_treatment(2) == 1) {
			compute_AxisPartners(gdata,0);
			cout << " involved " << gdata.rank << endl;
		}

		if(gdata.get_singularity_treatment(3) == 1) {
			cout << " doing this " << gdata.rank << endl;
			compute_AxisPartners(gdata,1);
		}


		// if(gdata.get_singularity_treatment(2) > 0) {
		if(gdata.get_singularity_treatment(2) == 1) {
			bc_Type[2] = 6;
		}

		// if(gdata.get_singularity_treatment(3) > 0) {
		if(gdata.get_singularity_treatment(3) == 1) {
			bc_Type[3] = 6;
		}
#else
		REAL eps(0.);
#endif
	// }
	// Ranks that do not participate in axis-stuff have to wait...
//	MPI_Barrier(gdata.comm3d);


	// Write specific boundary conditions if not on axis
	if(bc_Type[0] != 6) {
//		eps = 1.e-5*(gdata.global_xe[0]-gdata.global_xb[0]);
//		bc_Type[0] = -1;
//		if(gdata.xb[0] >= gdata.global_xb[0]-eps &&
//		   gdata.xb[0] <= gdata.global_xb[0]+eps) {
//			bc_Type[0] = int(value((char*)"bc_x_bot"));
//		}

		bc_Type[0] = -1;
		if(gdata.check_MpiLowest(0)) {
			bc_Type[0] = int(value((char*)"bc_x_bot"));
		}

	}

//	bc_Type[1] = -1;
//	if(gdata.xe[0] >= gdata.global_xe[0]-eps &&
//	   gdata.xe[0] <= gdata.global_xe[0]+eps) {
//		bc_Type[1] = int(value((char*)"bc_x_top"));
//	}

	bc_Type[1] = -1;
	if(gdata.check_MpiHighest(0)) {
		bc_Type[1] = int(value((char*)"bc_x_top"));
	}


	if(bc_Type[2] != 6) {
//		eps = 1.e-5*(gdata.global_xe[1] - gdata.global_xb[1]);
//		bc_Type[2] = -1;
//		if(gdata.xb[1] >= gdata.global_xb[1]-eps &&
//		   gdata.xb[1] <= gdata.global_xb[1]+eps) {
//			bc_Type[2] = int(value((char*)"bc_y_bot"));
//		}
		bc_Type[2] = -1;
		if(gdata.check_MpiLowest(1)) {
			bc_Type[2] = int(value((char*)"bc_y_bot"));
		}
	}

	if(bc_Type[3] != 6) {
//		bc_Type[3] = -1;
//		if(gdata.xe[1] >= gdata.global_xe[1]-eps &&
//		   gdata.xe[1] <= gdata.global_xe[1]+eps) {
//			bc_Type[3] = int(value((char*)"bc_y_top"));
//		}
		bc_Type[3] = -1;
		if(gdata.check_MpiHighest(1)) {
			bc_Type[3] = int(value((char*)"bc_y_top"));
		}
	}

//	eps = 1.e-5*(gdata.global_xe[2]-gdata.global_xb[2]);
//	bc_Type[4] = -1;
//	if(gdata.xb[2] >= gdata.global_xb[2]-eps &&
//	   gdata.xb[2] <= gdata.global_xb[2]+eps) {
//		bc_Type[4] = int(value((char*)"bc_z_bot"));
//	}
	bc_Type[4] = -1;
	if(gdata.check_MpiLowest(2)) {
		bc_Type[4] = int(value((char*)"bc_z_bot"));
	}

//	bc_Type[5] = -1;
//	if(gdata.xe[2] >= gdata.global_xe[2]-eps &&
//	   gdata.xe[2] <= gdata.global_xe[2]+eps) {
//		bc_Type[5] = int(value((char*)"bc_z_top"));
//	}
	bc_Type[5] = -1;
	if(gdata.check_MpiHighest(2)) {
		bc_Type[5] = int(value((char*)"bc_z_top"));
	}
#endif

	// Set all bc types to fields:
	int bc_TypeLow[DIM], bc_TypeHigh[DIM];
	for(int dir=0; dir<DIM; ++dir) {
		bc_TypeLow[dir] = bc_Type[2*dir];
		bc_TypeHigh[dir] = bc_Type[2*dir+1];
	}

	for(int q=0; q<n_om+N_P; ++q) {
		gdata.om[q].set_bcType(bc_TypeLow, bc_TypeHigh);
	}
#if (OMS_USER == TRUE)
	for(int q=0; q<n_omUser; ++q) {
		gdata.om_user[q].set_bcType(bc_TypeLow, bc_TypeHigh);
	}
#endif

	// For constant boundaries set all fields that are not directly
	// accessed by the user to extrapolation
	for(int dir=0; dir<DIM; ++dir) {
		for(int q=n_omInt; q<n_om+N_P; ++q) {
			if(gdata.om[q].get_bcTypeLow(dir) == 7) {
				gdata.om[q].set_bcTypeLow(dir, 2);
			}
			if(gdata.om[q].get_bcTypeHigh(dir) == 7) {
				gdata.om[q].set_bcTypeHigh(dir, 2);
			}
		}
#if (OMS_USER == TRUE)
		for(int q=n_omIntUser; q<n_omUser; ++q) {
			if(gdata.om_user[q].get_bcTypeLow(dir) == 7) {
				gdata.om_user[q].set_bcTypeLow(dir, 2);
			}
			if(gdata.om_user[q].get_bcTypeHigh(dir) == 7) {
				gdata.om_user[q].set_bcTypeHigh(dir, 2);
			}
		}
#endif
	}
	int use_fixed_bcs = 0;
	for(int q=0; q<n_om+N_P; ++q) {
//		if(gdata.rank==0) {
//			cout << " BC types for field " << q << " ";
//			cout << gdata.om[q].get_bcTypeLow(0) << " ";
//			cout << gdata.om[q].get_bcTypeLow(1) << " ";
//			cout << gdata.om[q].get_bcTypeLow(2) << " ";
//			cout << gdata.om[q].get_bcTypeHigh(0) << " ";
//			cout << gdata.om[q].get_bcTypeHigh(1) << " ";
//			cout << gdata.om[q].get_bcTypeHigh(2) << " ";
//			cout << endl;
//		}
		// Check type of boundary condition:
		if(gdata.om[q].get_bcTypeLow(0)==7) use_fixed_bcs = 1;
		if(gdata.om[q].get_bcTypeLow(1)==7) use_fixed_bcs = 1;
		if(gdata.om[q].get_bcTypeLow(2)==7) use_fixed_bcs = 1;
		if(gdata.om[q].get_bcTypeHigh(0)==7) use_fixed_bcs = 1;
		if(gdata.om[q].get_bcTypeHigh(1)==7) use_fixed_bcs = 1;
		if(gdata.om[q].get_bcTypeHigh(2)==7) use_fixed_bcs = 1;
	}
#ifdef parallel
	int use_fixed_bcs_global;
	MPI_Allreduce(&use_fixed_bcs, &use_fixed_bcs_global, 1, MPI_INT, MPI_MAX,
			MPI_COMM_WORLD);
	use_fixed_bcs = use_fixed_bcs_global;
#endif

#if (VEC_POT_BCS == TRUE)
	if(use_fixed_bcs) {
		cerr << " gridfunc::gridfunc " << endl;
		cerr << " Error: fixed BCs together with vector potential boundaries is not implemented " << endl;
		exit(-22);
	}
#endif

	// Prepare size of Send and Receive arrays
	// (currently only for case of edge gridding)

	if(gdata.get_EdgeGridding()) {

		int rim = 3;

		int iminSend = 1;
		int imaxSend = rim;

//		Send_x.resize(Index::set(iminSend, -rim, -rim),
//				Index::set(imaxSend, gdata.mx[1]+rim, gdata.mx[2]+rim));
//		SendArr_x.resize((imaxSend-iminSend+1)*(gdata.mx[1]+2*rim+1)*(gdata.mx[2]+2*rim+1));
//		SendArrSmall_x.resize((imaxSend-iminSend)*(gdata.mx[1]+2*rim+1)*(gdata.mx[2]+2*rim+1));
//
//		Send_y.resize(Index::set(-rim, iminSend, -rim),
//				Index::set(gdata.mx[0]+rim, imaxSend, gdata.mx[2]+rim));
//		SendArr_y.resize((gdata.mx[0]+2*rim+1)*(imaxSend-iminSend+1)*(gdata.mx[2]+2*rim+1));
//		SendArrSmall_y.resize((gdata.mx[0]+2*rim+1)*(imaxSend-iminSend)*(gdata.mx[2]+2*rim+1));
//
//		Send_z.resize(Index::set(-rim, -rim, iminSend),
//				Index::set(gdata.mx[0]+rim, gdata.mx[1]+rim, imaxSend));
//		SendArr_z.resize((gdata.mx[0]+2*rim+1)*(gdata.mx[1]+2*rim+1)*(imaxSend-iminSend+1));
//		SendArrSmall_z.resize((gdata.mx[0]+2*rim+1)*(gdata.mx[1]+2*rim+1)*(imaxSend-iminSend));
//
//		Send_empty.resize(Index::set(0,0,0), Index::set(0,0,0));
//		SendArr_empty.resize(0);

		int size_x = (imaxSend-iminSend+1)*(gdata.mx[1]+2*rim+1)*(gdata.mx[2]+2*rim+1);
		int size_y = (gdata.mx[0]+2*rim+1)*(imaxSend-iminSend+1)*(gdata.mx[2]+2*rim+1);
		int size_z = (gdata.mx[0]+2*rim+1)*(gdata.mx[1]+2*rim+1)*(imaxSend-iminSend+1);
		int size = std::max(size_x, std::max(size_y, size_z));

		SendArr_buff.resize(size);
		RecvArr_buff.resize(size);

		SendArr_xRL.resize(size_x);
		SendArr_xLR.resize(size_x);
		SendArr_yRL.resize(size_y);
		SendArr_yLR.resize(size_y);
		SendArr_zRL.resize(size_z);
		SendArr_zLR.resize(size_z);

		RecvArr_xRL.resize(size_x);
		RecvArr_xLR.resize(size_x);
		RecvArr_yRL.resize(size_y);
		RecvArr_yLR.resize(size_y);
		RecvArr_zRL.resize(size_z);
		RecvArr_zLR.resize(size_z);

		SendArrRL.resize(size);
		SendArrLR.resize(size);

		RecvArrRL.resize(size);
		RecvArrLR.resize(size);

		// imaxRecv = rim-1;
		int iminRecv = 1; // Corresponding to mx[0]+1
		int imaxRecv = rim;
//		Recv_x.resize(Index::set(iminRecv,-rim,-rim),
//				Index::set(imaxRecv,gdata.mx[1]+rim,gdata.mx[2]+rim));
//		RecvArr_x.resize((imaxSend-iminSend+1)*(gdata.mx[1]+2*rim+1)*(gdata.mx[2]+2*rim+1));
//		RecvArrSmall_x.resize((imaxSend-iminSend)*(gdata.mx[1]+2*rim+1)*(gdata.mx[2]+2*rim+1));
//
//		Recv_y.resize(Index::set(-rim,iminRecv,-rim),
//				Index::set(gdata.mx[0]+rim,imaxRecv,gdata.mx[2]+rim));
//		RecvArr_y.resize((gdata.mx[0]+2*rim+1)*(imaxSend-iminSend+1)*(gdata.mx[2]+2*rim+1));
//		RecvArrSmall_y.resize((gdata.mx[0]+2*rim+1)*(imaxSend-iminSend)*(gdata.mx[2]+2*rim+1));
//
//		Recv_z.resize(Index::set(-rim,-rim,iminRecv),
//				Index::set(gdata.mx[0]+rim,gdata.mx[1]+rim,imaxRecv));
//		RecvArr_z.resize((gdata.mx[0]+2*rim+1)*(gdata.mx[1]+2*rim+1)*(imaxSend-iminSend+1));
//		RecvArrSmall_z.resize((gdata.mx[0]+2*rim+1)*(gdata.mx[1]+2*rim+1)*(imaxSend-iminSend));
//
//		Recv_empty.resize(Index::set(0,0,0), Index::set(0,0,0));
//		RecvArr_empty.resize(0);

	} else 	{ // old gridding (mx)

		int rim = 3;

//		int iminSend = 1;
		int iminSend = 0;
		int imaxSend = rim;

		int size_x = (imaxSend-iminSend+1)*(gdata.mx[1]+2*rim+1)*(gdata.mx[2]+2*rim+1);
		int size_y = (gdata.mx[0]+2*rim+1)*(imaxSend-iminSend+1)*(gdata.mx[2]+2*rim+1);
		int size_z = (gdata.mx[0]+2*rim+1)*(gdata.mx[1]+2*rim+1)*(imaxSend-iminSend+1);
		int size = std::max(size_x, std::max(size_y, size_z));

		SendArr_buff.resize(size);
		RecvArr_buff.resize(size);

		SendArr_xRL.resize(size_x);
		SendArr_xLR.resize(size_x);
		SendArr_yRL.resize(size_y);
		SendArr_yLR.resize(size_y);
		SendArr_zRL.resize(size_z);
		SendArr_zLR.resize(size_z);

		RecvArr_xRL.resize(size_x);
		RecvArr_xLR.resize(size_x);
		RecvArr_yRL.resize(size_y);
		RecvArr_yLR.resize(size_y);
		RecvArr_zRL.resize(size_z);
		RecvArr_zLR.resize(size_z);

		SendArrRL.resize(size);
		SendArrLR.resize(size);

		RecvArrRL.resize(size);
		RecvArrLR.resize(size);
	}




//	cout << " Done prep " << gdata.rank << endl;
	// for(int dir=0; dir<DIM; ++dir) {
	// 	for(int q=0; q<n_om+N_P; ++q) {
	// 		if(bc_Type[2*dir] == 7) {
	// 		} else {
	// 			om[q].set
	// 		}
	// 	}
	// }

	// for(int ibound=0; ibound<2*DIM; ++ibound) {
	// 	if(bc_Type[ibound] == 7) {
	// 		for(int q=0; q<n_omInt; ++q) {
	// 			om.q
	// 		}
	// 	}
	// }
	cout << " My upper type: " << gdata.rank << " " << bc_Type[2] << " " << bc_Type[3] << endl;

}



int gridFunc::get_bc_Type(int num) 
{
	int bc_type_ret(0);
	if(num < 8) {
		bc_type_ret = bc_Type[num];
	} else {
		cerr << " no such boundary: " << num << endl;
		exit(-32);
		bc_type_ret = -99;
	}
	return bc_type_ret;
}



void gridFunc::prep_boundaries(Data &gdata, ProblemType &pt) {

	if(gdata.rank==0) {
		cout << " Storing values for fixed boundaries " << endl;
	}

	if(bc_Type[0] == 7) { // x, bottom

		int lbound[3] = {-B,-B,-B};
		int ubound[3] = {-1,gdata.mx[1]+B,gdata.mx[2]+B};


		for(int q=0; q<n_omInt; ++q) {
			bcVals_xb[q].resize(lbound, ubound);

			// Setting field at lower x bound:
			for(int iz=-B; iz<=gdata.mx[2]+B; ++iz) {
				for(int iy=-B; iy<=gdata.mx[1]+B; ++iy) {
					for(int ix=-B; ix<=-1; ++ix) {
						bcVals_xb[q](ix, iy, iz) = gdata.om[q](ix, iy, iz);
					}
				}
			}
		}
#if (OMS_USER == TRUE)

		for(int q=0; q<n_omIntUser; ++q) {
			bcVals_User_xb[q].resize(lbound, ubound);

			// Setting user field at lower x bound:
			for(int iz=-B; iz<=gdata.mx[2]+B; ++iz) {
				for(int iy=-B; iy<=gdata.mx[1]+B; ++iy) {
					for(int ix=-B; ix<=-1; ++ix) {
						bcVals_User_xb[q](ix, iy, iz) =
							gdata.om_user[q](ix, iy, iz);
					}
				}
			}
		}
#endif
	}

	if(bc_Type[1] == 7) { // x, top

		int lbound[3] = {gdata.mx[0]+1,-B,-B};
		int ubound[3] = {gdata.mx[0]+B,gdata.mx[1]+B,gdata.mx[2]+B};

		for(int q=0; q<n_omInt; ++q) {
			bcVals_xe[q].resize(lbound, ubound);

			// Setting field at upper x bound:
			for(int iz=-B; iz<=gdata.mx[2]+B; ++iz) {
				for(int iy=-B; iy<=gdata.mx[1]+B; ++iy) {
					for(int ix=gdata.mx[0]+1; ix<=gdata.mx[0]+B; ++ix) {
						bcVals_xe[q](ix, iy, iz) = gdata.om[q](ix, iy, iz);
					}
				}
			}
		}
#if (OMS_USER == TRUE)
		for(int q=0; q<n_omIntUser; ++q) {
			bcVals_User_xe[q].resize(lbound, ubound);

			// Setting user field at upper x bound:
			for(int iz=-B; iz<=gdata.mx[2]+B; ++iz) {
				for(int iy=-B; iy<=gdata.mx[1]+B; ++iy) {
					for(int ix=gdata.mx[0]+1; ix<=gdata.mx[0]+B; ++ix) {
						bcVals_User_xe[q](ix, iy, iz) = 
							gdata.om_user[q](ix, iy, iz);
					}
				}
			}
		}
#endif
	}

	
	if(bc_Type[2] == 7) { // y, bottom

		int lbound[3] = {-B,-B,-B};
		int ubound[3] = {gdata.mx[1]+B,-1,gdata.mx[2]+B};

		for(int q=0; q<n_omInt; ++q) {
			bcVals_yb[q].resize(lbound, ubound);

			// Setting field at lower y bound:
			for(int iz=-B; iz<=gdata.mx[2]+B; ++iz) {
				for(int iy=-B; iy<=-1; ++iy) {
					for(int ix=-B; ix<=gdata.mx[1]+B; ++ix) {
						bcVals_yb[q](ix, iy, iz) = gdata.om[q](ix, iy, iz);
					}
				}
			}
		}
#if (OMS_USER == TRUE)
		for(int q=0; q<n_omIntUser; ++q) {
			bcVals_User_yb[q].resize(lbound, ubound);

			// Setting user field at lower y bound:
			for(int iz=-B; iz<=gdata.mx[2]+B; ++iz) {
				for(int iy=-B; iy<=-1; ++iy) {
					for(int ix=-B; ix<=gdata.mx[1]+B; ++ix) {
						bcVals_User_yb[q](ix, iy, iz) = 
							gdata.om_user[q](ix, iy, iz);
					}
				}
			}
		}
#endif
	}
	
	if(bc_Type[3] == 7) { // y, top

		int lbound[3] = {-B,gdata.mx[1]+1,-B};
		int ubound[3] = {gdata.mx[0]+B,gdata.mx[1]+B,gdata.mx[2]+B};

		for(int q=0; q<n_omInt; ++q) {
			bcVals_ye[q].resize(lbound, ubound);

			// Setting field at uppter y bound:
			for(int iz=-B; iz<=gdata.mx[2]+B; ++iz) {
				for(int iy=gdata.mx[1]+1; iy<=gdata.mx[1]+B; ++iy) {
					for(int ix=-B; ix<=gdata.mx[0]+B; ++ix) {
						bcVals_ye[q](ix, iy, iz) = gdata.om[q](ix, iy, iz);
					}
				}
			}
		}
#if (OMS_USER == TRUE)
		for(int q=0; q<n_omIntUser; ++q) {
			bcVals_User_ye[q].resize(lbound, ubound);

			// Setting user field at uppter y bound:
			for(int iz=-B; iz<=gdata.mx[2]+B; ++iz) {
				for(int iy=gdata.mx[1]+1; iy<=gdata.mx[1]+B; ++iy) {
					for(int ix=-B; ix<=gdata.mx[0]+B; ++ix) {
						bcVals_User_ye[q](ix, iy, iz) = 
							gdata.om_user[q](ix, iy, iz);
					}
				}
			}
		}
#endif
	}


	if(bc_Type[4] == 7) { // z, bottom

		int lbound[3] = {-B,-B,-B};
		int ubound[3] = {gdata.mx[0]+B,gdata.mx[1]+B,-1};

		for(int q=0; q<n_omInt; ++q) {
			bcVals_zb[q].resize(lbound, ubound);

			// Setting field at lower z bound:
			for(int iz=-B; iz<=-1; ++iz) {
				for(int iy=-B; iy<=gdata.mx[1]+B; ++iy) {
					for(int ix=-B; ix<=gdata.mx[0]+B; ++ix) {
						bcVals_zb[q](ix, iy, iz) = gdata.om[q](ix, iy, iz);
					}
				}
			}
		}
#if (OMS_USER == TRUE)
		for(int q=0; q<n_omIntUser; ++q) {
			bcVals_User_zb[q].resize(lbound, ubound);
			
			// Setting user field at lower z bound:
			for(int iz=-B; iz<=-1; ++iz) {
				for(int iy=-B; iy<=gdata.mx[1]+B; ++iy) {
					for(int ix=-B; ix<=gdata.mx[0]+B; ++ix) {
						bcVals_User_zb[q](ix, iy, iz) = 
							gdata.om_user[q](ix, iy, iz);
					}
				}
			}
		}
#endif
	}

	if(bc_Type[5] == 7) { // z, top

		int lbound[3] = {-B,-B,gdata.mx[2]+1};
		int ubound[3] = {gdata.mx[0]+B,gdata.mx[1]+B,gdata.mx[2]+B};

		for(int q=0; q<n_omInt; ++q) {
			bcVals_ze[q].resize(lbound, ubound);

			// Setting field at upper z bound:
			for(int iz=gdata.mx[2]+1; iz<=gdata.mx[2]+B; ++iz) {
				for(int iy=-B; iy<=gdata.mx[1]+B; ++iy) {
					for(int ix=-B; ix<=gdata.mx[0]+B; ++ix) {
						bcVals_ze[q](ix, iy, iz) = gdata.om[q](ix, iy, iz);
					}
				}
			}
		}
#if (OMS_USER == TRUE)
		for(int q=0; q<n_omIntUser; ++q) {
			bcVals_User_ze[q].resize(lbound, ubound);

			// Setting user field at upper z bound:
			for(int iz=gdata.mx[2]+1; iz<=gdata.mx[2]+B; ++iz) {
				for(int iy=-B; iy<=gdata.mx[1]+B; ++iy) {
					for(int ix=-B; ix<=gdata.mx[0]+B; ++ix) {
						bcVals_User_ze[q](ix, iy, iz) = 
							gdata.om_user[q](ix, iy, iz);
					}
				}
			}
		}
#endif
	}
	

}



void gridFunc::boundary(Data &gdata, ProblemType &pr)
{
	/*
	  When doing the extrapolation:
	  First have to extrapolate velocity in order to use the latter for
	  other extrapolations
	*/

#if(FLUID_TYPE==CRONOS_MULTIFLUID)
	// Loop over all fluids:
	for(int iFluid=0; iFluid<gdata.fluids->get_numFluids(); ++iFluid) {
		CronosFluid fluid = gdata.fluids->fluids[iFluid];
#else
		CronosFluid fluid = gdata.fluid;
#endif
		int q_rho = fluid.get_q_rho_global();
		int q_sx = fluid.get_q_sx_global();
		int q_sy = fluid.get_q_sy_global();
		int q_sz = fluid.get_q_sz_global();


		boundary(gdata, pr, gdata.om[q_sx],B,q_sx);
		boundary(gdata, pr, gdata.om[q_sy],B,q_sy);
		boundary(gdata, pr, gdata.om[q_sz],B,q_sz);
		boundary(gdata, pr, gdata.om[q_rho],B,q_rho);

		//   for (int q = 4; q < N_OMEGA; q++) {
		//     boundary(gdata, Problem, gdata.om[q],B,q);
		//   }

		for (int q = 4; q < n_omInt; q++) {
			boundary(gdata, pr, gdata.om[q],B,q);
		}

#if (OMS_USER == TRUE)
		for (int q=0; q<n_omIntUser; q++) {
			boundary(gdata, pr, gdata.om_user[q], B, q);
		}
#endif

#if(FLUID_TYPE==CRONOS_MULTIFLUID)
	}
#endif
}


void gridFunc::boundary(Data &gdata, ProblemType &Problem,
                        NumMatrix<double,3> &omb, int rim)
{
	boundary(gdata, Problem, omb, rim,-1);
};



void gridFunc::boundary_old(Data &gdata, ProblemType &Problem,
                        NumMatrix<double,3> &omb, int rim, int q, int iFluid)
{

	NumMatrix<double,3> Recv;
	NumMatrix<double,3> Send;

	/*
	  Determine type of boundary condition for the individual directions:
	  (-1) -> MPI-periodic (only for internal use)
	  (0) -> Periodic
	  (1) -> Periodic (with data at mx to be overwritten)
	  (2) -> Extrapolation
	  (3) -> Outflow
	  (4) -> User defined local bcs (bc_User)
	  (5) -> Reflecting boundaries
	  (6) -> singularity axis boundary conditions
	  (7) -> constant boundaries (kept from initial conditions) - still under construction

	*/

	//-----------------------------------------------------------
	// First we determine which field we are actually dealing with
	// -----------------------------------------------------------

	if(q>-1) {
#if(FLUID_TYPE==CRONOS_MULTIFLUID)
		string fluidName = gdata.fluids->fluids[iFluid].get_Name();
		q = Problem.get_q(omb, fluidName);
#else
		q = Problem.get_q(omb);
#endif
	}


	//-----------------------------------------------------------
	// x-direction
	//-----------------------------------------------------------

	// Periodic AND MPI boundary conditions 
	bc_Periodic(gdata, Problem, omb, 0, q, rim);

	// Left:
	if(bc_Type[0] > 1) {
		if (bc_Type[0] == 2) {
			bc_Extrapolate(gdata, Problem, omb,0,0,q,rim);
		} else if (bc_Type[0] == 3) {
			bc_Outflow(gdata, Problem, omb,0,0,q,rim);
		} else if (bc_Type[0] == 4) {
			Problem.bc_User(gdata, omb,0,0,q,rim);
		} else if (bc_Type[0] == 5) {
			bc_Reflecting(gdata, Problem, omb,0,0,q,rim);
		} else if (bc_Type[0] == 6) {
			bc_Axis(gdata, Problem, omb, 0, 0, q,rim);
		} else if (bc_Type[0] == 7) {
			if(omb.get_bcTypeLow(0) == 2) {
				bc_Extrapolate(gdata, Problem, omb,0,0,q,rim);
			} else {
				bc_Fixed(gdata, Problem, omb,0,0,q,rim);
			}
		} else {
			cerr << " No such boundary condition 1 " << endl;
		}
	}

	// Right:
	if(bc_Type[1] > 1) {
		if (bc_Type[1] == 2){
			bc_Extrapolate(gdata, Problem, omb,0,1,q,rim);
		} else if (bc_Type[1] == 3) {
			bc_Outflow(gdata, Problem, omb,0,1,q,rim);
		} else if (bc_Type[1] == 4) {
			Problem.bc_User(gdata, omb,0,1,q,rim);
		} else if (bc_Type[1] == 5) {
			bc_Reflecting(gdata, Problem, omb,0,1,q,rim);
		} else if (bc_Type[1] == 7) {
			if(omb.get_bcTypeHigh(0) == 2) {
				bc_Extrapolate(gdata, Problem, omb,0,1,q,rim);
			} else {
				bc_Fixed(gdata, Problem, omb,0,1,q,rim);
			}
		} else {
			cerr << " No such boundary condition 2 " << endl;
		}
	}

	//-----------------------------------------------------------
	// y-direction
	//-----------------------------------------------------------

	// Periodic AND MPI boundary conditions 
	bc_Periodic(gdata, Problem, omb, 1, q, rim);

	if(bc_Type[2] > 1) {
		if(bc_Type[2] == 2) {
			bc_Extrapolate(gdata, Problem, omb,1,0,q,rim);
		} else if(bc_Type[2] == 3) {
			bc_Outflow(gdata, Problem, omb,1,0,q,rim);
		} else if(bc_Type[2] == 4) {
			Problem.bc_User(gdata, omb,1,0,q,rim);
		} else if(bc_Type[2] == 5) {
			bc_Reflecting(gdata, Problem, omb,1,0,q,rim);
		} else if (bc_Type[2] == 6) {
			bc_Axis(gdata, Problem, omb, 1, 0, q,rim);
		} else if (bc_Type[2] == 7) {
			if(omb.get_bcTypeLow(1) == 2) {
				bc_Extrapolate(gdata, Problem, omb,1,0,q,rim);
			} else {
				bc_Fixed(gdata, Problem, omb,1,0,q,rim);
			}
		} else {
			cerr << " No such boundary condition 3 " << " " << bc_Type[1]<< endl;
		}
	}

	if(bc_Type[3] > 1) {
		if(bc_Type[3] == 2) {
			bc_Extrapolate(gdata, Problem, omb,1,1,q,rim);
		} else if(bc_Type[3] == 3) {
			bc_Outflow(gdata, Problem, omb,1,1,q,rim);
		} else if(bc_Type[3] == 4) {
			Problem.bc_User(gdata, omb,1,1,q,rim);
		} else if(bc_Type[3] == 5) {
			bc_Reflecting(gdata, Problem, omb,1,1,q,rim);
		} else if (bc_Type[3] == 6) {
			bc_Axis(gdata, Problem, omb, 1, 1, q,rim);
		} else if (bc_Type[3] == 7) {
			if(omb.get_bcTypeHigh(1) == 2) {
				bc_Extrapolate(gdata, Problem, omb,1,1,q,rim);
			} else {
				bc_Fixed(gdata, Problem, omb,1,1,q,rim);
			}
		} else {
			cerr << " No such boundary condition 4 " << endl;
		}
	}

	//-----------------------------------------------------------
	// z-direction
	//-----------------------------------------------------------

	// Periodic AND MPI boundary conditions 
	bc_Periodic(gdata, Problem, omb, 2, q, rim);

	if(bc_Type[4] > 1) {
		if (bc_Type[4] == 2) {
			bc_Extrapolate(gdata, Problem, omb,2,0,q,rim);
		} else if (bc_Type[4] == 3) {
			bc_Outflow(gdata, Problem, omb,2,0,q,rim);
		} else if (bc_Type[4] == 4) {
			Problem.bc_User(gdata, omb,2,0,q,rim);
		} else if (bc_Type[4] == 5) {
			bc_Reflecting(gdata, Problem, omb,2,0,q,rim);
		} else if (bc_Type[4] == 7) {
			if(omb.get_bcTypeLow(2) == 2) {
				bc_Extrapolate(gdata, Problem, omb,2,0,q,rim);
			} else {
				bc_Fixed(gdata, Problem, omb,2,0,q,rim);
			}
		} else {
			cerr << " No such boundary condition 5 " << endl;
		}
	}
  
	if(bc_Type[5] > 1) {
		if (bc_Type[5] == 2) {
			bc_Extrapolate(gdata, Problem, omb,2,1,q,rim);
		} else if (bc_Type[5] == 3) {
			bc_Outflow(gdata, Problem, omb,2,1,q,rim);
		} else if (bc_Type[5] == 4) {
			Problem.bc_User(gdata, omb,2,1,q,rim);
		} else if (bc_Type[5] == 5) {
			bc_Reflecting(gdata, Problem, omb,2,1,q,rim);
		} else if (bc_Type[5] == 7) {
			if(omb.get_bcTypeHigh(2) == 2) {
				bc_Extrapolate(gdata, Problem, omb,2,1,q,rim);
			} else {
				bc_Fixed(gdata, Problem, omb,2,1,q,rim);
			}
		} else {
			cerr << " No such boundary condition 6 " << endl;
		}
	}

	// Do correction for magnetic field on axis after all
	// magnetic boundaries have been applied
//	if(q == q_Bz) {
	if(omb.getName()=="B_z") {

#if (GEOM == CYLINDRICAL)
		if(gdata.get_singularity_treatment(0) > 0) {
			do_AxisValCorrectionCyl(gdata, Problem);
		}
#else
		// Correction at lower theta boundary
		if(gdata.get_singularity_treatment(2) > 0) {
			do_AxisValCorrectionSph(gdata, Problem, 0);
		}
		// Correction at lower theta boundary
		// if(gdata.get_singularity_treatment(2) > 0) {
		// 	do_AxisValCorrectionSph(gdata, Problem, 2);
		// }
		if(gdata.get_singularity_treatment(3) > 0) {
			do_AxisValCorrectionSph(gdata, Problem, 1);
		}
#endif
//#ifdef parallel
//		MPI_Barrier(gdata.comm3d);
//#endif
	}

}



void gridFunc::boundary(Data &gdata, ProblemType &Problem,
                        NumMatrix<double,3> &omb, int rim, int q, int iFluid)
{

	/*
	  Determine type of boundary condition for the individual directions:
	  (-1) -> MPI-periodic (only for internal use)
	  (0) -> Periodic
	  (1) -> Periodic (with data at mx to be overwritten)
	  (2) -> Extrapolation
	  (3) -> Outflow
	  (4) -> User defined local bcs (bc_User)
	  (5) -> Reflecting boundaries
	  (6) -> singularity axis boundary conditions
	  (7) -> constant boundaries (kept from initial conditions) - still under construction

	*/

	//-----------------------------------------------------------
	// First we determine which field we are actually dealing with
	// -----------------------------------------------------------

	if(q>-1) {
#if(FLUID_TYPE==CRONOS_MULTIFLUID)
		string fluidName = gdata.fluids->fluids[iFluid].get_Name();
		q = Problem.get_q(omb, fluidName);
#else
		q = Problem.get_q(omb);
#endif
	}


	// Before the boundary conditions at the outer boundary of the domain are
	// addressed we begin with the periodic boundary conditions - in particularl
	// those for the MPI-parallel case

//	// In parallel case start with sending all necessary data packages
//#ifdef parallel
//	bc_MPISend(gdata, Problem, omb, q, rim);
//#endif

	//-----------------------------------------------------------
	// x-direction
	//-----------------------------------------------------------


	// Non-MPI periodic boundary conditions
#ifndef parallel
	bc_Periodic_serial(gdata, Problem, omb, 0, q, rim);
#else
	bc_MPISend(gdata, Problem, omb, 0, q, rim);
#endif



	if(bc_Type[0] > 1) {
		if (bc_Type[0] == 2) {
			bc_Extrapolate(gdata, Problem, omb,0,0,q,rim);
		} else if (bc_Type[0] == 3) {
			bc_Outflow(gdata, Problem, omb,0,0,q,rim);
		} else if (bc_Type[0] == 4) {
			Problem.bc_User(gdata, omb,0,0,q,rim);
		} else if (bc_Type[0] == 5) {
			bc_Reflecting(gdata, Problem, omb,0,0,q,rim);
		} else if (bc_Type[0] == 6) {
			bc_Axis(gdata, Problem, omb, 0, 0, q,rim);
		} else if (bc_Type[0] == 7) {
			if(omb.get_bcTypeLow(0) == 2) {
				bc_Extrapolate(gdata, Problem, omb,0,0,q,rim);
			} else {
				bc_Fixed(gdata, Problem, omb,0,0,q,rim);
			}
		} else {
			cerr << " No such boundary condition 1 " << endl;
		}
	}

	// Right:
	if(bc_Type[1] > 1) {
		if (bc_Type[1] == 2){
			bc_Extrapolate(gdata, Problem, omb,0,1,q,rim);
		} else if (bc_Type[1] == 3) {
			bc_Outflow(gdata, Problem, omb,0,1,q,rim);
		} else if (bc_Type[1] == 4) {
			Problem.bc_User(gdata, omb,0,1,q,rim);
		} else if (bc_Type[1] == 5) {
			bc_Reflecting(gdata, Problem, omb,0,1,q,rim);
		} else if (bc_Type[1] == 7) {
			if(omb.get_bcTypeHigh(0) == 2) {
				bc_Extrapolate(gdata, Problem, omb,0,1,q,rim);
			} else {
				bc_Fixed(gdata, Problem, omb,0,1,q,rim);
			}
		} else {
			cerr << " No such boundary condition 2 " << endl;
		}
	}

#ifdef parallel
	bc_MPIRecv(gdata, Problem, omb, 0, q, rim);
#endif

	//-----------------------------------------------------------
	// y-direction
	//-----------------------------------------------------------

	// Non-MPI periodic boundary conditions
#ifndef parallel
	bc_Periodic_serial(gdata, Problem, omb, 1, q, rim);
#else
	bc_MPISend(gdata, Problem, omb, 1, q, rim);
#endif

	if(bc_Type[2] > 1) {
		if(bc_Type[2] == 2) {
			bc_Extrapolate(gdata, Problem, omb,1,0,q,rim);
		} else if(bc_Type[2] == 3) {
			bc_Outflow(gdata, Problem, omb,1,0,q,rim);
		} else if(bc_Type[2] == 4) {
			Problem.bc_User(gdata, omb,1,0,q,rim);
		} else if(bc_Type[2] == 5) {
			bc_Reflecting(gdata, Problem, omb,1,0,q,rim);
		} else if (bc_Type[2] == 6) {
			bc_Axis(gdata, Problem, omb, 1, 0, q,rim);
		} else if (bc_Type[2] == 7) {
			if(omb.get_bcTypeLow(1) == 2) {
				bc_Extrapolate(gdata, Problem, omb,1,0,q,rim);
			} else {
				bc_Fixed(gdata, Problem, omb,1,0,q,rim);
			}
		} else {
			cerr << " No such boundary condition 3 " << " " << bc_Type[1]<< endl;
		}
	}

	if(bc_Type[3] > 1) {
		if(bc_Type[3] == 2) {
			bc_Extrapolate(gdata, Problem, omb,1,1,q,rim);
		} else if(bc_Type[3] == 3) {
			bc_Outflow(gdata, Problem, omb,1,1,q,rim);
		} else if(bc_Type[3] == 4) {
			Problem.bc_User(gdata, omb,1,1,q,rim);
		} else if(bc_Type[3] == 5) {
			bc_Reflecting(gdata, Problem, omb,1,1,q,rim);
		} else if (bc_Type[3] == 6) {
			bc_Axis(gdata, Problem, omb, 1, 1, q,rim);
		} else if (bc_Type[3] == 7) {
			if(omb.get_bcTypeHigh(1) == 2) {
				bc_Extrapolate(gdata, Problem, omb,1,1,q,rim);
			} else {
				bc_Fixed(gdata, Problem, omb,1,1,q,rim);
			}
		} else {
			cerr << " No such boundary condition 4 " << endl;
		}
	}

#ifdef parallel
	bc_MPIRecv(gdata, Problem, omb, 1, q, rim);
#endif

	//-----------------------------------------------------------
	// z-direction
	//-----------------------------------------------------------

	// Non-MPI periodic boundary conditions
#ifndef parallel
	bc_Periodic_serial(gdata, Problem, omb, 2, q, rim);
#else
	bc_MPISend(gdata, Problem, omb, 2, q, rim);
#endif

	if(bc_Type[4] > 1) {
		if (bc_Type[4] == 2) {
			bc_Extrapolate(gdata, Problem, omb,2,0,q,rim);
		} else if (bc_Type[4] == 3) {
			bc_Outflow(gdata, Problem, omb,2,0,q,rim);
		} else if (bc_Type[4] == 4) {
			Problem.bc_User(gdata, omb,2,0,q,rim);
		} else if (bc_Type[4] == 5) {
			bc_Reflecting(gdata, Problem, omb,2,0,q,rim);
		} else if (bc_Type[4] == 7) {
			if(omb.get_bcTypeLow(2) == 2) {
				bc_Extrapolate(gdata, Problem, omb,2,0,q,rim);
			} else {
				bc_Fixed(gdata, Problem, omb,2,0,q,rim);
			}
		} else {
			cerr << " No such boundary condition 5 " << endl;
		}
	}

	if(bc_Type[5] > 1) {
		if (bc_Type[5] == 2) {
			bc_Extrapolate(gdata, Problem, omb,2,1,q,rim);
		} else if (bc_Type[5] == 3) {
			bc_Outflow(gdata, Problem, omb,2,1,q,rim);
		} else if (bc_Type[5] == 4) {
			Problem.bc_User(gdata, omb,2,1,q,rim);
		} else if (bc_Type[5] == 5) {
			bc_Reflecting(gdata, Problem, omb,2,1,q,rim);
		} else if (bc_Type[5] == 7) {
			if(omb.get_bcTypeHigh(2) == 2) {
				bc_Extrapolate(gdata, Problem, omb,2,1,q,rim);
			} else {
				bc_Fixed(gdata, Problem, omb,2,1,q,rim);
			}
		} else {
			cerr << " No such boundary condition 6 " << endl;
		}
	}

#ifdef parallel
	bc_MPIRecv(gdata, Problem, omb, 2, q, rim);
#endif
//	// Finally write all that that was received withing earlier MPI communications
//#ifdef parallel
//	bc_MPIRecv(gdata, Problem, omb, q, rim);
//#endif


	// Do correction for magnetic field on axis after all
	// magnetic boundaries have been applied
//	if(q == q_Bz) {
	if(omb.getName()=="B_z") {

#if (GEOM == CYLINDRICAL)
		if(gdata.get_singularity_treatment(0) > 0) {
			do_AxisValCorrectionCyl(gdata, Problem);
		}
#else
		// Correction at lower theta boundary
		if(gdata.get_singularity_treatment(2) > 0) {
			do_AxisValCorrectionSph(gdata, Problem, 0);
		}
		// Correction at lower theta boundary
		// if(gdata.get_singularity_treatment(2) > 0) {
		// 	do_AxisValCorrectionSph(gdata, Problem, 2);
		// }
		if(gdata.get_singularity_treatment(3) > 0) {
			do_AxisValCorrectionSph(gdata, Problem, 1);
		}
#endif
//#ifdef parallel
//		MPI_Barrier(gdata.comm3d);
//#endif
	}

}

void gridFunc::boundary_MPI(Data &gdata, ProblemType &Problem,
                        NumMatrix<double,3> &omb, int rim, int q, int iFluid)
{
	// Periodic AND MPI boundary conditions in x-direction
	bc_Periodic(gdata, Problem, omb, 0, q, rim, true);

	// Periodic AND MPI boundary conditions in y-direction
	bc_Periodic(gdata, Problem, omb, 1, q, rim, true);

	// Periodic AND MPI boundary conditions in z-direction
	bc_Periodic(gdata, Problem, omb, 2, q, rim, true);
}


#ifdef parallel
void gridFunc::bc_SendRecv(Data &gdata, NumMatrix<REAL,3> &Send,
                           NumMatrix<REAL,3> &Recv, int dir, int above,
                           int q, bool do_Send, bool do_Receive)
{
	// Determine size of data:
	int sizeSend = ((Send.getHigh(0)-Send.getLow(0) + 1)*
	                (Send.getHigh(1)-Send.getLow(1) + 1)*
	                (Send.getHigh(2)-Send.getLow(2) + 1));

	int sizeRecv = ((Recv.getHigh(0)-Recv.getLow(0) + 1)*
	                (Recv.getHigh(1)-Recv.getLow(1) + 1)*
	                (Recv.getHigh(2)-Recv.getLow(2) + 1));


	int from, into;
	if(dir == 0) {
		if(above == 0) {
			from = gdata.right;
			into = gdata.left;
		} else {
			from = gdata.left;
			into = gdata.right;
		}
	} else if (dir == 1) {
		if(above == 0) {
			from = gdata.back;
			into = gdata.front;
		} else {
			from = gdata.front;
			into = gdata.back;
		}
	} else if (dir == 2) {
		if(above == 0) {
			from = gdata.top;
			into = gdata.bottom;
		} else {
			from = gdata.bottom;
			into = gdata.top;
		}
	} else {
		cerr << " No such dimension " << endl;
		exit(-33);
	}

	int numRequests = 0;
	if(do_Send) numRequests++;
	if(do_Receive) numRequests++;


	// Only do communication if necessary:
	if(from == gdata.rank && into == gdata.rank) {
		Recv = Send;
	} else {
		// Transfer of data:
		// initialize requests[4];
		//MPI_Request requests[numRequests];
		//MPI_Status statusrl[numRequests];
		std::vector<MPI_Request> requests(numRequests);
		std::vector<MPI_Status> statusrl(numRequests);
		for(int ireq=0; ireq<numRequests; ireq++) {
			requests[ireq] = MPI_REQUEST_NULL;
		}
    
		// receive data
		// message tag -- must not be less than 0!
		int ireq(0);
		if(do_Receive) {
			int tag = q + n_omInt*dir + (n_omInt+3)*from + 2;

			MPI_Irecv((double *)Recv, sizeRecv, MPI_DOUBLE, from , tag,
			          gdata.comm3d, &requests[ireq]);
			ireq++;
		}

		if(do_Send) {
			int tag = q + n_omInt*dir + (n_omInt+3)*gdata.rank + 2;
    
			MPI_Isend((double *)Send, sizeSend, MPI_DOUBLE, into, tag,
			          gdata.comm3d, &requests[ireq]);
		}


		/* wait for all communication to complete */
    
		if(numRequests > 0) {
			MPI_Waitall(numRequests, requests.data(), statusrl.data());
		}
    
	}
}

#else
void gridFunc::bc_SendRecv(Data &gdata, NumMatrix<REAL,3> &Send,
                           NumMatrix<REAL,3> &Recv, int dir, int above,
                           int q, bool do_Send, bool do_Receive)
{
	// (Nearly) Nothing to be done
	Recv = Send;
}
#endif




#ifdef parallel

void gridFunc::bc_SendRecv(Data &gdata, NumArray<REAL> &Send,
                           NumArray<REAL> &Recv, int dir, int size_send, int size_recv, int above,
                           int q, bool do_Send, bool do_Receive)
{
	int from, into;
	if(dir == 0) {
		if(above == 0) {
			from = gdata.right;
			into = gdata.left;
		} else {
			from = gdata.left;
			into = gdata.right;
		}
	} else if (dir == 1) {
		if(above == 0) {
			from = gdata.back;
			into = gdata.front;
		} else {
			from = gdata.front;
			into = gdata.back;
		}
	} else if (dir == 2) {
		if(above == 0) {
			from = gdata.top;
			into = gdata.bottom;
		} else {
			from = gdata.bottom;
			into = gdata.top;
		}
	} else {
		cerr << " No such dimension " << endl;
		exit(-33);
	}

	int numRequests = 0;
	if(do_Send) numRequests++;
	if(do_Receive) numRequests++;


	// Only do communication if necessary:
	if(from == gdata.rank && into == gdata.rank) {
		Recv = Send;
	} else {
		// Transfer of data:
		// initialize requests[4];
		//MPI_Request requests[numRequests];
		//MPI_Status statusrl[numRequests];
		std::vector<MPI_Request> requests(numRequests);
		std::vector<MPI_Status> statusrl(numRequests);
		for(int ireq=0; ireq<numRequests; ireq++) {
			requests[ireq] = MPI_REQUEST_NULL;
		}

		// receive data
		// message tag -- must not be less than 0!
		int ireq(0);
		if(do_Receive) {
			int tag = q + n_omInt*dir + (n_omInt+3)*from + 2;

			MPI_Irecv((double *)Recv, size_recv, MPI_DOUBLE, from , tag,
			          gdata.comm3d, &requests[ireq]);
			ireq++;
		}

		if(do_Send) {
			int tag = q + n_omInt*dir + (n_omInt+3)*gdata.rank + 2;

			MPI_Isend((double *)Send, size_send, MPI_DOUBLE, into, tag,
			          gdata.comm3d, &requests[ireq]);
		}


		/* wait for all communication to complete */

		if(numRequests > 0) {
			MPI_Waitall(numRequests, requests.data(), statusrl.data());
		}

	}
}
#else
void gridFunc::bc_SendRecv(Data &gdata, NumArray<REAL> &Send,
		NumArray<REAL> &Recv, int dir, int size_send, int size_recv, int above, int q, bool do_Send, bool do_Receive)
{
	// (Nearly) Nothing to be done
	for(int iter=0; iter<size_send; iter++) {
		Recv(iter) = Send(iter);
	}
}
#endif


#ifdef parallel

void gridFunc::bc_ISendRecv(Data &gdata, NumArray<REAL> &Send, NumArray<REAL> &Recv,
		MPI_Request *requests, int dir, int size_send, int size_recv, int above,
		int q, bool do_Send, bool do_Receive)
{
	int from, into;
	if(dir == 0) {
		if(above == 0) {
			from = gdata.right;
			into = gdata.left;
		} else {
			from = gdata.left;
			into = gdata.right;
		}
	} else if (dir == 1) {
		if(above == 0) {
			from = gdata.back;
			into = gdata.front;
		} else {
			from = gdata.front;
			into = gdata.back;
		}
	} else if (dir == 2) {
		if(above == 0) {
			from = gdata.top;
			into = gdata.bottom;
		} else {
			from = gdata.bottom;
			into = gdata.top;
		}
	} else {
		cerr << " No such dimension " << endl;
		exit(-33);
	}

	int numRequests = 0;
	if(do_Send) numRequests++;
	if(do_Receive) numRequests++;


	// initialize requests[4];
	for(int ireq=0; ireq<numRequests; ireq++) {
		requests[ireq] = MPI_REQUEST_NULL;
	}

	// Only do communication if necessary:
	if(from == gdata.rank && into == gdata.rank) {
		Recv = Send;
	} else {
		// Transfer of data:

		// receive data
		// message tag -- must not be less than 0!
		int ireq(0);
		if(do_Receive) {
			int tag = q + n_omInt*dir + (n_omInt+3)*from + 2;

			MPI_Irecv((double *)Recv, size_recv, MPI_DOUBLE, from , tag,
			          gdata.comm3d, &requests[ireq]);
			ireq++;
		}

		if(do_Send) {
			int tag = q + n_omInt*dir + (n_omInt+3)*gdata.rank + 2;

			MPI_Isend((double *)Send, size_send, MPI_DOUBLE, into, tag,
			          gdata.comm3d, &requests[ireq]);
		}


		/* wait for all communication to complete */

	}
}
#else
void gridFunc::bc_ISendRecv(Data &gdata, NumArray<REAL> &Send, NumArray<REAL> &Recv,
		int dir, int size_send, int size_recv, int above, int q, bool do_Send, bool do_Receive)
{
	// (Nearly) Nothing to be done
	for(int iter=0; iter<size_send; iter++) {
		Recv(iter) = Send(iter);
	}
}
#endif





void gridFunc::bc_Periodic_old(Data &gdata, ProblemType &Problem,
                           NumMatrix<REAL,3> &omb,
                           int dir, int q, int rim, bool internal_only) {

	NumMatrix<double,3> Recv;
	NumMatrix<double,3> Send;

	bool SendLeft(false), RecvRight(false);

	if(bc_Type[2*dir] < 2) SendLeft = true;
	if(bc_Type[2*dir+1] < 2) RecvRight = true;
	
	if(internal_only) {
#ifdef parallel
		if(gdata.check_MpiLowest(dir)) {
			SendLeft = false;
		}
		if(gdata.check_MpiHighest(dir)) {
			RecvRight = false;
		}
#else
		SendLeft = false;
		RecvRight = false;
#endif
	}

	// First part: from left to right (e.g., u(mx+1) = u(1) / u(0))
	if (SendLeft) { // Sender is Periodic
		
		int iminSend(1);
		int imaxSend(rim);
		int shift(0);
		
		// For new gridding scheme size is always the same:
		if(gdata.get_EdgeGridding()) {
			// iminSend = 0;
			// imaxSend = rim-1;
			iminSend = 1;
			imaxSend = rim;
			shift = 1;
		} else {		
			iminSend = 1;
			imaxSend = rim;
			if(bc_Type[2*dir] == 1 || gdata.time < 0.01*gdata.dt) {
				iminSend = 0;
			}
		}
		if(dir == 0) {
			Send.resize(Index::set(iminSend, -rim, -rim),
			            Index::set(imaxSend, gdata.mx[1]+rim,
			                       gdata.mx[2]+rim));

			for (int k = -rim; k <= gdata.mx[2]+rim; k++){
				for (int j = -rim; j <= gdata.mx[1]+rim; j++){
					for (int i = iminSend; i <= imaxSend; i++){
						Send(i,j,k) = omb(i-shift,j,k);    // x_min -> x_max
					}
				}
			}

		} else if (dir == 1) {
			Send.resize(Index::set(-rim, iminSend, -rim),
			            Index::set(gdata.mx[0]+rim, imaxSend,
			                       gdata.mx[2]+rim));

			for (int j = iminSend; j <= imaxSend; j++){
				for (int k = -rim; k <= gdata.mx[2]+rim; k++){
					for (int i = -rim; i <= gdata.mx[0]+rim; i++){
						Send(i,j,k) = omb(i,j-shift,k);    // y_min -> y_max
					}
				}
			}

		} else if (dir == 2) {
			Send.resize(Index::set(-rim, -rim, iminSend),
			            Index::set(gdata.mx[0]+rim, gdata.mx[1]+rim,
			                       imaxSend));
				
			for (int k = iminSend; k <= imaxSend; k++){
				for (int j = -rim; j <= gdata.mx[1]+rim; j++){
					for (int i = -rim; i <= gdata.mx[0]+rim; i++){
						Send(i,j,k) = omb(i,j,k-shift);    // z_max side
					}
				}
			}
				
		}
			
		if(bc_Type[2*dir] > -1) { // User specific periodic bcs
			Problem.bc_periodic(gdata, Send, omb, q, dir, false);
		}
	} else {
		Send.resize(Index::set(0,0,0), Index::set(0,0,0));
	}

		
	int iminRecv(1);
	int imaxRecv(rim);
	if(RecvRight) {  // Receiver is Periodic
		if(gdata.get_EdgeGridding()) {
			// iminRecv = 1; // Corresponding to mx[0]+1
			// iminRecv = 0; // Corresponding to mx[0]+1
			// imaxRecv = rim-1;
			iminRecv = 1; // Corresponding to mx[0]+1
			imaxRecv = rim;
		} else {
			if(bc_Type[2*dir+1] == 1 || gdata.time < 0.01*gdata.dt) {
				iminRecv = 0;
			}
		}
		if(dir == 0) {
			Recv.resize(Index::set(iminRecv,-rim,-rim),
			            Index::set(imaxRecv,gdata.mx[1]+rim,gdata.mx[2]+rim));
		} else if (dir == 1) {
			Recv.resize(Index::set(-rim,iminRecv,-rim),
			            Index::set(gdata.mx[0]+rim,imaxRecv,gdata.mx[2]+rim));
		} else if (dir == 2) {
			Recv.resize(Index::set(-rim,-rim,iminRecv),
			            Index::set(gdata.mx[0]+rim,gdata.mx[1]+rim,imaxRecv));
		}
	} else {
		Recv.resize(Index::set(0,0,0), Index::set(0,0,0));
	}
    

	if(SendLeft || RecvRight) { // Periodic
		bc_SendRecv(gdata, Send, Recv, dir, 0, q, SendLeft, RecvRight);
	}



		
	if(RecvRight) {
		if(dir == 0) {
			for (int k = -rim; k <= gdata.mx[2]+rim; k++){
				for (int j = -rim; j <= gdata.mx[1]+rim; j++){
					for (int i = iminRecv; i <= imaxRecv; i++){
						omb(gdata.mx[0]+i,j,k) = Recv(i,j,k);
					}
				}
			}
		} else if (dir == 1) {
			for (int j = iminRecv; j <= imaxRecv; j++){
				for (int k = -rim; k <= gdata.mx[2]+rim; k++){
					for (int i = -rim; i <= gdata.mx[0]+rim; i++){
						// y_min -> y_max
						omb(i,gdata.mx[1]+j,k) = Recv(i,j,k); 
					}
				}
			}
		} else if (dir == 2) {
			for (int k = iminRecv; k <= imaxRecv; k++){
				for (int j = -rim; j <= gdata.mx[1]+rim; j++){
					for (int i = -rim; i <= gdata.mx[0]+rim; i++){
						// z_min -> z_max
						omb(i,j,gdata.mx[2]+k) = Recv(i,j,k);
					}
				}
			}
		}
	}


	bool SendRight(false), RecvLeft(false);

	if(bc_Type[2*dir+1] < 2) SendRight = true;
	if(bc_Type[2*dir]   < 2) RecvLeft  = true;

	if(internal_only) {
#ifdef parallel
		if(gdata.check_MpiLowest(dir)) {
			RecvLeft = false;
		}
		if(gdata.check_MpiHighest(dir)) {
			SendRight = false;
		}
#else
		SendLeft = false;
		RecvRight = false;
#endif
	}

	// Second part: from right to left (e.g., u(-1) = u(mx-1) / u(mx))
	if(SendRight) { // Periodic

		int iminSend(-rim); // Corresponding to mx[0]+i
		int imaxSend(-1);
		int shift(0);

		// For new gridding scheme size is always the same:
		if(gdata.get_EdgeGridding()) {
			// iminSend = -rim+1;
			// imaxSend =  0;
			iminSend = -rim;
			imaxSend =  -1;
			if((bc_Type[2*dir+1] > -1 || gdata.time < 0.01*gdata.dt) &&
			   //if((bc_Type[2*dir+1] > -2 || gdata.time < 0.01*gdata.dt) &&
			   ((dir==0 && q==q_Bx) || (dir==1 && q==q_By) || 
			    (dir==2 && q==q_Bz))) {
				// imaxSend = -1;
				imaxSend = -2;
			}
			// shift=-1;
			shift=1;
		} else {		
			iminSend = -rim;
			imaxSend = -1;
			if((bc_Type[2*dir+1] > -1 || gdata.time < 0.01*gdata.dt) && 
			   ((dir==0 && q==q_Bx) || (dir==1 && q==q_By) ||
			    (dir==2 && q==q_Bz))) {
				imaxSend = -2;
			}
		}
		if(dir == 0) {
			Send.resize(Index::set(iminSend,-rim,-rim),
			            Index::set(imaxSend,gdata.mx[1]+rim,
			                       gdata.mx[2]+rim));

			for (int k = -rim; k <= gdata.mx[2]+rim; k++){
				for (int j = -rim; j <= gdata.mx[1]+rim; j++){
					for (int i = imaxSend; i >= iminSend; --i){
						// x_max -> x_min
						Send(i,j,k) = omb(gdata.mx[0]+i+shift,j,k); 
					}
				}
			}

		} else if (dir == 1) {
			Send.resize(Index::set(-rim,iminSend,-rim),
			            Index::set(gdata.mx[0]+rim,imaxSend,
			                       gdata.mx[2]+rim));

			for (int j = imaxSend; j >= iminSend; --j){
				for (int k = -rim; k <= gdata.mx[2]+rim; k++){
					for (int i = -rim; i <= gdata.mx[0]+rim; i++){
						// y_max -> y_min
						Send(i,j,k) = omb(i,gdata.mx[1]+j+shift,k);
					}
				}
			}

		} else if (dir == 2) {
			Send.resize(Index::set(-rim,-rim, iminSend),
			            Index::set(gdata.mx[0]+rim,gdata.mx[1]+rim,
			                       imaxSend));

			for (int k = imaxSend; k >= iminSend; --k){
				for (int j = -rim; j <= gdata.mx[1]+rim; j++){
					for (int i = -rim; i <= gdata.mx[0]+rim; i++){
						// z_max -> z_min
						Send(i,j,k) = omb(i,j,gdata.mx[2]+k+shift); 
					}
				}
			}

		}
    
		if(bc_Type[2*dir+1] > -1) {
			Problem.bc_periodic(gdata, Send, omb, q, dir, true);
		}
	} else {
		Send.resize(Index::set(0,0,0), Index::set(0,0,0));
	}
	
	imaxRecv = -1;
	if(RecvLeft) {
		if(gdata.get_EdgeGridding()) {
			//if((bc_Type[2*dir] > -2 || gdata.time < 0.01*gdata.dt) && 
			if((bc_Type[2*dir] > -1 || gdata.time < 0.01*gdata.dt) && 
			   ((dir==0 && q==q_Bx) || (dir==1 && q==q_By) ||
			    (dir==2 && q==q_Bz))) {
				imaxRecv = -2; // Corresponding to mx[0]+1
			}
		} else {
			if((bc_Type[2*dir] > -1 || gdata.time < 0.01*gdata.dt) && 
			   ((dir==0 && q==q_Bx) || (dir==1 && q==q_By) ||
			    (dir==2 && q==q_Bz))) {
				imaxRecv = -2;
			}
		}
		if(dir == 0) {
			Recv.resize(Index::set(-rim,-rim,-rim),
			            Index::set(imaxRecv,gdata.mx[1]+rim,
			                       gdata.mx[2]+rim));
		} else if (dir == 1) {
			Recv.resize(Index::set(-rim,-rim,-rim),
			            Index::set(gdata.mx[0]+rim,imaxRecv,
			                       gdata.mx[2]+rim));
		} else if (dir == 2) {
			Recv.resize(Index::set(-rim,-rim,-rim),
			            Index::set(gdata.mx[0]+rim,gdata.mx[1]+rim,
			                       imaxRecv));
		}
	} else {
		Recv.resize(Index::set(0,0,0), Index::set(0,0,0));
	}
		
	if(SendRight || RecvLeft) {
		bc_SendRecv(gdata, Send, Recv, dir, 1, q, SendRight, RecvLeft);
	}
  

	if(RecvLeft) {
		if(dir == 0) {
			for (int k = -rim; k <= gdata.mx[2]+rim; k++){
				for (int j = -rim; j <= gdata.mx[1]+rim; j++){
					for (int i = imaxRecv; i >= -rim; --i){
						omb(i,j,k) = Recv(i,j,k);
					}
				}
			}
		} else if (dir == 1) {
			for (int j = imaxRecv; j >= -rim; --j){
				for (int k = -rim; k <= gdata.mx[2]+rim; k++){
					for (int i = -rim; i <= gdata.mx[0]+rim; i++){
						omb(i,j,k) = Recv(i,j,k);
					}
				}
			}
		} else if (dir == 2) {
			for (int k = imaxRecv; k >= -rim; --k){
				for (int j = -rim; j <= gdata.mx[1]+rim; j++){
					for (int i = -rim; i <= gdata.mx[0]+rim; i++){
						omb(i,j,k) = Recv(i,j,k);
					}
				}
			}
		}
	}

}

#ifndef parallel
void gridFunc::bc_Periodic_serial(Data &gdata, ProblemType &Problem,
		NumMatrix<REAL,3> &omb, int dir, int q, int rim) {

	int shift(0);
	if(gdata.get_EdgeGridding()) {
		shift = 1;
	}

	int i_min, i_max;

	// Check if direction is periodic
	if((bc_Type[2*dir] < 2) && (bc_Type[2*dir+1] < 2)) {

		// x-direction
		if(dir == 0) {
			// From left to right:
			get_bcBuffRange(gdata, i_min, i_max, 0, 1, q, rim);
			for (int k = -rim; k <= gdata.mx[2]+rim; k++){
				for (int j = -rim; j <= gdata.mx[1]+rim; j++){
					for (int i = i_min; i <= i_max; i++){
						omb(gdata.mx[0]+i,j,k) = omb(i-shift,j,k);
					}
				}
			}
			// From right to left:
			get_bcBuffRange(gdata, i_min, i_max, 0, 0, q, rim);
			for (int k = -rim; k <= gdata.mx[2]+rim; k++){
				for (int j = -rim; j <= gdata.mx[1]+rim; j++){
					for (int i = i_max; i >= i_min; --i){
						omb(i,j,k) = omb(gdata.mx[0]+i+shift,j,k);
					}
				}
			}

			// y-direction
		} else if (dir == 1) {
			// From left to right:
			get_bcBuffRange(gdata, i_min, i_max, 1, 1, q, rim);
			for (int j = i_min; j <= i_max; j++){
				for (int k = -rim; k <= gdata.mx[2]+rim; k++){
					for (int i = -rim; i <= gdata.mx[0]+rim; i++){
						// y_min -> y_max
						omb(i,gdata.mx[1]+j,k) = omb(i,j-shift,k);
					}
				}
			}
			// From right to left:
			get_bcBuffRange(gdata, i_min, i_max, 1, 0, q, rim);
			for (int j = i_max; j >= -rim; --j){
				for (int k = -rim; k <= gdata.mx[2]+rim; k++){
					for (int i = -rim; i <= gdata.mx[0]+rim; i++){
						omb(i,j,k) = omb(i,gdata.mx[1]+j+shift,k);
					}
				}
			}

			// z-direction
		} else {
			// From left to right:
			get_bcBuffRange(gdata, i_min, i_max, 2, 1, q, rim);
			for (int k = i_min; k <= i_max; k++){
				for (int j = -rim; j <= gdata.mx[1]+rim; j++){
					for (int i = -rim; i <= gdata.mx[0]+rim; i++){
						// z_min -> z_max
						omb(i,j,gdata.mx[2]+k) = omb(i,j,k-shift);
					}
				}
			}
			// From right to left:
			get_bcBuffRange(gdata, i_min, i_max, 2, 0, q, rim);
			for (int k = i_max; k >= -rim; --k){
				for (int j = -rim; j <= gdata.mx[1]+rim; j++){
					for (int i = -rim; i <= gdata.mx[0]+rim; i++){
						omb(i,j,k) = omb(i,j,gdata.mx[2]+k+shift);
					}
				}
			}
		}

	}
}
#endif

void gridFunc::bc_Periodic(Data &gdata, ProblemType &Problem,
                           NumMatrix<REAL,3> &omb,
                           int dir, int q, int rim, bool internal_only) {



	bool SendLeft(false), RecvRight(false);

	if(bc_Type[2*dir] < 2) SendLeft = true;
	if(bc_Type[2*dir+1] < 2) RecvRight = true;

	if(internal_only) {
#ifdef parallel
		if(gdata.check_MpiLowest(dir)) {
			SendLeft = false;
		}
		if(gdata.check_MpiHighest(dir)) {
			RecvRight = false;
		}
#else
		SendLeft = false;
		RecvRight = false;
#endif
	}

	int size_send = 0;
	// First part: from left to right (e.g., u(mx+1) = u(1) / u(0))
	if (SendLeft) { // Sender is Periodic

		int iminSend(1);
		int imaxSend(rim);
		int shift(0);

		// For new gridding scheme size is always the same:
		if(gdata.get_EdgeGridding()) {
			// iminSend = 0;
			// imaxSend = rim-1;
			iminSend = 1;
			imaxSend = rim;
			shift = 1;
		} else {
			iminSend = 1;
			imaxSend = rim;
			if(bc_Type[2*dir] == 1 || gdata.time < 0.01*gdata.dt) {
				iminSend = 0;
			}
		}
		if(dir == 0) {
			int iter=0;
			for (int k = -rim; k <= gdata.mx[2]+rim; k++){
				for (int j = -rim; j <= gdata.mx[1]+rim; j++){
					for (int i = iminSend; i <= imaxSend; i++){
						SendArr_buff(iter) = omb(i-shift,j,k);    // x_min -> x_max
						iter++;
					}
				}
			}
			size_send = iter;


		} else if (dir == 1) {
			int iter=0;
			for (int j = iminSend; j <= imaxSend; j++){
				for (int k = -rim; k <= gdata.mx[2]+rim; k++){
					for (int i = -rim; i <= gdata.mx[0]+rim; i++){
						SendArr_buff(iter) = omb(i,j-shift,k);    // y_min -> y_max
						iter++;
					}
				}
			}
			size_send = iter;

		} else if (dir == 2) {
			int iter=0;
			for (int k = iminSend; k <= imaxSend; k++){
				for (int j = -rim; j <= gdata.mx[1]+rim; j++){
					for (int i = -rim; i <= gdata.mx[0]+rim; i++){
						SendArr_buff(iter) = omb(i,j,k-shift);    // z_max side
						iter++;
					}
				}
			}
			size_send = iter;

		}

		if(bc_Type[2*dir] > -1) { // User specific periodic bcs
		  // COMMENT: not implemented for 1D storage
			// cerr << " Not implemented, yet" << endl;
			// exit(3);
//			Problem.bc_periodic(gdata, *Send, omb, q, dir, false);
		}

	} else {

		size_send = 0;
//		Send.resize(Index::set(0,0,0), Index::set(0,0,0));
	}


	int iminRecv(1);
	int imaxRecv(rim);
	int size_recv(0);
	if(RecvRight) {  // Receiver is Periodic
		if(gdata.get_EdgeGridding()) {
			// iminRecv = 1; // Corresponding to mx[0]+1
			// iminRecv = 0; // Corresponding to mx[0]+1
			// imaxRecv = rim-1;
			iminRecv = 1; // Corresponding to mx[0]+1
			imaxRecv = rim;
		} else {
			if(bc_Type[2*dir+1] == 1 || gdata.time < 0.01*gdata.dt) {
				iminRecv = 0;
			}
		}
		if(dir == 0) {
			size_recv = (imaxRecv-iminRecv+1)*(gdata.mx[1]+2*rim+1)*(gdata.mx[2]+2*rim+1);
		} else if (dir == 1) {
			size_recv = (gdata.mx[0]+2*rim+1)*(imaxRecv-iminRecv+1)*(gdata.mx[2]+2*rim+1);
		} else if (dir == 2) {
			size_recv = (gdata.mx[0]+2*rim+1)*(gdata.mx[1]+2*rim+1)*(imaxRecv-iminRecv+1);
		}
	} else {
		size_recv = 0;
	}



	if(SendLeft || RecvRight) { // Periodic
		bc_SendRecv(gdata, SendArr_buff, RecvArr_buff, dir, size_send, size_recv, 0, q, SendLeft, RecvRight);
	}




	if(RecvRight) {
		if(dir == 0) {
			int iter = 0;
			for (int k = -rim; k <= gdata.mx[2]+rim; k++){
				for (int j = -rim; j <= gdata.mx[1]+rim; j++){
					for (int i = iminRecv; i <= imaxRecv; i++){
						omb(gdata.mx[0]+i,j,k) = RecvArr_buff(iter);
						iter++;
					}
				}
			}
		} else if (dir == 1) {
			int iter = 0;
			for (int j = iminRecv; j <= imaxRecv; j++){
				for (int k = -rim; k <= gdata.mx[2]+rim; k++){
					for (int i = -rim; i <= gdata.mx[0]+rim; i++){
						// y_min -> y_max
						omb(i,gdata.mx[1]+j,k) = RecvArr_buff(iter);
						iter++;
					}
				}
			}
		} else if (dir == 2) {
			int iter = 0;
			for (int k = iminRecv; k <= imaxRecv; k++){
				for (int j = -rim; j <= gdata.mx[1]+rim; j++){
					for (int i = -rim; i <= gdata.mx[0]+rim; i++){
						// z_min -> z_max
						omb(i,j,gdata.mx[2]+k) = RecvArr_buff(iter);
						iter++;
					}
				}
			}
		}
	}


	bool SendRight(false), RecvLeft(false);

	if(bc_Type[2*dir+1] < 2) SendRight = true;
	if(bc_Type[2*dir]   < 2) RecvLeft  = true;

	if(internal_only) {
#ifdef parallel
		if(gdata.check_MpiLowest(dir)) {
			RecvLeft = false;
		}
		if(gdata.check_MpiHighest(dir)) {
			SendRight = false;
		}
#else
		SendLeft = false;
		RecvRight = false;
#endif
	}

	// Second part: from right to left (e.g., u(-1) = u(mx-1) / u(mx))
	if(SendRight) { // Periodic

		int iminSend(-rim); // Corresponding to mx[0]+i
		int imaxSend(-1);
		int shift(0);

		// For new gridding scheme size is always the same:
		if(gdata.get_EdgeGridding()) {
			// iminSend = -rim+1;
			// imaxSend =  0;
			iminSend = -rim;
			imaxSend =  -1;
			if((bc_Type[2*dir+1] > -1 || gdata.time < 0.01*gdata.dt) &&
			   //if((bc_Type[2*dir+1] > -2 || gdata.time < 0.01*gdata.dt) &&
			   ((dir==0 && q==q_Bx) || (dir==1 && q==q_By) ||
			    (dir==2 && q==q_Bz))) {
				// imaxSend = -1;
				imaxSend = -2;
			}
			// shift=-1;
			shift=1;
		} else {
			iminSend = -rim;
			imaxSend = -1;
			if((bc_Type[2*dir+1] > -1 || gdata.time < 0.01*gdata.dt) &&
			   ((dir==0 && q==q_Bx) || (dir==1 && q==q_By) ||
			    (dir==2 && q==q_Bz))) {
				imaxSend = -2;
			}
		}
		if(dir == 0) {
			int iter=0;
			for (int k = -rim; k <= gdata.mx[2]+rim; k++){
				for (int j = -rim; j <= gdata.mx[1]+rim; j++){
					for (int i = imaxSend; i >= iminSend; --i){
						// x_max -> x_min
						SendArr_buff(iter) = omb(gdata.mx[0]+i+shift,j,k);
						iter++;
					}
				}
			}
			size_send = iter;

		} else if (dir == 1) {
			int iter=0;
			for (int j = imaxSend; j >= iminSend; --j){
				for (int k = -rim; k <= gdata.mx[2]+rim; k++){
					for (int i = -rim; i <= gdata.mx[0]+rim; i++){
						// y_max -> y_min
						SendArr_buff(iter) = omb(i,gdata.mx[1]+j+shift,k);
						iter++;
					}
				}
			}
			size_send = iter;

		} else if (dir == 2) {
			int iter=0;
			for (int k = imaxSend; k >= iminSend; --k){
				for (int j = -rim; j <= gdata.mx[1]+rim; j++){
					for (int i = -rim; i <= gdata.mx[0]+rim; i++){
						// z_max -> z_min
						SendArr_buff(iter) = omb(i,j,gdata.mx[2]+k+shift);
						iter++;
					}
				}
			}
			size_send = iter;

		}

//		if(bc_Type[2*dir+1] > -1) {
//			Problem.bc_periodic(gdata, Send, omb, q, dir, true);
//		}
	} else {
		size_send = 0;
	}

	imaxRecv = -1;
	if(RecvLeft) {
		if(gdata.get_EdgeGridding()) {
			//if((bc_Type[2*dir] > -2 || gdata.time < 0.01*gdata.dt) &&
			if((bc_Type[2*dir] > -1 || gdata.time < 0.01*gdata.dt) &&
					((dir==0 && q==q_Bx) || (dir==1 && q==q_By) ||
							(dir==2 && q==q_Bz))) {
				imaxRecv = -2; // Corresponding to mx[0]+1
			}
		} else {
			if((bc_Type[2*dir] > -1 || gdata.time < 0.01*gdata.dt) &&
					((dir==0 && q==q_Bx) || (dir==1 && q==q_By) ||
							(dir==2 && q==q_Bz))) {
				imaxRecv = -2;
			}
		}
		if(dir == 0) {
			size_recv = (imaxRecv+rim+1)*(gdata.mx[1]+2*rim+1)*(gdata.mx[2]+2*rim+1);
		} else if (dir == 1) {
			size_recv = (gdata.mx[0]+2*rim+1)*(imaxRecv+rim+1)*(gdata.mx[2]+2*rim+1);
		} else if (dir == 2) {
			size_recv = (gdata.mx[0]+2*rim+1)*(gdata.mx[1]+2*rim+1)*(imaxRecv+rim+1);
		}
	} else {
		size_recv = 0;
	}


	if(SendRight || RecvLeft) {
		bc_SendRecv(gdata, SendArr_buff, RecvArr_buff, dir, size_send, size_recv, 1, q, SendRight, RecvLeft);
	}


	if(RecvLeft) {
		if(dir == 0) {
			int iter = 0;
			for (int k = -rim; k <= gdata.mx[2]+rim; k++){
				for (int j = -rim; j <= gdata.mx[1]+rim; j++){
					for (int i = imaxRecv; i >= -rim; --i){
						omb(i,j,k) = RecvArr_buff(iter);
						iter++;
					}
				}
			}
		} else if (dir == 1) {
			int iter = 0;
			for (int j = imaxRecv; j >= -rim; --j){
				for (int k = -rim; k <= gdata.mx[2]+rim; k++){
					for (int i = -rim; i <= gdata.mx[0]+rim; i++){
						omb(i,j,k) = RecvArr_buff(iter);
						iter++;
					}
				}
			}
		} else if (dir == 2) {
			int iter = 0;
			for (int k = imaxRecv; k >= -rim; --k){
				for (int j = -rim; j <= gdata.mx[1]+rim; j++){
					for (int i = -rim; i <= gdata.mx[0]+rim; i++){
						omb(i,j,k) = RecvArr_buff(iter);;
						iter++;
					}
				}
			}
		}
	}

}


#ifdef parallel
void gridFunc::bc_MPISend(Data &gdata, ProblemType &Problem, NumMatrix<REAL,3> &omb,
		int dir, int q, int rim) {



	bool SendLeft(false), RecvRight(false);

	if(bc_Type[2*dir] < 2) SendLeft = true;
	if(bc_Type[2*dir+1] < 2) RecvRight = true;

#ifndef parallel
	SendLeft = false;
	RecvRight = false;
#endif

	int size_send_LR = 0;
	// First part: from left to right (e.g., u(mx+1) = u(1) / u(0))
	if (SendLeft) { // Sender is Periodic

		int iminSend(1);
		int imaxSend(rim);
		int shift(0);

		// For new gridding scheme size is always the same:
		if(gdata.get_EdgeGridding()) {
			// iminSend = 0;
			// imaxSend = rim-1;
			iminSend = 1;
			imaxSend = rim;
			shift = 1;
		} else {
			iminSend = 1;
			imaxSend = rim;
			if(bc_Type[2*dir] == 1 || gdata.time < 0.01*gdata.dt) {
				iminSend = 0;
			}
		}
		if(dir == 0) {
			int iter=0;
			for (int k = -rim; k <= gdata.mx[2]+rim; k++){
				for (int j = -rim; j <= gdata.mx[1]+rim; j++){
					for (int i = iminSend; i <= imaxSend; i++){
						SendArrLR(iter) = omb(i-shift,j,k);    // x_min -> x_max
						iter++;
					}
				}
			}
			size_send_LR = iter;


		} else if (dir == 1) {
			int iter=0;
			for (int j = iminSend; j <= imaxSend; j++){
				for (int k = -rim; k <= gdata.mx[2]+rim; k++){
					for (int i = -rim; i <= gdata.mx[0]+rim; i++){
						SendArrLR(iter) = omb(i,j-shift,k);    // y_min -> y_max
						iter++;
					}
				}
			}
			size_send_LR = iter;

		} else if (dir == 2) {
			int iter=0;
			for (int k = iminSend; k <= imaxSend; k++){
				for (int j = -rim; j <= gdata.mx[1]+rim; j++){
					for (int i = -rim; i <= gdata.mx[0]+rim; i++){
						SendArrLR(iter) = omb(i,j,k-shift);    // z_max side
						iter++;
					}
				}
			}
			size_send_LR = iter;

		}

		if(bc_Type[2*dir] > -1) { // User specific periodic bcs
		  // COMMENT: not implemented for 1D storage
		  //			cerr << " Not implemented, yet" << endl;
		  //	exit(3);
//			Problem.bc_periodic(gdata, *Send, omb, q, dir, false);
		}

	} else {

		size_send_LR = 0;
//		Send.resize(Index::set(0,0,0), Index::set(0,0,0));
	}


	int iminRecv(1);
	int imaxRecv(rim);
	int size_recv_LR(0);
	if(RecvRight) {  // Receiver is Periodic
		if(gdata.get_EdgeGridding()) {
			// iminRecv = 1; // Corresponding to mx[0]+1
			// iminRecv = 0; // Corresponding to mx[0]+1
			// imaxRecv = rim-1;
			iminRecv = 1; // Corresponding to mx[0]+1
			imaxRecv = rim;
		} else {
			if(bc_Type[2*dir+1] == 1 || gdata.time < 0.01*gdata.dt) {
				iminRecv = 0;
			}
		}
		if(dir == 0) {
			size_recv_LR = (imaxRecv-iminRecv+1)*(gdata.mx[1]+2*rim+1)*(gdata.mx[2]+2*rim+1);
		} else if (dir == 1) {
			size_recv_LR = (gdata.mx[0]+2*rim+1)*(imaxRecv-iminRecv+1)*(gdata.mx[2]+2*rim+1);
		} else if (dir == 2) {
			size_recv_LR = (gdata.mx[0]+2*rim+1)*(gdata.mx[1]+2*rim+1)*(imaxRecv-iminRecv+1);
		}
	} else {
		size_recv_LR = 0;
	}



	if(SendLeft || RecvRight) { // Periodic
#ifdef parallel
		bc_ISendRecv(gdata, SendArrLR, RecvArrLR, requests_LR, dir,
				size_send_LR, size_recv_LR, 0, q, SendLeft, RecvRight);
#else
//		bc_ISendRecv(gdata, SendArrLR, RecvArrLR, dir, requestsLR,
//				size_send_LR, size_recv_LR, 0, q, SendLeft, RecvRight);
#endif
	}




//	if(RecvRight) {
//		if(dir == 0) {
//			int iter = 0;
//			for (int k = -rim; k <= gdata.mx[2]+rim; k++){
//				for (int j = -rim; j <= gdata.mx[1]+rim; j++){
//					for (int i = iminRecv; i <= imaxRecv; i++){
//						omb(gdata.mx[0]+i,j,k) = RecvArr_buff(iter);
//						iter++;
//					}
//				}
//			}
//		} else if (dir == 1) {
//			int iter = 0;
//			for (int j = iminRecv; j <= imaxRecv; j++){
//				for (int k = -rim; k <= gdata.mx[2]+rim; k++){
//					for (int i = -rim; i <= gdata.mx[0]+rim; i++){
//						// y_min -> y_max
//						omb(i,gdata.mx[1]+j,k) = RecvArr_buff(iter);
//						iter++;
//					}
//				}
//			}
//		} else if (dir == 2) {
//			int iter = 0;
//			for (int k = iminRecv; k <= imaxRecv; k++){
//				for (int j = -rim; j <= gdata.mx[1]+rim; j++){
//					for (int i = -rim; i <= gdata.mx[0]+rim; i++){
//						// z_min -> z_max
//						omb(i,j,gdata.mx[2]+k) = RecvArr_buff(iter);
//						iter++;
//					}
//				}
//			}
//		}
//	}


	bool SendRight(false), RecvLeft(false);

	if(bc_Type[2*dir+1] < 2) SendRight = true;
	if(bc_Type[2*dir]   < 2) RecvLeft  = true;

#ifndef parallel
	SendLeft = false;
	RecvRight = false;
#endif

	int size_send_RL(0);
	// Second part: from right to left (e.g., u(-1) = u(mx-1) / u(mx))
	if(SendRight) { // Periodic

		int iminSend(-rim); // Corresponding to mx[0]+i
		int imaxSend(-1);
		int shift(0);

		// For new gridding scheme size is always the same:
		if(gdata.get_EdgeGridding()) {
			// iminSend = -rim+1;
			// imaxSend =  0;
			iminSend = -rim;
			imaxSend =  -1;
			if((bc_Type[2*dir+1] > -1 || gdata.time < 0.01*gdata.dt) &&
			   //if((bc_Type[2*dir+1] > -2 || gdata.time < 0.01*gdata.dt) &&
			   ((dir==0 && q==q_Bx) || (dir==1 && q==q_By) ||
			    (dir==2 && q==q_Bz))) {
				// imaxSend = -1;
				imaxSend = -2;
			}
			// shift=-1;
			shift=1;
		} else {
			iminSend = -rim;
			imaxSend = -1;
			if((bc_Type[2*dir+1] > -1 || gdata.time < 0.01*gdata.dt) &&
			   ((dir==0 && q==q_Bx) || (dir==1 && q==q_By) ||
			    (dir==2 && q==q_Bz))) {
				imaxSend = -2;
			}
		}
		if(dir == 0) {
			int iter=0;
			for (int k = -rim; k <= gdata.mx[2]+rim; k++){
				for (int j = -rim; j <= gdata.mx[1]+rim; j++){
					for (int i = imaxSend; i >= iminSend; --i){
						// x_max -> x_min
						SendArrRL(iter) = omb(gdata.mx[0]+i+shift,j,k);
						iter++;
					}
				}
			}
			size_send_RL = iter;

		} else if (dir == 1) {
			int iter=0;
			for (int j = imaxSend; j >= iminSend; --j){
				for (int k = -rim; k <= gdata.mx[2]+rim; k++){
					for (int i = -rim; i <= gdata.mx[0]+rim; i++){
						// y_max -> y_min
						SendArrRL(iter) = omb(i,gdata.mx[1]+j+shift,k);
						iter++;
					}
				}
			}
			size_send_RL = iter;

		} else if (dir == 2) {
			int iter=0;
			for (int k = imaxSend; k >= iminSend; --k){
				for (int j = -rim; j <= gdata.mx[1]+rim; j++){
					for (int i = -rim; i <= gdata.mx[0]+rim; i++){
						// z_max -> z_min
						SendArrRL(iter) = omb(i,j,gdata.mx[2]+k+shift);
						iter++;
					}
				}
			}
			size_send_RL = iter;

		}

//		if(bc_Type[2*dir+1] > -1) {
//			Problem.bc_periodic(gdata, Send, omb, q, dir, true);
//		}
	} else {
		size_send_RL = 0;
	}

	int size_recv_RL(0);
	imaxRecv = -1;
	if(RecvLeft) {
		if(gdata.get_EdgeGridding()) {
			//if((bc_Type[2*dir] > -2 || gdata.time < 0.01*gdata.dt) &&
			if((bc_Type[2*dir] > -1 || gdata.time < 0.01*gdata.dt) &&
					((dir==0 && q==q_Bx) || (dir==1 && q==q_By) ||
							(dir==2 && q==q_Bz))) {
				imaxRecv = -2; // Corresponding to mx[0]+1
			}
		} else {
			if((bc_Type[2*dir] > -1 || gdata.time < 0.01*gdata.dt) &&
					((dir==0 && q==q_Bx) || (dir==1 && q==q_By) ||
							(dir==2 && q==q_Bz))) {
				imaxRecv = -2;
			}
		}
		if(dir == 0) {
			size_recv_RL = (imaxRecv+rim+1)*(gdata.mx[1]+2*rim+1)*(gdata.mx[2]+2*rim+1);
		} else if (dir == 1) {
			size_recv_RL = (gdata.mx[0]+2*rim+1)*(imaxRecv+rim+1)*(gdata.mx[2]+2*rim+1);
		} else if (dir == 2) {
			size_recv_RL = (gdata.mx[0]+2*rim+1)*(gdata.mx[1]+2*rim+1)*(imaxRecv+rim+1);
		}
	} else {
		size_recv_RL = 0;
	}


	if(SendRight || RecvLeft) {
#ifdef parallel
		bc_ISendRecv(gdata, SendArrRL, RecvArrRL, requests_RL, dir,
				size_send_RL, size_recv_RL, 1, q, SendRight, RecvLeft);
#else
#endif
	}


//	if(RecvLeft) {
//		if(dir == 0) {
//			int iter = 0;
//			for (int k = -rim; k <= gdata.mx[2]+rim; k++){
//				for (int j = -rim; j <= gdata.mx[1]+rim; j++){
//					for (int i = imaxRecv; i >= -rim; --i){
//						omb(i,j,k) = RecvArr_buff(iter);
//						iter++;
//					}
//				}
//			}
//		} else if (dir == 1) {
//			int iter = 0;
//			for (int j = imaxRecv; j >= -rim; --j){
//				for (int k = -rim; k <= gdata.mx[2]+rim; k++){
//					for (int i = -rim; i <= gdata.mx[0]+rim; i++){
//						omb(i,j,k) = RecvArr_buff(iter);
//						iter++;
//					}
//				}
//			}
//		} else if (dir == 2) {
//			int iter = 0;
//			for (int k = imaxRecv; k >= -rim; --k){
//				for (int j = -rim; j <= gdata.mx[1]+rim; j++){
//					for (int i = -rim; i <= gdata.mx[0]+rim; i++){
//						omb(i,j,k) = RecvArr_buff(iter);;
//						iter++;
//					}
//				}
//			}
//		}
//	}

}
#endif


void gridFunc::get_bcBuffRange(Data &gdata, int &range_min, int & range_max,
		int dir, int leftToRight, int q, int rim ) {
	if(leftToRight) {

		range_min = 1;
		range_max = rim;

		if(!gdata.get_EdgeGridding()) {
			range_min = 1;
			if(bc_Type[2*dir] == 1 || gdata.time < 0.01*gdata.dt) {
				range_min = 0;
			}
		}

	} else {

		range_min = -rim;
		range_max = -1;

		if((bc_Type[2*dir+1] > -1 || gdata.time < 0.01*gdata.dt) &&
				((dir==0 && q==q_Bx) || (dir==1 && q==q_By) || (dir==2 && q==q_Bz))) {
			range_max = -2;
		}

	}
}

#ifdef parallel
void gridFunc::bc_MPISend(Data &gdata, ProblemType &Problem,
		NumMatrix<REAL,3> &omb, int q, int rim) {



	bool SendLeft(false), RecvRight(false);

	// First part: from left to right (e.g., u(mx+1) = u(1) / u(0))

	// Size of send arrays
	int iminSend(1);
	int imaxSend(rim);
	int shift(0);

	// Size of recv arrays
	int iminRecv(1);
	int imaxRecv(rim);

	// For new gridding scheme size is always the same:
	if(gdata.get_EdgeGridding()) {
		// iminSend = 0;
		// imaxSend = rim-1;
		shift = 1;
	}





	// Begin with x-dimesion
	int size_send_xLR = 0;
	if(bc_Type[0] < 2) {
		SendLeft = true;

		get_bcBuffRange(gdata, iminSend, imaxSend, 0, 1, q, rim);
//		if(!gdata.get_EdgeGridding()) {
//			iminSend = 1;
//			if(bc_Type[0] == 1 || gdata.time < 0.01*gdata.dt) {
//				iminSend = 0;
//			}
//		}

		int iter=0;
		for (int k = -rim; k <= gdata.mx[2]+rim; k++){
			for (int j = -rim; j <= gdata.mx[1]+rim; j++){
				for (int i = iminSend; i <= imaxSend; i++){
					SendArr_xLR(iter) = omb(i-shift,j,k);    // x_min -> x_max
					iter++;
				}
			}
		}
		size_send_xLR = iter;

	} else {

		size_send_xLR = 0;

	}

	int size_recv_xLR = 0;
	if(bc_Type[1] < 2) {
		RecvRight = true;
		get_bcBuffRange(gdata, iminRecv, imaxRecv, 0, 1, q, rim);
//		if(!gdata.get_EdgeGridding()) {
//			if(bc_Type[1] == 1 || gdata.time < 0.01*gdata.dt) {
//				iminRecv = 0;
//			}
//		}
		size_recv_xLR = (imaxRecv-iminRecv+1)*(gdata.mx[1]+2*rim+1)*(gdata.mx[2]+2*rim+1);

	} else {
		size_recv_xLR = 0;
	}


	if(SendLeft || RecvRight) { // Periodic
#ifdef parallel
 		bc_ISendRecv(gdata, SendArr_xLR, RecvArr_xLR, requests_xLR,
				0, size_send_xLR, size_recv_xLR, 0, q, SendLeft, RecvRight);
#else
		bc_ISendRecv(gdata, SendArr_xLR, RecvArr_xLR,
				0, size_send_xLR, size_recv_xLR, 0, q, SendLeft, RecvRight);
#endif
	}


	// Next: y-dimension
	SendLeft = false;
	RecvRight = false;
	int size_send_yLR = 0;
	if(bc_Type[2] < 2) {
		SendLeft = true;
		get_bcBuffRange(gdata, iminSend, imaxSend, 1, 1, q, rim);
//		if(!gdata.get_EdgeGridding()) {
//			iminSend = 1;
//			if(bc_Type[2] == 1 || gdata.time < 0.01*gdata.dt) {
//				iminSend = 0;
//			}
//		}

		int iter=0;
		for (int j = iminSend; j <= imaxSend; j++){
			for (int k = -rim; k <= gdata.mx[2]+rim; k++){
				for (int i = -rim; i <= gdata.mx[0]+rim; i++){
					SendArr_yLR(iter) = omb(i,j-shift,k);    // y_min -> y_max
					iter++;
				}
			}
		}
		size_send_yLR = iter;

	} else {

		size_send_yLR = 0;

	}

	int size_recv_yLR = 0;
	if(bc_Type[3] < 2) {
		RecvRight = true;
		get_bcBuffRange(gdata, iminRecv, imaxRecv, 1, 1, q, rim);
//		if(!gdata.get_EdgeGridding()) {
//			if(bc_Type[3] == 1 || gdata.time < 0.01*gdata.dt) {
//				iminRecv = 0;
//			}
//		}
		size_recv_yLR = (gdata.mx[0]+2*rim+1)*(imaxRecv-iminRecv+1)*(gdata.mx[2]+2*rim+1);

	} else {
		size_recv_yLR = 0;
	}

	if(SendLeft || RecvRight) { // Periodic
#ifdef parallel
		bc_ISendRecv(gdata, SendArr_yLR, RecvArr_yLR, requests_yLR,
				1, size_send_yLR, size_recv_yLR, 0, q, SendLeft, RecvRight);
#else
		bc_ISendRecv(gdata, SendArr_buff, RecvArr_buff,
				1, size_send_yLR, size_recv_yLR, 0, q, SendLeft, RecvRight);
#endif
	}


	// Finally: z-dimension
	SendLeft = false;
	RecvRight = false;
	int size_send_zLR = 0;
	if(bc_Type[4] < 2) {
		SendLeft = true;
		get_bcBuffRange(gdata, iminSend, imaxSend, 2, 1, q, rim);
//		if(!gdata.get_EdgeGridding()) {
//			iminSend = 1;
//			if(bc_Type[2] == 1 || gdata.time < 0.01*gdata.dt) {
//				iminSend = 0;
//			}
//		}

		int iter=0;
		for (int k = iminSend; k <= imaxSend; k++){
			for (int j = -rim; j <= gdata.mx[1]+rim; j++){
				for (int i = -rim; i <= gdata.mx[0]+rim; i++){
					SendArr_zLR(iter) = omb(i,j,k-shift);    // z_max side
					iter++;
				}
			}
		}
		size_send_zLR = iter;

	} else {

		size_send_zLR = 0;
	}

	int size_recv_zLR = 0;
	if(bc_Type[5] < 2) {
		RecvRight = true;
		get_bcBuffRange(gdata, iminRecv, imaxRecv, 2, 1, q, rim);
//		if(!gdata.get_EdgeGridding()) {
//			if(bc_Type[5] == 1 || gdata.time < 0.01*gdata.dt) {
//				iminRecv = 0;
//			}
//		}
		size_recv_zLR = (gdata.mx[0]+2*rim+1)*(gdata.mx[1]+2*rim+1)*(imaxRecv-iminRecv+1);

	} else {
		size_recv_zLR = 0;
	}


	if(SendLeft || RecvRight) { // Periodic
#ifdef parallel
		bc_ISendRecv(gdata, SendArr_zLR, RecvArr_zLR, requests_zLR,
				2, size_send_zLR, size_recv_zLR, 0, q, SendLeft, RecvRight);
#else
		bc_ISendRecv(gdata, SendArr_zLR, RecvArr_zLR,
				2, size_send_zLR, size_recv_zLR, 0, q, SendLeft, RecvRight);
#endif
	}

//	if(bc_Type[0] > -1) { // User specific periodic bcs
			  // COMMENT: not implemented for 1D storage

//		cerr << " Not implemented, yet" << endl;
//		exit(3);
//		//			Problem.bc_periodic(gdata, *Send, omb, q, dir, false);
//	}
//	if(bc_Type[2] > -1) { // User specific periodic bcs
//		cerr << " Not implemented, yet" << endl;
//		exit(3);
//		//			Problem.bc_periodic(gdata, *Send, omb, q, dir, false);
//	}
//	if(bc_Type[4] > -1) { // User specific periodic bcs
//		cerr << " Not implemented, yet" << endl;
//		exit(3);
//		//			Problem.bc_periodic(gdata, *Send, omb, q, dir, false);
//	}




	// Second part: from right to left  (e.g., u(-1) = u(mx-1) / u(mx))



	bool SendRight(false), RecvLeft(false);

	iminSend = -rim; // Corresponding to mx[0]+i
	imaxSend = -1;
	shift = 0;
	imaxRecv = -1;
	if(gdata.get_EdgeGridding()) {
		shift = 1;
	}


	// Begin with x-dimension
	int size_send_xRL = 0;
	if(bc_Type[1] < 2) {
		SendRight = true;
		get_bcBuffRange(gdata, iminSend, imaxSend, 0, 0, q, rim);
//		if(gdata.get_EdgeGridding()) {
//			imaxSend =  -1;
//			if((bc_Type[1] > -1 || gdata.time < 0.01*gdata.dt) && q==q_Bx) {
//				imaxSend = -2;
//			}
//			shift=1;
//		} else {
//			imaxSend = -1;
//			if((bc_Type[1] > -1 || gdata.time < 0.01*gdata.dt) && q==q_Bx) {
//				imaxSend = -2;
//			}
//		}

		int iter=0;
		for (int k = -rim; k <= gdata.mx[2]+rim; k++){
			for (int j = -rim; j <= gdata.mx[1]+rim; j++){
				for (int i = imaxSend; i >= iminSend; --i){
					// x_max -> x_min
					SendArr_xRL(iter) = omb(gdata.mx[0]+i+shift,j,k);
					iter++;
				}
			}
		}
		size_send_xRL = iter;
	} else {
		SendRight = false;
		size_send_xRL = 0;
	}


	int size_recv_xRL = 0;
	if(bc_Type[0] < 2) {
		RecvLeft = true;
		get_bcBuffRange(gdata, iminRecv, imaxRecv, 0, 0, q, rim);
//		if(gdata.get_EdgeGridding()) {
//			if((bc_Type[0] > -1 || gdata.time < 0.01*gdata.dt) && q==q_Bx) {
//				imaxRecv = -2; // Corresponding to mx[0]+1
//			}
//		} else {
//			if((bc_Type[0] > -1 || gdata.time < 0.01*gdata.dt) && q==q_Bx) {
//				imaxRecv = -2;
//			}
//		}
		size_recv_xRL = (imaxRecv+rim+1)*(gdata.mx[1]+2*rim+1)*(gdata.mx[2]+2*rim+1);
	} else {
		RecvLeft = false;
		size_recv_xRL = 0;
	}

	if(SendRight || RecvLeft) { // Periodic
#ifdef parallel
		bc_ISendRecv(gdata, SendArr_xRL, RecvArr_xRL, requests_xRL,
				0, size_send_xRL, size_recv_xRL, 1, q, SendRight, RecvLeft);
#else
		bc_ISendRecv(gdata, SendArr_xRL, RecvArr_xRL,
				0, size_send_xRL, size_recv_xRL, 1, q, SendRight, RecvLeft);
#endif
	}



	// Proceed to y-dimension
	SendRight = false;
	RecvLeft = false;
	int size_send_yRL = 0;
	if(bc_Type[3] < 2) {
		SendRight = true;
		get_bcBuffRange(gdata, iminSend, imaxSend, 1, 0, q, rim);
//		if(gdata.get_EdgeGridding()) {
//			imaxSend =  -1;
//			if((bc_Type[3] > -1 || gdata.time < 0.01*gdata.dt) && q==q_By) {
//				imaxSend = -2;
//			}
//			shift=1;
//		} else {
//			imaxSend = -1;
//			if((bc_Type[3] > -1 || gdata.time < 0.01*gdata.dt) && q==q_By) {
//				imaxSend = -2;
//			}
//		}

		int iter=0;
		for (int j = imaxSend; j >= iminSend; --j){
			for (int k = -rim; k <= gdata.mx[2]+rim; k++){
				for (int i = -rim; i <= gdata.mx[0]+rim; i++){
					// y_max -> y_min
					SendArr_yRL(iter) = omb(i,gdata.mx[1]+j+shift,k);
					iter++;
				}
			}
		}
		size_send_yRL = iter;
	} else {
		SendRight = false;
		size_send_yRL = 0;
	}


	int size_recv_yRL = 0;
	if(bc_Type[2] < 2) {
		RecvLeft = true;
		get_bcBuffRange(gdata, iminRecv, imaxRecv, 1, 0, q, rim);
//		if(gdata.get_EdgeGridding()) {
//			if((bc_Type[2] > -1 || gdata.time < 0.01*gdata.dt) && q==q_By) {
//				imaxRecv = -2; // Corresponding to mx[0]+1
//			}
//		} else {
//			if((bc_Type[2] > -1 || gdata.time < 0.01*gdata.dt) && q==q_By) {
//				imaxRecv = -2;
//			}
//		}
		size_recv_yRL = (gdata.mx[0]+2*rim+1)*(imaxRecv+rim+1)*(gdata.mx[2]+2*rim+1);
	} else {
		RecvLeft = false;
		size_recv_yRL = 0;
	}

	if(SendRight || RecvLeft) { // Periodic
#ifdef parallel
		bc_ISendRecv(gdata, SendArr_yRL, RecvArr_yRL, requests_yRL,
				1, size_send_yRL, size_recv_yRL, 1, q, SendRight, RecvLeft);
#else
		bc_ISendRecv(gdata, SendArr_yRL, RecvArr_yRL,
				1, size_send_yRL, size_recv_yRL, 1, q, SendRight, RecvLeft);
#endif
	}



	// Finally: z-dimension
	SendRight = false;
	RecvLeft = false;
	int size_send_zRL = 0;
	if(bc_Type[5] < 2) {
		SendRight = true;
		get_bcBuffRange(gdata, iminSend, imaxSend, 2, 0, q, rim);
//		if(gdata.get_EdgeGridding()) {
//			imaxSend =  -1;
//			if((bc_Type[5] > -1 || gdata.time < 0.01*gdata.dt) && q==q_Bz) {
//				imaxSend = -2;
//			}
//			shift=1;
//		} else {
//			imaxSend = -1;
//			if((bc_Type[5] > -1 || gdata.time < 0.01*gdata.dt) && q==q_Bz) {
//				imaxSend = -2;
//			}
//		}

		int iter=0;
		for (int k = imaxSend; k >= iminSend; --k){
			for (int j = -rim; j <= gdata.mx[1]+rim; j++){
				for (int i = -rim; i <= gdata.mx[0]+rim; i++){
					// z_max -> z_min
					SendArr_zRL(iter) = omb(i,j,gdata.mx[2]+k+shift);
					iter++;
				}
			}
		}
		size_send_zRL = iter;
	} else {
		SendRight = false;
		size_send_zRL = 0;
	}


	int size_recv_zRL = 0;
	if(bc_Type[4] < 2) {
		RecvLeft = true;
		get_bcBuffRange(gdata, iminRecv, imaxRecv, 2, 0, q, rim);
//		if(gdata.get_EdgeGridding()) {
//			if((bc_Type[4] > -1 || gdata.time < 0.01*gdata.dt) && q==q_Bz) {
//				imaxRecv = -2; // Corresponding to mx[0]+1
//			}
//		} else {
//			if((bc_Type[4] > -1 || gdata.time < 0.01*gdata.dt) && q==q_Bz) {
//				imaxRecv = -2;
//			}
//		}
		size_recv_zRL = (gdata.mx[0]+2*rim+1)*(gdata.mx[1]+2*rim+1)*(imaxRecv+rim+1);
	} else {
		RecvLeft = false;
		size_recv_zRL = 0;
	}

	if(SendRight || RecvLeft) { // Periodic
#ifdef parallel
		bc_ISendRecv(gdata, SendArr_zRL, RecvArr_zRL, requests_zRL,
				2, size_send_zRL, size_recv_zRL, 1, q, SendRight, RecvLeft);
#else
		bc_ISendRecv(gdata, SendArr_zRL, RecvArr_zRL,
				2, size_send_zRL, size_recv_zRL, 1, q, SendRight, RecvLeft);
#endif
	}

}
#endif



#ifdef parallel
void gridFunc::bc_MPIRecv(Data &gdata, ProblemType &Problem,
		NumMatrix<REAL,3> &omb, int dir, int q, int rim) {

	// First part: from left to right
	bool SendLeft(false), RecvRight(false);

	if(bc_Type[2*dir] < 2) SendLeft = true;
	if(bc_Type[2*dir+1] < 2) RecvRight = true;

	int iminRecv(1);
	int imaxRecv(rim);
	int size_recv(0);
	if(RecvRight) {  // Receiver is Periodic
		if(gdata.get_EdgeGridding()) {
			// iminRecv = 1; // Corresponding to mx[0]+1
			// iminRecv = 0; // Corresponding to mx[0]+1
			// imaxRecv = rim-1;
			iminRecv = 1; // Corresponding to mx[0]+1
			imaxRecv = rim;
		} else {
			if(bc_Type[2*dir+1] == 1 || gdata.time < 0.01*gdata.dt) {
				iminRecv = 0;
			}
		}
	}

	// Ensure that communication from left to right has completely finished:
	if(SendLeft || RecvRight) { // Periodic
		int numRequests=0;
		if(SendLeft) numRequests++;
		if(RecvRight) numRequests++;
		//MPI_Status statusrl[numRequests];
		std::vector<MPI_Status> statusrl(numRequests);
		MPI_Waitall(numRequests, requests_LR, statusrl.data());
	}


	if(RecvRight) {
//		int numRequests=1;
//
//		// Ensure that communication has finished
//		if(SendLeft) numRequests++;
//		MPI_Status statusrl[numRequests];
//		MPI_Waitall(numRequests, requests_LR, statusrl);


		if(dir == 0) {
			int iter = 0;
			for (int k = -rim; k <= gdata.mx[2]+rim; k++){
				for (int j = -rim; j <= gdata.mx[1]+rim; j++){
					for (int i = iminRecv; i <= imaxRecv; i++){
						omb(gdata.mx[0]+i,j,k) = RecvArrLR(iter);
						iter++;
					}
				}
			}
		} else if (dir == 1) {
			int iter = 0;
			for (int j = iminRecv; j <= imaxRecv; j++){
				for (int k = -rim; k <= gdata.mx[2]+rim; k++){
					for (int i = -rim; i <= gdata.mx[0]+rim; i++){
						// y_min -> y_max
						omb(i,gdata.mx[1]+j,k) = RecvArrLR(iter);
						iter++;
					}
				}
			}
		} else if (dir == 2) {
			int iter = 0;
			for (int k = iminRecv; k <= imaxRecv; k++){
				for (int j = -rim; j <= gdata.mx[1]+rim; j++){
					for (int i = -rim; i <= gdata.mx[0]+rim; i++){
						// z_min -> z_max
						omb(i,j,gdata.mx[2]+k) = RecvArrLR(iter);
						iter++;
					}
				}
			}
		}
	}

	// Second part
	bool SendRight(false), RecvLeft(false);

	if(bc_Type[2*dir+1] < 2) SendRight = true;
	if(bc_Type[2*dir]   < 2) RecvLeft  = true;

	imaxRecv = -1;
	if(RecvLeft) {
		if(gdata.get_EdgeGridding()) {
			//if((bc_Type[2*dir] > -2 || gdata.time < 0.01*gdata.dt) &&
			if((bc_Type[2*dir] > -1 || gdata.time < 0.01*gdata.dt) &&
					((dir==0 && q==q_Bx) || (dir==1 && q==q_By) ||
							(dir==2 && q==q_Bz))) {
				imaxRecv = -2; // Corresponding to mx[0]+1
			}
		} else {
			if((bc_Type[2*dir] > -1 || gdata.time < 0.01*gdata.dt) &&
					((dir==0 && q==q_Bx) || (dir==1 && q==q_By) ||
							(dir==2 && q==q_Bz))) {
				imaxRecv = -2;
			}
		}
	}


	// Ensure that communication from left to right has completely finished:
	if(SendRight || RecvLeft) { // Periodic
		int numRequests=0;
		if(SendRight) numRequests++;
		if(RecvLeft) numRequests++;
		//MPI_Status statusrl[numRequests];
		std::vector<MPI_Status> statusrl(numRequests);
		MPI_Waitall(numRequests, requests_RL, statusrl.data());
	}


	if(RecvLeft) {
//		int numRequests=1;
//
//		// Ensure that communication has finished
//		if(SendRight) numRequests++;
//		MPI_Status statusrl[numRequests];
//		MPI_Waitall(numRequests, requests_RL, statusrl);

		if(dir == 0) {
			int iter = 0;
			for (int k = -rim; k <= gdata.mx[2]+rim; k++){
				for (int j = -rim; j <= gdata.mx[1]+rim; j++){
					for (int i = imaxRecv; i >= -rim; --i){
						omb(i,j,k) = RecvArrRL(iter);
						iter++;
					}
				}
			}
		} else if (dir == 1) {
			int iter = 0;
			for (int j = imaxRecv; j >= -rim; --j){
				for (int k = -rim; k <= gdata.mx[2]+rim; k++){
					for (int i = -rim; i <= gdata.mx[0]+rim; i++){
						omb(i,j,k) = RecvArrRL(iter);
						iter++;
					}
				}
			}
		} else if (dir == 2) {
			int iter = 0;
			for (int k = imaxRecv; k >= -rim; --k){
				for (int j = -rim; j <= gdata.mx[1]+rim; j++){
					for (int i = -rim; i <= gdata.mx[0]+rim; i++){
						omb(i,j,k) = RecvArrRL(iter);
						iter++;
					}
				}
			}
		}
	}

}
#endif


#ifdef parallel
void gridFunc::bc_MPIRecv(Data &gdata, ProblemType &Problem,
		NumMatrix<REAL,3> &omb, int q, int rim) {

	// First write data transferred from left to right


	// x-direction
	int iminRecv, imaxRecv;
	if(bc_Type[1] < 2) {
		int numRequests=1;

		// Wait for communication to finish:
		if(bc_Type[0] < 2) numRequests++;
		//MPI_Status statusrl[numRequests];
		std::vector<MPI_Status> statusrl(numRequests);
		MPI_Waitall(numRequests, requests_xLR, statusrl.data());

		get_bcBuffRange(gdata, iminRecv, imaxRecv, 0, 1, q, rim);
		int iter = 0;
		for (int k = -rim; k <= gdata.mx[2]+rim; k++){
			for (int j = -rim; j <= gdata.mx[1]+rim; j++){
				for (int i = iminRecv; i <= imaxRecv; i++){
					omb(gdata.mx[0]+i,j,k) = RecvArr_xLR(iter);
					iter++;
				}
			}
		}
	}

	// y-direction
	if(bc_Type[3] < 2) {
		int numRequests=1;

		// Wait for communication to finish:
		if(bc_Type[2] < 2) numRequests++;
		//MPI_Status statusrl[numRequests];
		std::vector<MPI_Status> statusrl(numRequests);
		MPI_Waitall(numRequests, requests_yLR, statusrl.data());

		get_bcBuffRange(gdata, iminRecv, imaxRecv, 1, 1, q, rim);
		int iter = 0;
		for (int j = iminRecv; j <= imaxRecv; j++){
			for (int k = -rim; k <= gdata.mx[2]+rim; k++){
				for (int i = -rim; i <= gdata.mx[0]+rim; i++){
					// y_min -> y_max
					omb(i,gdata.mx[1]+j,k) = RecvArr_yLR(iter);
					iter++;
				}
			}
		}
	}


	// z-direction
	if(bc_Type[5] < 2) {
		int numRequests=1;

		// Wait for communication to finish:
		if(bc_Type[4] < 2) numRequests++;
		//MPI_Status statusrl[numRequests];
		std::vector<MPI_Status> statusrl(numRequests);
		MPI_Waitall(numRequests, requests_zLR, statusrl.data());

		get_bcBuffRange(gdata, iminRecv, imaxRecv, 2, 1, q, rim);
		int iter = 0;
		for (int k = iminRecv; k <= imaxRecv; k++){
			for (int j = -rim; j <= gdata.mx[1]+rim; j++){
				for (int i = -rim; i <= gdata.mx[0]+rim; i++){
					// z_min -> z_max
					omb(i,j,gdata.mx[2]+k) = RecvArr_zLR(iter);
					iter++;
				}
			}
		}
	}



	// Now write data transferred from right to left

	// x-direction
	if(bc_Type[0] < 2) {
		int numRequests=1;

		// Wait for communication to finish:
		if(bc_Type[1] < 2) numRequests++;
		//MPI_Status statusrl[numRequests];
		std::vector<MPI_Status> statusrl(numRequests);
		MPI_Waitall(numRequests, requests_xRL, statusrl.data());

		get_bcBuffRange(gdata, iminRecv, imaxRecv, 0, 0, q, rim);
		int iter = 0;
		for (int k = -rim; k <= gdata.mx[2]+rim; k++){
			for (int j = -rim; j <= gdata.mx[1]+rim; j++){
				for (int i = imaxRecv; i >= -rim; --i){
					omb(i,j,k) = RecvArr_xRL(iter);
					iter++;
				}
			}
		}
	}


	// y-direction
	if(bc_Type[2] < 2) {
		int numRequests=1;

		// Wait for communication to finish:
		if(bc_Type[3] < 2) numRequests++;
		//MPI_Status statusrl[numRequests];
		std::vector<MPI_Status> statusrl(numRequests);
		MPI_Waitall(numRequests, requests_yRL, statusrl.data());

		get_bcBuffRange(gdata, iminRecv, imaxRecv, 1, 0, q, rim);
		int iter = 0;
		for (int j = imaxRecv; j >= -rim; --j){
			for (int k = -rim; k <= gdata.mx[2]+rim; k++){
				for (int i = -rim; i <= gdata.mx[0]+rim; i++){
					omb(i,j,k) = RecvArr_yRL(iter);
					iter++;
				}
			}
		}
	}


	// z-direction
	if(bc_Type[4] < 2) {
		int numRequests=1;

		// Wait for communication to finish:
		if(bc_Type[5] < 2) numRequests++;
		//MPI_Status statusrl[numRequests];
		std::vector<MPI_Status> statusrl(numRequests);
		MPI_Waitall(numRequests, requests_zRL, statusrl.data());


		get_bcBuffRange(gdata, iminRecv, imaxRecv, 2, 0, q, rim);
		int iter = 0;
		for (int k = imaxRecv; k >= -rim; --k){
			for (int j = -rim; j <= gdata.mx[1]+rim; j++){
				for (int i = -rim; i <= gdata.mx[0]+rim; i++){
					omb(i,j,k) = RecvArr_zRL(iter);;
					iter++;
				}
			}
		}
	}

}
#endif


//void gridFunc::bc_Periodic_new(Data &gdata, ProblemType &Problem,
//                           NumMatrix<REAL,3> &omb,
//                           int dir, int q, int rim, bool internal_only) {
//
//
//	NumArray<double> *Recv;
//	NumArray<double> *Send;
//
//	bool SendLeft(false), RecvRight(false);
//
//	if(bc_Type[2*dir] < 2) SendLeft = true;
//	if(bc_Type[2*dir+1] < 2) RecvRight = true;
//
//	if(internal_only) {
//#ifdef parallel
//		if(gdata.check_MpiLowest(dir)) {
//			SendLeft = false;
//		}
//		if(gdata.check_MpiHighest(dir)) {
//			RecvRight = false;
//		}
//#else
//		SendLeft = false;
//		RecvRight = false;
//#endif
//	}
//
//	// First part: from left to right (e.g., u(mx+1) = u(1) / u(0))
//	if (SendLeft) { // Sender is Periodic
//
//		int iminSend(1);
//		int imaxSend(rim);
//		int shift(0);
//
//		// For new gridding scheme size is always the same:
//		if(gdata.get_EdgeGridding()) {
//			// iminSend = 0;
//			// imaxSend = rim-1;
//			iminSend = 1;
//			imaxSend = rim;
//			shift = 1;
//		} else {
//			iminSend = 1;
//			imaxSend = rim;
//			if(bc_Type[2*dir] == 1 || gdata.time < 0.01*gdata.dt) {
//				iminSend = 0;
//			}
//		}
//		if(dir == 0) {
//			Send = &SendArr_x;
//			int iter=0;
//			for (int k = -rim; k <= gdata.mx[2]+rim; k++){
//				for (int j = -rim; j <= gdata.mx[1]+rim; j++){
//					for (int i = iminSend; i <= imaxSend; i++){
//						Send(iter) = omb(i-shift,j,k);    // x_min -> x_max
//						iter++;
//					}
//				}
//			}
//
//
//		} else if (dir == 1) {
//			Send = &SendArr_y;
//			int iter=0;
//			for (int j = iminSend; j <= imaxSend; j++){
//				for (int k = -rim; k <= gdata.mx[2]+rim; k++){
//					for (int i = -rim; i <= gdata.mx[0]+rim; i++){
//						Send(iter) = omb(i,j-shift,k);    // y_min -> y_max
//						iter++;
//					}
//				}
//			}
//
//		} else if (dir == 2) {
//			Send = &SendArr_z;
//			int iter=0;
//			for (int k = iminSend; k <= imaxSend; k++){
//				for (int j = -rim; j <= gdata.mx[1]+rim; j++){
//					for (int i = -rim; i <= gdata.mx[0]+rim; i++){
//						Send(iter) = omb(i,j,k-shift);    // z_max side
//						iter++;
//					}
//				}
//			}
//
//		}
//
//		if(bc_Type[2*dir] > -1) { // User specific periodic bcs
//			cerr << " Not implemented, yet" << endl;
//			exit(3);
////			Problem.bc_periodic(gdata, *Send, omb, q, dir, false);
//		}
//
//	} else {
//
//		Send = &SendArr_empty;
//
////		Send.resize(Index::set(0,0,0), Index::set(0,0,0));
//	}
//
//
//	int iminRecv(1);
//	int imaxRecv(rim);
//	if(RecvRight) {  // Receiver is Periodic
//		if(gdata.get_EdgeGridding()) {
//			// iminRecv = 1; // Corresponding to mx[0]+1
//			// iminRecv = 0; // Corresponding to mx[0]+1
//			// imaxRecv = rim-1;
//			iminRecv = 1; // Corresponding to mx[0]+1
//			imaxRecv = rim;
//		} else {
//			if(bc_Type[2*dir+1] == 1 || gdata.time < 0.01*gdata.dt) {
//				iminRecv = 0;
//			}
//		}
//
//		if(dir == 0) {
//			Recv = &RecvArr_x;
//		} else if (dir == 1) {
//			Recv = &RecvArr_y;
//		} else if (dir == 2) {
//			Recv = &RecvArr_z;
//		}
//
//	} else {
//		Recv = &Recv_empty;
////		Recv.resize(Index::set(0,0,0), Index::set(0,0,0));
//	}
//
//
//	if(SendLeft || RecvRight) { // Periodic
//		bc_SendRecv(gdata, *Send, *Recv, dir, 0, q, SendLeft, RecvRight);
//	}
//
//
//
//
//	if(RecvRight) {
//		if(dir == 0) {
//			for (int k = -rim; k <= gdata.mx[2]+rim; k++){
//				for (int j = -rim; j <= gdata.mx[1]+rim; j++){
//					for (int i = iminRecv; i <= imaxRecv; i++){
//						omb(gdata.mx[0]+i,j,k) = Recv_x(i,j,k);
//					}
//				}
//			}
//		} else if (dir == 1) {
//			for (int j = iminRecv; j <= imaxRecv; j++){
//				for (int k = -rim; k <= gdata.mx[2]+rim; k++){
//					for (int i = -rim; i <= gdata.mx[0]+rim; i++){
//						// y_min -> y_max
//						omb(i,gdata.mx[1]+j,k) = Recv_y(i,j,k);
//					}
//				}
//			}
//		} else if (dir == 2) {
//			for (int k = iminRecv; k <= imaxRecv; k++){
//				for (int j = -rim; j <= gdata.mx[1]+rim; j++){
//					for (int i = -rim; i <= gdata.mx[0]+rim; i++){
//						// z_min -> z_max
//						omb(i,j,gdata.mx[2]+k) = Recv_z(i,j,k);
//					}
//				}
//			}
//		}
//	}
//
//
//	bool SendRight(false), RecvLeft(false);
//
//	if(bc_Type[2*dir+1] < 2) SendRight = true;
//	if(bc_Type[2*dir]   < 2) RecvLeft  = true;
//
//	if(internal_only) {
//#ifdef parallel
//		if(gdata.check_MpiLowest(dir)) {
//			RecvLeft = false;
//		}
//		if(gdata.check_MpiHighest(dir)) {
//			SendRight = false;
//		}
//#else
//		SendLeft = false;
//		RecvRight = false;
//#endif
//	}
//
//	// Second part: from right to left (e.g., u(-1) = u(mx-1) / u(mx))
//	if(SendRight) { // Periodic
//
//		int iminSend(-rim); // Corresponding to mx[0]+i
//		int imaxSend(-1);
//		int shift(0);
//
//		// For new gridding scheme size is always the same:
//		if(gdata.get_EdgeGridding()) {
//			// iminSend = -rim+1;
//			// imaxSend =  0;
//			iminSend = -rim;
//			imaxSend =  -1;
//			if((bc_Type[2*dir+1] > -1 || gdata.time < 0.01*gdata.dt) &&
//			   //if((bc_Type[2*dir+1] > -2 || gdata.time < 0.01*gdata.dt) &&
//			   ((dir==0 && q==q_Bx) || (dir==1 && q==q_By) ||
//			    (dir==2 && q==q_Bz))) {
//				// imaxSend = -1;
//				imaxSend = -2;
//			}
//			// shift=-1;
//			shift=1;
//		} else {
//			iminSend = -rim;
//			imaxSend = -1;
//			if((bc_Type[2*dir+1] > -1 || gdata.time < 0.01*gdata.dt) &&
//			   ((dir==0 && q==q_Bx) || (dir==1 && q==q_By) ||
//			    (dir==2 && q==q_Bz))) {
//				imaxSend = -2;
//			}
//		}
//		if(dir == 0) {
//			Send.resize(Index::set(iminSend,-rim,-rim),
//			            Index::set(imaxSend,gdata.mx[1]+rim,
//			                       gdata.mx[2]+rim));
//
//			for (int k = -rim; k <= gdata.mx[2]+rim; k++){
//				for (int j = -rim; j <= gdata.mx[1]+rim; j++){
//					for (int i = imaxSend; i >= iminSend; --i){
//						// x_max -> x_min
//						Send(i,j,k) = omb(gdata.mx[0]+i+shift,j,k);
//					}
//				}
//			}
//
//		} else if (dir == 1) {
//			Send.resize(Index::set(-rim,iminSend,-rim),
//			            Index::set(gdata.mx[0]+rim,imaxSend,
//			                       gdata.mx[2]+rim));
//
//			for (int j = imaxSend; j >= iminSend; --j){
//				for (int k = -rim; k <= gdata.mx[2]+rim; k++){
//					for (int i = -rim; i <= gdata.mx[0]+rim; i++){
//						// y_max -> y_min
//						Send(i,j,k) = omb(i,gdata.mx[1]+j+shift,k);
//					}
//				}
//			}
//
//		} else if (dir == 2) {
//			Send.resize(Index::set(-rim,-rim, iminSend),
//			            Index::set(gdata.mx[0]+rim,gdata.mx[1]+rim,
//			                       imaxSend));
//
//			for (int k = imaxSend; k >= iminSend; --k){
//				for (int j = -rim; j <= gdata.mx[1]+rim; j++){
//					for (int i = -rim; i <= gdata.mx[0]+rim; i++){
//						// z_max -> z_min
//						Send(i,j,k) = omb(i,j,gdata.mx[2]+k+shift);
//					}
//				}
//			}
//
//		}
//
//		if(bc_Type[2*dir+1] > -1) {
//			Problem.bc_periodic(gdata, Send, omb, q, dir, true);
//		}
//	} else {
//		Send.resize(Index::set(0,0,0), Index::set(0,0,0));
//	}
//
//	imaxRecv = -1;
//	if(RecvLeft) {
//		if(gdata.get_EdgeGridding()) {
//			//if((bc_Type[2*dir] > -2 || gdata.time < 0.01*gdata.dt) &&
//			if((bc_Type[2*dir] > -1 || gdata.time < 0.01*gdata.dt) &&
//			   ((dir==0 && q==q_Bx) || (dir==1 && q==q_By) ||
//			    (dir==2 && q==q_Bz))) {
//				imaxRecv = -2; // Corresponding to mx[0]+1
//			}
//		} else {
//			if((bc_Type[2*dir] > -1 || gdata.time < 0.01*gdata.dt) &&
//			   ((dir==0 && q==q_Bx) || (dir==1 && q==q_By) ||
//			    (dir==2 && q==q_Bz))) {
//				imaxRecv = -2;
//			}
//		}
//		if(dir == 0) {
//			Recv.resize(Index::set(-rim,-rim,-rim),
//			            Index::set(imaxRecv,gdata.mx[1]+rim,
//			                       gdata.mx[2]+rim));
//		} else if (dir == 1) {
//			Recv.resize(Index::set(-rim,-rim,-rim),
//			            Index::set(gdata.mx[0]+rim,imaxRecv,
//			                       gdata.mx[2]+rim));
//		} else if (dir == 2) {
//			Recv.resize(Index::set(-rim,-rim,-rim),
//			            Index::set(gdata.mx[0]+rim,gdata.mx[1]+rim,
//			                       imaxRecv));
//		}
//	} else {
//		Recv.resize(Index::set(0,0,0), Index::set(0,0,0));
//	}
//
//	if(SendRight || RecvLeft) {
//		bc_SendRecv(gdata, Send, Recv, dir, 1, q, SendRight, RecvLeft);
//	}
//
//
//	if(RecvLeft) {
//		if(dir == 0) {
//			for (int k = -rim; k <= gdata.mx[2]+rim; k++){
//				for (int j = -rim; j <= gdata.mx[1]+rim; j++){
//					for (int i = imaxRecv; i >= -rim; --i){
//						omb(i,j,k) = Recv(i,j,k);
//					}
//				}
//			}
//		} else if (dir == 1) {
//			for (int j = imaxRecv; j >= -rim; --j){
//				for (int k = -rim; k <= gdata.mx[2]+rim; k++){
//					for (int i = -rim; i <= gdata.mx[0]+rim; i++){
//						omb(i,j,k) = Recv(i,j,k);
//					}
//				}
//			}
//		} else if (dir == 2) {
//			for (int k = imaxRecv; k >= -rim; --k){
//				for (int j = -rim; j <= gdata.mx[1]+rim; j++){
//					for (int i = -rim; i <= gdata.mx[0]+rim; i++){
//						omb(i,j,k) = Recv(i,j,k);
//					}
//				}
//			}
//		}
//	}
//
//}
//void gridFunc::bc_Periodic_new(Data &gdata, ProblemType &Problem,
//                           NumMatrix<REAL,3> &omb,
//                           int dir, int q, int rim, bool internal_only) {
//
//
//	NumMatrix<double,3> *Recv;
//	NumMatrix<double,3> *Send;
//
//	bool SendLeft(false), RecvRight(false);
//
//	if(bc_Type[2*dir] < 2) SendLeft = true;
//	if(bc_Type[2*dir+1] < 2) RecvRight = true;
//
//	if(internal_only) {
//#ifdef parallel
//		if(gdata.check_MpiLowest(dir)) {
//			SendLeft = false;
//		}
//		if(gdata.check_MpiHighest(dir)) {
//			RecvRight = false;
//		}
//#else
//		SendLeft = false;
//		RecvRight = false;
//#endif
//	}
//
//	// First part: from left to right (e.g., u(mx+1) = u(1) / u(0))
//	if (SendLeft) { // Sender is Periodic
//
//		int iminSend(1);
//		int imaxSend(rim);
//		int shift(0);
//
//		// For new gridding scheme size is always the same:
//		if(gdata.get_EdgeGridding()) {
//			// iminSend = 0;
//			// imaxSend = rim-1;
//			iminSend = 1;
//			imaxSend = rim;
//			shift = 1;
//		} else {
//			iminSend = 1;
//			imaxSend = rim;
//			if(bc_Type[2*dir] == 1 || gdata.time < 0.01*gdata.dt) {
//				iminSend = 0;
//			}
//		}
//		if(dir == 0) {
//			for (int k = -rim; k <= gdata.mx[2]+rim; k++){
//				for (int j = -rim; j <= gdata.mx[1]+rim; j++){
//					for (int i = iminSend; i <= imaxSend; i++){
//						Send_x(i,j,k) = omb(i-shift,j,k);    // x_min -> x_max
//					}
//				}
//			}
//			Send = &Send_x;
//
//
//		} else if (dir == 1) {
//			for (int j = iminSend; j <= imaxSend; j++){
//				for (int k = -rim; k <= gdata.mx[2]+rim; k++){
//					for (int i = -rim; i <= gdata.mx[0]+rim; i++){
//						Send_y(i,j,k) = omb(i,j-shift,k);    // y_min -> y_max
//					}
//				}
//			}
//			Send = &Send_y;
//
//		} else if (dir == 2) {
//			for (int k = iminSend; k <= imaxSend; k++){
//				for (int j = -rim; j <= gdata.mx[1]+rim; j++){
//					for (int i = -rim; i <= gdata.mx[0]+rim; i++){
//						Send_z(i,j,k) = omb(i,j,k-shift);    // z_max side
//					}
//				}
//			}
//			Send = &Send_z;
//
//		}
//
//		if(bc_Type[2*dir] > -1) { // User specific periodic bcs
//			Problem.bc_periodic(gdata, *Send, omb, q, dir, false);
//		}
//
//	} else {
//
//		Send = &Send_empty;
//
////		Send.resize(Index::set(0,0,0), Index::set(0,0,0));
//	}
//
//
//	int iminRecv(1);
//	int imaxRecv(rim);
//	if(RecvRight) {  // Receiver is Periodic
//		if(gdata.get_EdgeGridding()) {
//			// iminRecv = 1; // Corresponding to mx[0]+1
//			// iminRecv = 0; // Corresponding to mx[0]+1
//			// imaxRecv = rim-1;
//			iminRecv = 1; // Corresponding to mx[0]+1
//			imaxRecv = rim;
//		} else {
//			if(bc_Type[2*dir+1] == 1 || gdata.time < 0.01*gdata.dt) {
//				iminRecv = 0;
//			}
//		}
//
//		if(dir == 0) {
//			Recv = &Recv_x;
//		} else if (dir == 1) {
//			Recv = &Recv_y;
//		} else if (dir == 2) {
//			Recv = &Recv_z;
//		}
//
//	} else {
//		Recv = &Recv_empty;
////		Recv.resize(Index::set(0,0,0), Index::set(0,0,0));
//	}
//
//
//	if(SendLeft || RecvRight) { // Periodic
//		bc_SendRecv(gdata, *Send, *Recv, dir, 0, q, SendLeft, RecvRight);
//	}
//
//
//
//
//	if(RecvRight) {
//		if(dir == 0) {
//			for (int k = -rim; k <= gdata.mx[2]+rim; k++){
//				for (int j = -rim; j <= gdata.mx[1]+rim; j++){
//					for (int i = iminRecv; i <= imaxRecv; i++){
//						omb(gdata.mx[0]+i,j,k) = Recv_x(i,j,k);
//					}
//				}
//			}
//		} else if (dir == 1) {
//			for (int j = iminRecv; j <= imaxRecv; j++){
//				for (int k = -rim; k <= gdata.mx[2]+rim; k++){
//					for (int i = -rim; i <= gdata.mx[0]+rim; i++){
//						// y_min -> y_max
//						omb(i,gdata.mx[1]+j,k) = Recv_y(i,j,k);
//					}
//				}
//			}
//		} else if (dir == 2) {
//			for (int k = iminRecv; k <= imaxRecv; k++){
//				for (int j = -rim; j <= gdata.mx[1]+rim; j++){
//					for (int i = -rim; i <= gdata.mx[0]+rim; i++){
//						// z_min -> z_max
//						omb(i,j,gdata.mx[2]+k) = Recv_z(i,j,k);
//					}
//				}
//			}
//		}
//	}
//
//
//	bool SendRight(false), RecvLeft(false);
//
//	if(bc_Type[2*dir+1] < 2) SendRight = true;
//	if(bc_Type[2*dir]   < 2) RecvLeft  = true;
//
//	if(internal_only) {
//#ifdef parallel
//		if(gdata.check_MpiLowest(dir)) {
//			RecvLeft = false;
//		}
//		if(gdata.check_MpiHighest(dir)) {
//			SendRight = false;
//		}
//#else
//		SendLeft = false;
//		RecvRight = false;
//#endif
//	}
//
//	// Second part: from right to left (e.g., u(-1) = u(mx-1) / u(mx))
//	if(SendRight) { // Periodic
//
//		int iminSend(-rim); // Corresponding to mx[0]+i
//		int imaxSend(-1);
//		int shift(0);
//
//		// For new gridding scheme size is always the same:
//		if(gdata.get_EdgeGridding()) {
//			// iminSend = -rim+1;
//			// imaxSend =  0;
//			iminSend = -rim;
//			imaxSend =  -1;
//			if((bc_Type[2*dir+1] > -1 || gdata.time < 0.01*gdata.dt) &&
//			   //if((bc_Type[2*dir+1] > -2 || gdata.time < 0.01*gdata.dt) &&
//			   ((dir==0 && q==q_Bx) || (dir==1 && q==q_By) ||
//			    (dir==2 && q==q_Bz))) {
//				// imaxSend = -1;
//				imaxSend = -2;
//			}
//			// shift=-1;
//			shift=1;
//		} else {
//			iminSend = -rim;
//			imaxSend = -1;
//			if((bc_Type[2*dir+1] > -1 || gdata.time < 0.01*gdata.dt) &&
//			   ((dir==0 && q==q_Bx) || (dir==1 && q==q_By) ||
//			    (dir==2 && q==q_Bz))) {
//				imaxSend = -2;
//			}
//		}
//		if(dir == 0) {
//			Send.resize(Index::set(iminSend,-rim,-rim),
//			            Index::set(imaxSend,gdata.mx[1]+rim,
//			                       gdata.mx[2]+rim));
//
//			for (int k = -rim; k <= gdata.mx[2]+rim; k++){
//				for (int j = -rim; j <= gdata.mx[1]+rim; j++){
//					for (int i = imaxSend; i >= iminSend; --i){
//						// x_max -> x_min
//						Send(i,j,k) = omb(gdata.mx[0]+i+shift,j,k);
//					}
//				}
//			}
//
//		} else if (dir == 1) {
//			Send.resize(Index::set(-rim,iminSend,-rim),
//			            Index::set(gdata.mx[0]+rim,imaxSend,
//			                       gdata.mx[2]+rim));
//
//			for (int j = imaxSend; j >= iminSend; --j){
//				for (int k = -rim; k <= gdata.mx[2]+rim; k++){
//					for (int i = -rim; i <= gdata.mx[0]+rim; i++){
//						// y_max -> y_min
//						Send(i,j,k) = omb(i,gdata.mx[1]+j+shift,k);
//					}
//				}
//			}
//
//		} else if (dir == 2) {
//			Send.resize(Index::set(-rim,-rim, iminSend),
//			            Index::set(gdata.mx[0]+rim,gdata.mx[1]+rim,
//			                       imaxSend));
//
//			for (int k = imaxSend; k >= iminSend; --k){
//				for (int j = -rim; j <= gdata.mx[1]+rim; j++){
//					for (int i = -rim; i <= gdata.mx[0]+rim; i++){
//						// z_max -> z_min
//						Send(i,j,k) = omb(i,j,gdata.mx[2]+k+shift);
//					}
//				}
//			}
//
//		}
//
//		if(bc_Type[2*dir+1] > -1) {
//			Problem.bc_periodic(gdata, Send, omb, q, dir, true);
//		}
//	} else {
//		Send.resize(Index::set(0,0,0), Index::set(0,0,0));
//	}
//
//	imaxRecv = -1;
//	if(RecvLeft) {
//		if(gdata.get_EdgeGridding()) {
//			//if((bc_Type[2*dir] > -2 || gdata.time < 0.01*gdata.dt) &&
//			if((bc_Type[2*dir] > -1 || gdata.time < 0.01*gdata.dt) &&
//			   ((dir==0 && q==q_Bx) || (dir==1 && q==q_By) ||
//			    (dir==2 && q==q_Bz))) {
//				imaxRecv = -2; // Corresponding to mx[0]+1
//			}
//		} else {
//			if((bc_Type[2*dir] > -1 || gdata.time < 0.01*gdata.dt) &&
//			   ((dir==0 && q==q_Bx) || (dir==1 && q==q_By) ||
//			    (dir==2 && q==q_Bz))) {
//				imaxRecv = -2;
//			}
//		}
//		if(dir == 0) {
//			Recv.resize(Index::set(-rim,-rim,-rim),
//			            Index::set(imaxRecv,gdata.mx[1]+rim,
//			                       gdata.mx[2]+rim));
//		} else if (dir == 1) {
//			Recv.resize(Index::set(-rim,-rim,-rim),
//			            Index::set(gdata.mx[0]+rim,imaxRecv,
//			                       gdata.mx[2]+rim));
//		} else if (dir == 2) {
//			Recv.resize(Index::set(-rim,-rim,-rim),
//			            Index::set(gdata.mx[0]+rim,gdata.mx[1]+rim,
//			                       imaxRecv));
//		}
//	} else {
//		Recv.resize(Index::set(0,0,0), Index::set(0,0,0));
//	}
//
//	if(SendRight || RecvLeft) {
//		bc_SendRecv(gdata, Send, Recv, dir, 1, q, SendRight, RecvLeft);
//	}
//
//
//	if(RecvLeft) {
//		if(dir == 0) {
//			for (int k = -rim; k <= gdata.mx[2]+rim; k++){
//				for (int j = -rim; j <= gdata.mx[1]+rim; j++){
//					for (int i = imaxRecv; i >= -rim; --i){
//						omb(i,j,k) = Recv(i,j,k);
//					}
//				}
//			}
//		} else if (dir == 1) {
//			for (int j = imaxRecv; j >= -rim; --j){
//				for (int k = -rim; k <= gdata.mx[2]+rim; k++){
//					for (int i = -rim; i <= gdata.mx[0]+rim; i++){
//						omb(i,j,k) = Recv(i,j,k);
//					}
//				}
//			}
//		} else if (dir == 2) {
//			for (int k = imaxRecv; k >= -rim; --k){
//				for (int j = -rim; j <= gdata.mx[1]+rim; j++){
//					for (int i = -rim; i <= gdata.mx[0]+rim; i++){
//						omb(i,j,k) = Recv(i,j,k);
//					}
//				}
//			}
//		}
//	}
//
//}

void gridFunc::ExchangePositions(Data &gdata, int top) {
	// First we need to find if we have a corresponding opposite rank
	// (should moslty be the case - only when deviding by odd number
	// would we expect problems)

	// cerr << " Exchanging " << gdata.rank << " " << gdata.coords[0] << endl;

#if (GEOM != CARTESIAN)
	PhiPos[top].resize(Index::set(-B),
	                   Index::set(gdata.global_mx[gdata.phi_dir]+B));

#ifdef parallel
	if(gdata.nproc[gdata.phi_dir] > 1) {
#if(GEOM == CYLINDRICAL)
		if(gdata.coords[0] == 0) {
#elif (GEOM == SPHERICAL)
		   if(top==0 &&gdata.coords[1]==0 || top==1 && gdata.coords[1]==gdata.nproc[1]-1) {
#endif

			// Prepare array holding nproc entries in phi-direction
			NumMatrix<REAL, 1> PhiPos_loc[gdata.nproc[gdata.phi_dir]];
			
#if(GEOM == CYLINDRICAL)
			int rankOwn = gdata.rank_constz;
#elif (GEOM == SPHERICAL)
			int rankOwn = gdata.rank_constxy;
#endif

			// Resize and clear local array
			PhiPos_loc[rankOwn].resize(Index::set(-B), 
			                           Index::set(gdata.mx[gdata.phi_dir]+B));
			PhiPos_loc[rankOwn].clear();

			// Fill array with phi-data
			for(int ip = -B; ip <= gdata.mx[gdata.phi_dir]+B; ++ip) { 
				double phi = gdata.getCen(gdata.phi_dir, ip);
			
				// PhiPos_loc[gdata.rank_constz](ip) = phi;
				// MyData(ip) = phi;
				PhiPos_loc[rankOwn](ip) = phi;
				
			}

 
			// Prepare mpi-communication:
			int size = (PhiPos_loc[rankOwn].getHigh(0) - 
			            PhiPos_loc[rankOwn].getLow(0) + 1);
		
			int req_num = 2*(gdata.nproc[gdata.phi_dir]-1);
			MPI_Request requests[req_num];
			MPI_Status statusrl[req_num];
			
			for (int i=0; i<req_num; i++){
				requests[i] = MPI_REQUEST_NULL;
			}
			
			// Loop over all processes with same z at origin:
			int number = 0;
			
			//		NumMatrix<REAL, 1> PhiPos_loc[gdata.nproc[gdata.phi_dir]];
// #if(GEOM == CYLINDRICAL)
// 			int rankOwn = gdata.rank_constz;
// #elif (GEOM == SPHERICAL)
// 			int rankOwn = gdata.rank_constx;
// #endif
			for(int irank=0; irank<gdata.nproc[gdata.phi_dir]; ++irank) {
			
				int rankother;
#if(GEOM == CYLINDRICAL)
				rankother = gdata.rankpos_z(gdata.coords[0],irank);
				MPI_Comm my_comm = gdata.comm_constz;
#elif (GEOM == SPHERICAL)
				rankother = gdata.rankpos_xy(irank);
				MPI_Comm my_comm = gdata.comm_constxy;
#endif

				// Treat data from other ranks
				if(rankOwn != rankother) {

					PhiPos_loc[irank].resize(Index::set(-B), 
					                         Index::set(gdata.mx[gdata.phi_dir]+B));
					PhiPos_loc[irank].clear();


					int tag = rankother; // message tag
					
					MPI_Irecv((double *)PhiPos_loc[irank], size, MPI_DOUBLE,
					          rankother, tag, my_comm, &requests[number]);
					number++;
					
					tag = gdata.rank_constxy;
					
					// Send data:
					MPI_Isend((double *)PhiPos_loc[rankOwn], size, 
					          MPI_DOUBLE, rankother, tag, my_comm, 
					          &requests[number]);
					number++;

				// } else {
				// 	// If this ever becomes a problem I have to store the local rank
				// 	// separately.
				// 	if(irank != rankOwn) {
				// 		cerr << " Mixup with local ranks along phi " << endl;
				// 		cerr << " Exiting " << endl;
				// 		exit(3);
				// 	}
				}
				// } else {
				// 	// Is this correct??

				// 	PhiPos_loc[irank].resize(Index::set(-B), 
				// 	                         Index::set(gdata.mx[gdata.phi_dir]+B));
				// 	PhiPos_loc[irank].clear();
				// 	for(int ip = -B; ip <= gdata.mx[gdata.phi_dir]+B; ++ip) { 
				// 		double phi = gdata.getCen(gdata.phi_dir, ip);
						
				// 		// PhiPos_loc[gdata.rank_constz](ip) = phi;
				// 		// MyData(ip) = phi;
				// 		PhiPos_loc[irank](ip) = phi;
				// 	}

				// 	//					PhiPos_loc[irank] = irank;
				// }
			
			}

			if(req_num > 0) {
				MPI_Waitall(req_num, requests, statusrl);
			}

			for(int irank=0; irank<gdata.nproc[gdata.phi_dir]; ++irank) {
				int shift = irank*gdata.get_RankWidth(gdata.phi_dir);
				for(int ip = -B; ip <= gdata.mx[gdata.phi_dir]+B; ++ip) { 
					PhiPos[top](ip+shift) = PhiPos_loc[irank](ip);
				}
			}
			
		   }			
			

		} else {
#endif

	// Compute Phi-positions also for serial run
	for(int ip = -B; ip <= gdata.mx[gdata.phi_dir]+B; ++ip) { 
		double phi = gdata.getCen(gdata.phi_dir, ip);
		
		PhiPos[top](ip) = phi;
		
	}
	

#ifdef parallel

	}

#endif

#endif // not cartesian

		}


#ifdef parallel

void gridFunc::get_AxisData(Data &gdata, 
                            NumMatrix<REAL,3> &myBound,
                            NumMatrix<REAL,3> &global_bcData,
                            int top) {
	// For parallel case: exchange data an build global array of
	// boundary data. For serial data - just pipe through?? No! That
	// would double the array data - argh So maybe avoid that
	// altogether. Maybe for the serial case something like & operator
	// in the external routine?

	// Array holding the global version of the data (full phi-extent)
	// NumMatrix<REAL,3> global_bcData;

#if (GEOM != CARTESIAN)

#if (GEOM == CYLINDRICAL) // phi in y-direction
	int lbound[DIM] = { 0, -B, -B};
	int ubound_global[DIM] = { 2, gdata.global_mx[1]+B, gdata.mx[2]+B };
	int ubound_local[DIM]  = { 2, gdata.mx[1]+B, gdata.mx[2]+B };

	if(myBound.getName() == "B_x") {
		lbound[0] = -1;
		ubound_global[0] = 1;
		ubound_local[0] = 1;
	}
#elif (GEOM == SPHERICAL) // phi in z-direction

	int lbound[DIM] = { -B, 0, -B};
	int ubound_global[DIM] = {gdata.mx[0]+B, 2, gdata.global_mx[2]+B };
	int ubound_local[DIM]  = {gdata.mx[0]+B, 2, gdata.mx[2]+B };

	if(myBound.getName() == "B_y") {
		lbound[1] = -1;
		ubound_global[1] = 1;
		ubound_local[1] = 1;
	}

	// At upper boundary:
	if(top==1) {
		int mxtet = gdata.mx[1];
		lbound[1] = mxtet-B+1;
		ubound_global[1] = mxtet;
		ubound_local[1] = mxtet;

		if(myBound.getName() == "B_y") {
			lbound[1] -= 1;
			ubound_global[1] -= 1;
			ubound_local[1] -= 1;
		} else if (myBound.getName() == "A_x" || myBound.getName() == "A_z") {
		  ubound_global[1] += 1;
		  ubound_local[1] += 1;
		}
	}
#endif

#if (GEOM != CARTESIAN)
	global_bcData.resize(lbound, ubound_global);
#endif

	// So far only for cylindrical coords
#if (GEOM == CYLINDRICAL)

	// cout << " My coords: " << gdata.coords[0] << " " << gdata.rank << " " << gdata.rank_constz << " " << gdata.rankpos(gdata.coords[0],1, gdata.coords[2]) << " " << gdata.coords[2] << endl;

	// MPI_Finalize();
	// exit(3);

	// Only applicable for axis-adjacent grid points
	if(gdata.coords[0] == 0) {

		int ownRank(gdata.rank_constz);

		// Assiging local boundary data array:
		NumMatrix<REAL,3> local_bcData[gdata.nproc[gdata.phi_dir]];
		local_bcData[ownRank].resize(lbound, ubound_local);

		for(int iz=lbound[2]; iz<=ubound_local[2]; ++iz) {
			for(int iy=lbound[1]; iy<=ubound_local[1]; ++iy) {
				for(int ix=lbound[0]; ix<=ubound_local[0]; ++ix) {
					local_bcData[ownRank](ix, iy, iz) = myBound(ix, iy, iz);
				}
			}
		}

		// Preparing MPI-communication
		int size = ((local_bcData[ownRank].getHigh(0) -
		             local_bcData[ownRank].getLow(0) + 1)*
		            (local_bcData[ownRank].getHigh(1) -
		             local_bcData[ownRank].getLow(1) + 1)*
		            (local_bcData[ownRank].getHigh(2) -
		             local_bcData[ownRank].getLow(2) + 1));
		
		// number of MPI communication events (one less, because no
		// self-communication)
		int req_num = 2*(gdata.nproc[gdata.phi_dir]-1);
		MPI_Request requests[req_num];
		MPI_Status statusrl[req_num];

		// Initialise requests
		for (int i=0; i<req_num; i++){
			requests[i] = MPI_REQUEST_NULL;
		}

		// Loop over all processes with same z at origin:
		int number = 0;

#if(GEOM == CYLINDRICAL)
		int rankOwn = gdata.rank_constz;
		MPI_Comm my_comm = gdata.comm_constz;
#elif (GEOM == SPHERICAL)
		int rankOwn = gdata.rank_constxy;
		MPI_Comm my_comm = gdata.comm_constxy;
#endif

		for(int irank=0; irank<gdata.nproc[gdata.phi_dir]; ++irank) {

			int rankother;
			// obtain position of rank irank
			if(gdata.phi_dir == 1) {
				rankother = gdata.rankpos_z(gdata.coords[0],irank);
			}

			// Get other data
			// Avoid self-communication

			// if(gdata.rank == 2) {
			// 	cerr << " The others: " << rankother << endl;
			// }

			if(rankOwn != rankother) {

				// Prepare array to hold other boundary data
				local_bcData[irank].resize(lbound, ubound_local);

				int tag = rankother; // message tag

				MPI_Irecv((double *)local_bcData[irank], size, MPI_DOUBLE,
				          rankother, tag, my_comm, &requests[number]);
				number++;

				tag = gdata.rank;

				// Send data:
				MPI_Isend((double *)local_bcData[ownRank], size, MPI_DOUBLE, 
				          rankother, tag, my_comm, &requests[number]);
				number++;

			}

		}

		MPI_Waitall(req_num, requests, statusrl);

		for(int irank=0; irank<gdata.nproc[gdata.phi_dir]; ++irank) {
			store_AxisData(gdata, lbound, ubound_local, irank,
			               local_bcData[irank], global_bcData);
		}
	
		// if(gdata.rank == 1 && myBound.getName() == "v_z") {
		// 	cout << " other rank " << 3 << endl;
		// 	for(int j = -B; j <= gdata.mx[1]+B; ++j) {
		// 		cout << j << " fu " << (local_bcData[0](1,j,2)-30)/pi << " " << (local_bcData[1](1,j,2)-30)/pi << " " << gdata.get_z(2) << " " << ownRank << endl;
		// 	}
			
		// }




	} // if coords[0] == 0
	// MPI_Finalize();
	// exit(3);
#else // Now spherical(?!)
	if((top==0 && gdata.coords[1] == 0 )||
	   (top==1 && gdata.coords[1]==gdata.nproc[1]-1)) {

		// Communicate with ranks with same x-values
		int ownRank(gdata.rank_constxy);

		// Assiging local boundary data array:
		NumMatrix<REAL,3> local_bcData[gdata.nproc[gdata.phi_dir]];
		local_bcData[ownRank].resize(lbound, ubound_local);

		for(int iz=lbound[2]; iz<=ubound_local[2]; ++iz) {
			for(int iy=lbound[1]; iy<=ubound_local[1]; ++iy) {
				for(int ix=lbound[0]; ix<=ubound_local[0]; ++ix) {
					local_bcData[ownRank](ix, iy, iz) = myBound(ix, iy, iz);
				}
			}
		}

		// Preparing MPI-communication
		int size = ((local_bcData[ownRank].getHigh(0) - local_bcData[ownRank].getLow(0) + 1)*
				(local_bcData[ownRank].getHigh(1) - local_bcData[ownRank].getLow(1) + 1)*
				(local_bcData[ownRank].getHigh(2) - local_bcData[ownRank].getLow(2) + 1));

		// number of MPI communication events (one less, because no
		// self-communication)
		int req_num = 2*(gdata.nproc[gdata.phi_dir]-1);

		if(req_num>0) {
			MPI_Request requests[req_num];
			MPI_Status statusrl[req_num];

			// Initialise requests
			for (int i=0; i<req_num; i++){
				requests[i] = MPI_REQUEST_NULL;
			}

			// Loop over all processes with same z at origin:
			int number = 0;

#if(GEOM == CYLINDRICAL)
			int rankOwn = gdata.rank_constz;
			MPI_Comm my_comm = gdata.comm_constz;
#elif (GEOM == SPHERICAL)
			int rankOwn = gdata.rank_constxy;
			MPI_Comm my_comm = gdata.comm_constxy;
#endif

			for(int irank=0; irank<gdata.nproc[gdata.phi_dir]; ++irank) {

				int rankother;
				// obtain position of rank irank
				// AT ME: does this have to be 2?
				if(gdata.phi_dir == 2) {
					rankother = gdata.rankpos_xy(irank);
				}

				// Get other data
				// Avoid self-communication

				if(rankOwn != rankother) {

					// Prepare array to hold other boundary data
					local_bcData[irank].resize(lbound, ubound_local);

					int tag = rankother; // message tag

					MPI_Irecv((double *)local_bcData[irank], size, MPI_DOUBLE,
							rankother, tag, my_comm, &requests[number]);
					number++;

					tag = gdata.rank_constxy;

					// Send data:
					MPI_Isend((double *)local_bcData[ownRank], size, MPI_DOUBLE,
							rankother, tag, my_comm, &requests[number]);
					number++;

				}

			}

			MPI_Waitall(req_num, requests, statusrl);

		}

		for(int irank=0; irank<gdata.nproc[gdata.phi_dir]; ++irank) {
			store_AxisData(gdata, lbound, ubound_local, irank,
					local_bcData[irank], global_bcData);
		}
	}

	// exit(3);
#endif
#endif // not Cartesian

}


void gridFunc::store_AxisData(Data &gdata, int lbound[DIM],
                              int ubound_local[DIM], int rank_local,
                              NumMatrix<REAL,3> &local_bcData,
                              NumMatrix<REAL,3> &global_bcData) {

	// Use data by other rank (or own) to store in global data array

#if (GEOM == CYLINDRICAL) // phi in y-direction
	
	//	int shift = rank_local*gdata.mx[1];
	int shift = rank_local*gdata.get_RankWidth(1);

	// set data into global array
	for(int iz=lbound[2]; iz<=ubound_local[2]; ++iz) {
		for(int iy=lbound[1]; iy<=ubound_local[1]; ++iy) {
			for(int ix=lbound[0]; ix<=ubound_local[0]; ++ix) {
				global_bcData(ix, iy+shift, iz) = local_bcData(ix, iy, iz);
			}
		}
	}

#elif (GEOM == SPHERICAL) // phi in z-direction

	int shift = rank_local*gdata.get_RankWidth(2);

	// set data into global array
	for(int iz=lbound[2]; iz<=ubound_local[2]; ++iz) {
		for(int iy=lbound[1]; iy<=ubound_local[1]; ++iy) {
			for(int ix=lbound[0]; ix<=ubound_local[0]; ++ix) {
				global_bcData(ix, iy, iz+shift) = local_bcData(ix, iy, iz);
			}
		}
	}

#endif

#if (FALSE)
	// General Version:
	// Axis for cyl coords
	if(gdata.phi_dir == 1) { // for cyl coords
		// Compute shift in phi-direction (depends on local rank)
		// in this case rank_local == rank_constz

		int shift[3] = {0,0,0};

		shift[gdata.phi_dir] = rank_local*gdata.mx[gdata.phi_dir];

		// set data into global array
		for(int iz=-B; iz<=gdata.mx[2]+B; ++iz) {
			for(int iy=-B; iy<=gdata.mx[1]+B; ++iy) {
				for(int ix=-B; ix<=gdata.mx[0]+B; ++ix) {
					global_bcData(ix+shift[0], iy+shift[1], iz+shift[2]) = 
						local_bcData(ix, iy, iz);
				}
			}
		}

	}
#endif

}


void gridFunc::ExchangeData(Data &gdata) {
	// First we need to find if we have a corresponding opposite rank
	// (should moslty be the case - only when deviding by odd number
	// would we expect problems)

}
#endif


void gridFunc::compute_AxisPartners(Data &gdata, int top) {
	// Computing the partners for axis boundary conditions - so far
	// only for cylindrical coordinates

#if (GEOM != CARTESIAN)
	ExchangePositions(gdata, top);

	// Array over full phi-domain
#ifdef parallel
	AxisPartnersPhi[top].resize(Index::set(-B),
			Index::set(gdata.global_mx[gdata.phi_dir]+B));
#else
	AxisPartnersPhi[top].resize(Index::set(-B),
			Index::set(gdata.mx[gdata.phi_dir]+B));
#endif

	double Pi2 = 2.*M_PI;
	double eps = (1.e-5)*(gdata.xe[gdata.phi_dir]-
	                      gdata.xb[gdata.phi_dir])/gdata.global_mx[gdata.phi_dir];

	// loop over azimuthal domain:
	// look for mapping p.s
	for (int ip = -B; ip <= gdata.global_mx[gdata.phi_dir]+B; ++ip) {

		AxisPartnersPhi[top](ip) = -100;   // default value to flag a failure
		// 'phi' = required partner position

		double phi = PhiPos[top](ip) + M_PI;
		phi -= Pi2*static_cast<int>(phi/Pi2);    // keep phi in [0,2pi]

		// if (phi >= gdata.xb[gdata.phi_dir] and   // desired mp inside grid
		//     phi <= gdata.xe[gdata.phi_dir]) {   //  => loop candidates
		if (phi >= gdata.global_xb[gdata.phi_dir] and  // desired mp inside grid
		    phi <= gdata.global_xe[gdata.phi_dir]) {   //  => loop candidates

			// loop candidates
			for (int imp = 0; imp <= gdata.global_mx[gdata.phi_dir]; ++imp) { 
				//if (abs(gdata.getCen(gdata.phi_dir, imp) - phi) < eps) {
				if (abs(PhiPos[top](imp) - phi) < eps) {
					AxisPartnersPhi[top](ip) = imp; // found one
					break;
				}
			}

			if (AxisPartnersPhi[top](ip) == -100) {  // mp needed, but none found
				cerr << "error: singularity present, "
				     << "but cells cannot be mapped." << endl;
				exit(1);
			}
		} else {
			AxisPartnersPhi[top](ip) = -10;  // this means "no mp needed"
		}


		// Workaround for single grid point in phi-direction:
		if(gdata.global_mx[gdata.phi_dir]==0) {
			AxisPartnersPhi[top](ip) = 0;
		}
//		 if(gdata.rank==0) {
//		 	cout << " My partner: " << ip << " " << AxisPartnersPhi[top](ip) << " ";
//		 	cout << phi << " " << (phi-M_PI)/M_PI;
//		 	cout << endl;
//		 }
	}
#endif	

	// Additional convenience functions holding sin and cos of positions
	// along the grid
	for(unsigned int dir=0; dir<3; ++dir) {
		cosPos[dir].resize(Index::set(-B), Index::set(gdata.mx[dir]+B));
		sinPos[dir].resize(Index::set(-B), Index::set(gdata.mx[dir]+B));
		// Fill values of arrays
		for(int iPos=-B; iPos<=gdata.mx[dir]+B; ++iPos) {
			sinPos[dir](iPos) = sin(gdata.get_pos(dir, iPos, 0));
			cosPos[dir](iPos) = cos(gdata.get_pos(dir, iPos, 0));
		}
	}

}



int gridFunc::get_AxisPartners(Data &gdata, int top, int iz) {
	int shift = 0;
#ifdef parallel
	int dir = 0;
#if (GEOM == CYLINDRICAL)
	dir = 1;
#elif (GEOM == SPHERICAL)
	dir = 2;
#endif
	shift = gdata.get_RankShift(dir);
#endif

	int i_phi_other = AxisPartnersPhi[top](iz+shift);
	return i_phi_other;
}




void gridFunc::bc_Axis(Data &gdata, ProblemType &Problem,
		NumMatrix<REAL,3> &omb, int dir, int top, int q, int rim) {
#if (GEOM == CYLINDRICAL)
	if(dir==0) {
		bc_AxisCyl(gdata, Problem, omb, q, rim);
	}
#elif (GEOM == 3)
	if(dir==1) {
		bc_AxisSph(gdata, Problem, omb, top, q, rim);
	}
#endif
}


void gridFunc::bc_AxisCyl(Data &gdata, ProblemType &Problem,
		NumMatrix<REAL,3> &omb, int q, int rim) {

	// Special treatment of coordinate axis - so far only for
	// cylindrical coordinates.

#if (GEOM == CYLINDRICAL)
	// Check for necessity of sign inversion at the origin:
	string qname = omb.getName(); 
	string subsp = qname.substr(qname.size()-2,2); // last two letters
	short sign = (subsp == "_x" || subsp == "_y") ? -1 : +1; // sign inv.

	// if(q == 7) {
	// 	cout << " Partner: " << AxisPartnersPhi(5) << " ";
	// 	cout << (gdata.mx[1]+1)/2 + 5 << endl;
	// }

#ifdef parallel
	// Parallel version -- need to work with extra global array
	NumMatrix<REAL,3> omb_global; // additional global array
	get_AxisData(gdata, omb, omb_global);
#else
	NumMatrix<REAL,3> & omb_global = omb;
#endif

	for (int k = -rim; k <= gdata.mx[2]+rim; ++k) {
		for (int j = -rim; j <= gdata.mx[1]+rim; ++j) {
//			if (AxisPartnersPhi[0](j) != -10) {
			if (get_AxisPartners(gdata,0,j) != -10) {
				for (int i = -rim; i <= -1; i++) {
					if (qname == "B_x") {
						omb(i,j,k) = -1*omb_global(-2-i,get_AxisPartners(gdata,0,j),k);
//						omb(i,j,k) = -1*omb_global(-2-i,AxisPartnersPhi[0](j),k);
					} else {
						omb(i,j,k) = sign*omb_global(-1-i,get_AxisPartners(gdata,0,j),k);
//						omb(i,j,k) = sign*omb_global(-1-i,AxisPartnersPhi[0](j),k);
					}

				}
			}
		}
	}

#endif

}


void gridFunc::bc_AxisSph(Data &gdata, ProblemType &Problem,
		NumMatrix<REAL,3> &omb, int top, int q, int rim) {
#if (GEOM == 3)

#ifdef parallel
	// Parallel version -- need to work with extra global array
	NumMatrix<REAL,3> omb_global; // additional global array
	get_AxisData(gdata, omb, omb_global, top);
#else
	NumMatrix<REAL,3> & omb_global = omb;
#endif

	string qname = omb.getName();
	string subsp = qname.substr(qname.size()-2,2); // last two letters

	const int mxtet = gdata.mx[1];
	short si = (subsp == "_y" or subsp == "_z") ? -1 : +1;
	for (int iz = -rim; iz <= gdata.mx[2]+rim; ++iz) {
		if (get_AxisPartners(gdata,top,iz) != -10) {
//		if (AxisPartnersPhi[top](iz) != -10) {
			for (int ix = -rim; ix <= gdata.mx[0]+rim; ++ix) {
				if (top == 0 and gdata.get_singularity_treatment(2)==1) {
					for (int iy = -1; iy >= -rim; --iy) {
						if (qname == "B_y") {
							omb(ix,iy,iz) = -1*omb_global(ix,-iy-2,get_AxisPartners(gdata, top, iz));
//							omb(ix,iy,iz) = -1*omb_global(ix,-iy-2,AxisPartnersPhi[top](iz));
						} else if (qname == "A_x") { // means A_r
						  if(iy>-3) {
							omb(ix,iy,iz) = omb_global(ix,-iy,get_AxisPartners(gdata, top, iz));
//							omb(ix,iy,iz) = omb_global(ix,-iy,AxisPartnersPhi[top](iz));
						  }
						} else if (qname == "A_z") { // means A_phi
						  if(iy>-3) {
							omb(ix,iy,iz) = omb_global(ix,-iy,get_AxisPartners(gdata, top, iz));
//							omb(ix,iy,iz) = omb_global(ix,-iy,AxisPartnersPhi[top](iz));
//							omb(ix,iy,iz) = -omb_global(ix,-iy,AxisPartnersPhi[top](iz));
						  }
						} else {
							omb(ix,iy,iz) = si*omb_global(ix,-iy-1,get_AxisPartners(gdata, top, iz));
//							omb(ix,iy,iz) = si*omb_global(ix,-iy-1,AxisPartnersPhi[top](iz));
						}
					}
				} else if (top == 1 and gdata.get_singularity_treatment(3)==1) {
					for (int iy = 1; iy <= rim; ++iy) {
						if (qname == "B_y") {
							omb(ix,mxtet+iy,iz) = -1*omb_global(ix,mxtet-iy,get_AxisPartners(gdata, top, iz));
//							omb(ix,mxtet+iy,iz) = -1*omb_global(ix,mxtet-iy,AxisPartnersPhi[top](iz));
						} else if (qname == "A_x") { // means A_r
						  omb(ix,mxtet+iy,iz) = omb_global(ix,mxtet+2-iy,get_AxisPartners(gdata, top, iz));
//						  omb(ix,mxtet+iy,iz) = omb_global(ix,mxtet+2-iy, AxisPartnersPhi[top](iz));
//						  omb(ix,mxtet+iy,iz) = -omb_global(ix,mxtet+2-iy,AxisPartnersPhi[top](iz));
						} else if (qname == "A_z") { // means A_phi
							omb(ix,mxtet+iy,iz) = omb_global(ix,mxtet+2-iy,get_AxisPartners(gdata, top, iz));
//							omb(ix,mxtet+iy,iz) = omb_global(ix,mxtet+2-iy,AxisPartnersPhi[top](iz));
//							omb(ix,mxtet+iy,iz) = -1*omb_global(ix,mxtet+2-iy,AxisPartnersPhi[top](iz));
						} else {
							omb(ix,mxtet+iy,iz) = si*omb_global(ix,mxtet+1-iy,get_AxisPartners(gdata, top, iz));
//							omb(ix,mxtet+iy,iz) = si*omb_global(ix,mxtet+1-iy,AxisPartnersPhi[top](iz));
						}
					}
				}
			} // ix loop
		}
	}
#endif

}


void gridFunc::do_AxisValCorrectionCyl(Data &gdata, ProblemType &Problem) {
#if (GEOM == CYLINDRICAL)
	int jnum = gdata.get_RankWidth(1);
	for (int iz = -B; iz <= gdata.mx[2]+B; ++iz) {
		REAL Bcx_sum = 0.;  // new values for Cartesian Bx, By
		REAL Bcy_sum = 0.;  //  on the singular axis
		int N_terms(0);
			//		for (int j = 0; j <= gdata.mx[1]; j++) {
			// Do not use phi=0 twice!
		if(gdata.get_singularity_treatment(0) == 1) { // at lower r-bound
			for (int iy = 0; iy < jnum; iy++) {
				// br0_j, bp0_j: projected to (r,phi)=(0,phi_j) via avg.
				REAL br0_j = (+gdata.om[q_Bx]( 0,iy,  iz)
				              +gdata.om[q_Bx](-2,iy,  iz))/2.;
				REAL bp0_j = (+gdata.om[q_By]( 0,iy-1,iz)
				              +gdata.om[q_By]( 0,iy  ,iz)
				              +gdata.om[q_By](-1,iy-1,iz)
				              +gdata.om[q_By](-1,iy  ,iz))/4.;
				REAL phi = gdata.getCen(1, iy);
				Bcx_sum += br0_j * cos(phi) - bp0_j * sin(phi);
				Bcy_sum += br0_j * sin(phi) + bp0_j * cos(phi);
				// if(gdata.rank==0 && k==1) {
				// 	cout << " axis: " << j << " " << br0_j << " " << bp0_j << " " << Bcx_sum << " " << Bcy_sum << endl;
				// 	cout << "    " << gdata.om[q_By]( 0,j-1,k) << " ";
				// 	cout << gdata.om[q_By]( 0,j  ,k) << " ";
				// 	cout << gdata.om[q_By](-1,j  ,k) << " ";
				// 	cout << gdata.om[q_By](-1,j-1,k) << " ";
				// 	cout << endl;
				// }
				N_terms++;
			}

		} else {
			Bcx_sum = 0.;
			Bcy_sum = 0.;
			N_terms = 0.;
		}
#ifdef parallel
		gdata.MpiSum(2, Bcx_sum); // Sum only for ranks with same z
		gdata.MpiSum(2, Bcy_sum); // Sum only for ranks with same z
		gdata.MpiSum(2, N_terms); // Sum only for ranks with same z
#endif
		REAL Bcx_ave(Bcx_sum/N_terms);
		REAL Bcy_ave(Bcy_sum/N_terms);

		// if((gdata.rank==0 && k==gdata.mx[2]) ||
		//    (gdata.rank==1 && k==-1)) {
		// 	cout << setiosflags(ios::scientific) << setprecision(12);
		// 	cout << " corr val: " << gdata.rank << " " << Bcx_ave << " " <<	Bcy_ave << endl;
		// 	cout << resetiosflags(ios::scientific) << setprecision(6);
		// }

		if(gdata.get_singularity_treatment(0) == 1) { // at lower r-bound
			for (int iy = -B; iy <= gdata.mx[1]+B; iy++) { // set Br(0)s
				REAL phi = gdata.getCen(1, iy);
				REAL Br0_ave = ( Bcx_ave * cos(phi) +
				                 Bcy_ave * sin(phi) );
				gdata.om[q_Bx](-1,iy,iz) = Br0_ave;
			} // j
		} // check if r==0
	} // k
	// exit(3);


#endif
}
void gridFunc::do_AxisValCorrectionSph(Data &gdata, ProblemType &Problem, bool top) {
	//! Correction for mag field values on polar axis

//#if (GEOM == SPHERICAL)
	// Check if singularity is on current rank
	if(gdata.get_singularity_treatment(2+top) == 1) {
		REAL cost0 = 1.-2.*top;      // cos({0|Pi}) = {+1|-1}
		int jS = top*(gdata.mx[1]+1)-1; // {-1|mxtet} (sing. index)
		for (int ix = -B; ix <= gdata.mx[0]+B; ++ix) {
			REAL Bcx_sum = 0.;  // new values for Cartesian Bx, By
			REAL Bcy_sum = 0.;  //  on the singular axis
			for (int iz = 0; iz <= gdata.mx[2]; ++iz) {
				// bt0_k, bp0_k: projected to (tet,phi)=(0,phi_k) via avg.
				REAL bt0_k = (+gdata.om[q_By](ix,jS-1, iz  )
						+gdata.om[q_By](ix,jS+1, iz  ))*0.5;
				REAL bp0_k = (+gdata.om[q_Bz](ix,jS  , iz-1)
						+gdata.om[q_Bz](ix,jS  , iz  )
						+gdata.om[q_Bz](ix,jS+1, iz-1)
						+gdata.om[q_Bz](ix,jS+1, iz  ))*0.25;

				Bcx_sum += bt0_k * cost0 * cosPos[2](iz) - bp0_k * sinPos[2](iz);
				Bcy_sum += bt0_k * cost0 * sinPos[2](iz) + bp0_k * cosPos[2](iz);
			}
			int N_terms = gdata.mx[2]+1;

#ifdef parallel
			// ToDo: Communicate [Bc_sum, N_terms] as above
			gdata.MpiSum(0, Bcx_sum); // Sum only for ranks with same x
			gdata.MpiSum(0, Bcy_sum); // Sum only for ranks with same x
			gdata.MpiSum(0, N_terms); // Sum only for ranks with same x

#endif // [parallel]
			REAL Bcx_ave(Bcx_sum/N_terms);
			REAL Bcy_ave(Bcy_sum/N_terms);

			for (int iz = -B; iz <= gdata.mx[2]+B; iz++) { // set Bt(0)s
				REAL Bt0_ave = ( Bcx_ave * cost0 * cosPos[2](iz) +
						Bcy_ave * cost0 * sinPos[2](iz) );
				//         NB:  -Bcz_sum * sint0 ~ sin({0|Pi}) = 0
				gdata.om[q_By](ix,jS,iz) = Bt0_ave;
				// gdata.om[q_By](i,jS,k) = (gdata.om[q_By](i,jS-1,k)+
				//                           gdata.om[q_By](i,jS+1,k))/2.;
			} // k
		} // i
	} else { // Other ranks have to communicate, too
#ifdef parallel
		double dummy(0.);
		int dummy_int(0);
		for (int ix = -B; ix <= gdata.mx[0]+B; ++ix) {
			dummy = 0.;
			gdata.MpiSum(0, dummy); // Sum only for ranks with same x
			dummy = 0.;
			gdata.MpiSum(0, dummy); // Sum only for ranks with same x
			dummy_int = 0;
			gdata.MpiSum(0, dummy_int); // Sum only for ranks with same x
		}
#endif
	}

//#endif
}



void gridFunc::bc_Extrapolate(Data &gdata, ProblemType &Problem, 
                              NumMatrix<REAL,3> &omb,
                              int dir, int above, int q, int rim)
{

	string qname = omb.getName();

	if(dir == 0) {
		if(above == 0) {
			int imin = -1;
//			if(q == Problem.q_Bx) {
			if(qname=="B_x") {
				imin = -2;
			}
			for (int k = -rim; k <= gdata.mx[2]+rim; k++){
				for (int j = -rim; j <= gdata.mx[1]+rim; ++j){
					for (int i = imin; i >= -rim; --i){
						omb(i,j,k) = omb(i+1,j,k);    // x_min
						// omb(i,j,k) = 2.*omb(i+1,j,k) - omb(i+2,j,k);  // x_min
					}
				}
			}
		} else if (above == 1) {
			int imax = 1;
			if(omb.getName() == "A_y" || omb.getName() == "A_z") {
				imax = 2;
			}
			for (int k = -rim; k <= gdata.mx[2]+rim; k++){
				for (int j = -rim; j <= gdata.mx[1]+rim; j++){
					for (int i = imax; i <= +rim; i++){
						omb(gdata.mx[0]+i,j,k) = omb(gdata.mx[0]+i-1,j,k); //  x_max
						// omb(gdata.mx[0]+i,j,k) = 2.*omb(gdata.mx[0]+i-1,j,k) - omb(gdata.mx[0]+i-2,j,k); //  x_max
					}
				}
			}
		}
	}
	
	// if(q == 7) {
	// 	cout << " Dir: " << dir << " " << above << " " << omb(0,0,0) << endl;
	// }

	if(dir == 1) {
		if(above == 0) {
			int jmin = -1;
//			if(q == Problem.q_By) {
			if(qname=="B_y") {
				jmin = -2;
			}
			for (int k = -rim; k <= gdata.mx[2]+rim; k++){
				for (int i = -rim; i <= gdata.mx[0]+rim; i++){
					for (int j = jmin; j >= -rim; --j){
						omb(i,j,k) = omb(i,j+1,k);    // y_min
						// omb(i,j,k) = 2.*omb(i,j+1,k) - omb(i,j+2,k);    // y_min
						// if(i==396 && k==0) {
						// 	cout << endl << " ext: ";
						// 	cout << j << " ";
						// 	cout << omb(i,j,k) << " ";
						// 	cout << omb(i,j+1,k) << " ";
						// 	cout << omb(i,j+2,k) << endl;
						// }
					}
				}
			}
		} else if (above == 1) {
			int jmax = 1;
			if(omb.getName() == "A_x" || omb.getName() == "A_z") {
				jmax = 2;
			}
			for (int k = -rim; k <= gdata.mx[2]+rim; k++){
				for (int i = -rim; i <= gdata.mx[0]+rim; i++){
					for (int j = jmax; j <= rim; j++){
						omb(i,gdata.mx[1]+j,k) = omb(i,gdata.mx[1]+j-1,k); //  y_max
						// omb(i,gdata.mx[1]+j,k) = 2.*omb(i,gdata.mx[1]+j-1,k) - omb(i,gdata.mx[1]+j-2,k); //  y_max
					}
				}
			}
		}
	}


	if(dir == 2) {
		if(above == 0) {
			int kmin = -1;
//			if(q == Problem.q_Bz) {
			if(qname=="B_z") {
				kmin = -2;
			}
			for (int j = -rim; j <= gdata.mx[1]+rim; j++){
				for (int i = -rim; i <= gdata.mx[0]+rim; i++){
					for (int k = kmin; k >= -rim; --k){
						omb(i,j,k) = omb(i,j,k+1);     // z_min
						// omb(i,j,k) = 2.*omb(i,j,k+1) - omb(i,j,k+2);     // z_min
					}
				}
			}
		} else if (above == 1) {
			int kmax = 1;
			if(omb.getName() == "A_x" || omb.getName() == "A_y") {
				kmax = 2;
			}
			for (int j = -rim; j <= gdata.mx[1]+rim; j++){
				for (int i = -rim; i <= gdata.mx[0]+rim; i++){
					for (int k = kmax; k <= rim; k++){
						omb(i,j,gdata.mx[2]+k) = omb(i,j,gdata.mx[2]+k-1);    // z_max
						// omb(i,j,gdata.mx[2]+k) = 2.*omb(i,j,gdata.mx[2]+k-1) - omb(i,j,gdata.mx[2]+k-2);    // z_max
					}
				}
			}
		}
	}
}




void gridFunc::bc_Outflow(Data &gdata, ProblemType &pr,
                          NumMatrix<REAL,3> &omb,
                          int dir, int above, int q, int rim)
{

	string qname = omb.getName();
	bool is_mag = false;
	if(qname=="B_x" || qname=="B_y" || qname=="B_z") {
		is_mag = true;
	}
	bool is_sx = false;
	if(qname=="s_x") {
		is_sx = true;
	}
	bool is_sy = false;
	if(qname=="s_y") {
		is_sy = true;
	}
	bool is_sz = false;
	if(qname=="s_z") {
		is_sz = true;
	}

	if(dir == 0) {
		if(above == 0) {
			for (int k = -rim; k <= gdata.mx[2]+rim; k++){
				for (int j = -rim; j <= gdata.mx[1]+rim; ++j){
					for (int i = -1; i >= -rim; --i){
						if(gdata.om[pr.q_sx](0,j,k) < 0.) {
//							if(q >= pr.q_Bx && q <= pr.q_Bz) {
							if(is_mag) {
								omb(i,j,k) = (3.*omb(i+1,j,k) -
								              3.*omb(i+2,j,k) +
								              omb(i+3,j,k));
							} else {
								omb(i,j,k) = omb(i+1,j,k);  // x_min
							}
						} else {
//							if(q == pr.q_sx) {
							if(is_sx) {
								omb(i,j,k) = -omb(-i-1,j,k);  // x_min
							} else {
								omb(i,j,k) =  omb(-i-1,j,k);  // x_min
							}
						}
					}
				}
			}
		} else if (above == 1) {
			for (int k = -rim; k <= gdata.mx[2]+rim; k++){
				for (int j = -rim; j <= gdata.mx[1]+rim; j++){
					for (int i = 1; i <= +rim; i++){
						if(gdata.om[pr.q_sx](gdata.mx[0],j,k) > 0.) {
//							if(q >= pr.q_Bx && q <= pr.q_Bz) {
							if(is_mag) {
								omb(gdata.mx[0]+i,j,k) = (3.*omb(gdata.mx[0]+i-1,j,k) -
								                          3.*omb(gdata.mx[0]+i-2,j,k) +
								                          omb(gdata.mx[0]+i-3,j,k));
							} else {
								omb(gdata.mx[0]+i,j,k) = omb(gdata.mx[0]+i-1,j,k); //  x_max
							}
						} else {
//							if(q == pr.q_sx) {
							if(is_sx) {
								omb(gdata.mx[0]+i,j,k) = -omb(gdata.mx[0]-i+1,j,k);
							} else {
								omb(gdata.mx[0]+i,j,k) =  omb(gdata.mx[0]-i+1,j,k);
							}
						}
					}
				}
			}
		}
	}


	if(dir == 1) {
		if(above == 0) {
			for (int k = -rim; k <= gdata.mx[2]+rim; k++){
				for (int i = -rim; i <= gdata.mx[0]+rim; i++){
					for (int j = -1; j >= -rim; --j){
						if(gdata.om[pr.q_sy](i,0,k) < 0.) {
//							if(q >= pr.q_Bx && q <= pr.q_Bz) {
							if(is_mag) {
								omb(i,j,k) = (3.*omb(i,j+1,k) -
								              3.*omb(i,j+2,k) +
								              omb(i,j+3,k));
							} else {
								omb(i,j,k) = omb(i,j+1,k);    // y_min
							}
						} else {
//							if(q == pr.q_sy) {
							if(is_sy) {
								omb(i,j,k) = -omb(i,-j-1,k);
							} else {
								omb(i,j,k) =  omb(i,-j-1,k);
							}
						}
					}
				}
			}
		} else if (above == 1) {
			for (int k = -rim; k <= gdata.mx[2]+rim; k++){
				for (int i = -rim; i <= gdata.mx[0]+rim; i++){
					for (int j = 1; j <= rim; j++){
						if(gdata.om[pr.q_sy](i,gdata.mx[1],k) > 0.) {
//							if(q >= pr.q_Bx && q <= pr.q_Bz) {
							if(is_mag) {
								omb(i,gdata.mx[1]+j,k) = (3.*omb(i,gdata.mx[1]+j-1,k) -
								                          3.*omb(i,gdata.mx[1]+j-2,k) +
								                          omb(i,gdata.mx[1]+j-3,k));
							} else {
								omb(i,gdata.mx[1]+j,k) = omb(i,gdata.mx[1]+j-1,k); //  y_max
							}
						} else {
//							if(q == pr.q_sy) {
							if(is_sy) {
								omb(i,gdata.mx[1]+j,k) = -omb(i,gdata.mx[1]-j+1,k);
							} else {
								omb(i,gdata.mx[1]+j,k) =  omb(i,gdata.mx[1]-j+1,k);
							}
						}
					}
				}
			}
		}
	}


	if(dir == 2) {
		int kmin = -1;
//		if(q == pr.q_Bz) {
		if(qname=="B_x") {
			kmin = -2;
		}
		if(above == 0) {
			for (int j = -rim; j <= gdata.mx[1]+rim; j++){
				for (int i = -rim; i <= gdata.mx[0]+rim; i++){
					for (int k = kmin; k >= -rim; --k){
						if(gdata.om[pr.q_sz](i,j,0) < 0.) {
//							if(q >= pr.q_Bx && q <= pr.q_Bz) {
							if(is_mag) {
								omb(i,j,k) = (3.*omb(i,j,k+1) -
								              3.*omb(i,j,k+2) +
								              omb(i,j,k+3));
							} else {
								omb(i,j,k) = omb(i,j,k+1);     // z_min
							}
						} else {
//							if(q == pr.q_sz) {
							if(is_sz) {
								omb(i,j,k) = -omb(i,j,-k-1);
							} else {
								omb(i,j,k) = omb(i,j,-k-1);
							}
						}
					}
				}
			}
			//       }
		} else if (above == 1) {
			for (int j = -rim; j <= gdata.mx[1]+rim; j++){
				for (int i = -rim; i <= gdata.mx[0]+rim; i++){
					for (int k = 1; k <= rim; k++){
						if(gdata.om[pr.q_sz](i,j,gdata.mx[2]) > 0.) {
//							if(q >= pr.q_Bx && q <= pr.q_Bz) {
							if(is_mag) {
								omb(i,j,gdata.mx[2]+k) = (3.*omb(i,j,gdata.mx[2]+k-1) -
								                          3.*omb(i,j,gdata.mx[2]+k-2) +
								                          omb(i,j,gdata.mx[2]+k-3));
							} else {
								omb(i,j,gdata.mx[2]+k) = omb(i,j,gdata.mx[2]+k-1);    // z_max
							}
						} else {
//							if(q == pr.q_sz) {
							if(is_sz) {
								omb(i,j,gdata.mx[2]+k) = -omb(i,j,gdata.mx[2]-k+1);
							} else {
								omb(i,j,gdata.mx[2]+k) = omb(i,j,gdata.mx[2]-k+1);
							}
						}
					}
				}
			}
		}
	}
}




void gridFunc::bc_Reflecting(Data &gdata, ProblemType &pr,
                             NumMatrix<REAL,3> &omb,
                             int dir, int above, int q, int rim)
{
	/*
	  Reflecting boundaries: inversion of velocity at the
	  boundaries. All other quantities are extrapolated via zero
	  gradient.
	*/

	string qname = omb.getName();
	bool is_sx = false;
	if(qname=="s_x") {
		is_sx = true;
	}
	bool is_sy = false;
	if(qname=="s_y") {
		is_sy = true;
	}
	bool is_sz = false;
	if(qname=="s_z") {
		is_sz = true;
	}

	if(dir == 0) {
		if(above == 0) {
			int imin = -1;
//			if(q == pr.q_Bx) {
			if(qname=="B_x") {
				imin = -2;
			}
			for (int k = -rim; k <= gdata.mx[2]+rim; k++){
				for (int j = -rim; j <= gdata.mx[1]+rim; ++j){
					for (int i = imin; i >= -rim; --i){
//						if(q == pr.q_sx) {
						if(is_sx) {
							omb(i,j,k) = -omb(-i-1,j,k);  // x_min
						} else {
							omb(i,j,k) = omb(i+1,j,k);  // x_min
						}
					}
				}
			}
		} else if (above == 1) {
			for (int k = -rim; k <= gdata.mx[2]+rim; k++){
				for (int j = -rim; j <= gdata.mx[1]+rim; j++){
					for (int i = 1; i <= +rim; i++){
//						if(q == pr.q_sx) {
						if(is_sx) {
							omb(gdata.mx[0]+i,j,k) = -omb(gdata.mx[0]-i+1,j,k);
						} else {
							omb(gdata.mx[0]+i,j,k) =  omb(gdata.mx[0]+i-1,j,k);
						}
					}
				}
			}
		}
	}


	if(dir == 1) {
		if(above == 0) {
			int jmin = -1;
//			if(q == pr.q_By) {
			if(qname=="B_y") {
				jmin = -2;
			}
			for (int k = -rim; k <= gdata.mx[2]+rim; k++){
				for (int i = -rim; i <= gdata.mx[0]+rim; i++){
					for (int j = jmin; j >= -rim; --j){
//						if(q == pr.q_sy) {
						if(is_sy) {
							omb(i,j,k) = -omb(i,-j-1,k);
						} else {
							omb(i,j,k) =  omb(i,j+1,k);
						}
					}
				}
			}
		} else if (above == 1) {
			for (int k = -rim; k <= gdata.mx[2]+rim; k++){
				for (int i = -rim; i <= gdata.mx[0]+rim; i++){
					for (int j = 1; j <= rim; j++){
//						if(q == pr.q_sy) {
						if(is_sy) {
							omb(i,gdata.mx[1]+j,k) = -omb(i,gdata.mx[1]-j+1,k);
						} else {
							omb(i,gdata.mx[1]+j,k) =  omb(i,gdata.mx[1]+j-1,k);
						}
					}
				}
			}
		}
	}


	if(dir == 2) {
		int kmin = -1;
//		if(q == pr.q_Bz) {
		if(qname=="B_z") {
			kmin = -2;
		}
		if(above == 0) {
			for (int j = -rim; j <= gdata.mx[1]+rim; j++){
				for (int i = -rim; i <= gdata.mx[0]+rim; i++){
					for (int k = kmin; k >= -rim; --k){
//						if(q == pr.q_sz) {
						if(is_sz) {
							omb(i,j,k) = -omb(i,j,-k-1);
						} else {
							omb(i,j,k) = omb(i,j,k-1);
						}
					}
				}
			}
			//       }
		} else if (above == 1) {
			for (int j = -rim; j <= gdata.mx[1]+rim; j++){
				for (int i = -rim; i <= gdata.mx[0]+rim; i++){
					for (int k = 1; k <= rim; k++){
//						if(q == pr.q_sz) {
						if(is_sz) {
							omb(i,j,gdata.mx[2]+k) = -omb(i,j,gdata.mx[2]-k+1);
						} else {
							omb(i,j,gdata.mx[2]+k) = omb(i,j,gdata.mx[2]+k-1);
						}
					}
				}
			}
		}
	}
}



void gridFunc::bc_Fixed(Data &gdata, ProblemType &pr,
                        NumMatrix<REAL,3> &omb,
                        int dir, int above, int q, int rim)
{
	/*
	  Fixed type boundary conditions. The initial values are stored
	  into an array and they are just copied into the boundaries for
	  every call of this bc.
	 */

	// Still missing some details: destinction between generic and
	// non-generic vars

	bool generic_vars = true;
#if (OMS_USER == TRUE)
	if(gdata.is_userField(omb)) {
		generic_vars = false;
	}
#endif

	if(dir == 0) {
		if(above == 0) {

			int imin[3] = {-rim, -rim, -rim};
			int imax[3] = {-1, gdata.mx[1]+rim, gdata.mx[2]+rim};

			// Do not overwrite B on face between computing domain and
			// ghost cells.
			if(q==pr.q_Bx) {
				imax[0] = -2;
			}

			if(generic_vars) {
				NumMatrix<REAL,3> & bcVals = bcVals_xb[q];
				// Doing bcs:
				bc_Fixed_general(omb, bcVals, imin, imax);
#if (OMS_USER == TRUE)
			} else {
				NumMatrix<REAL,3> & bcVals = bcVals_User_xb[q];
				// Doing bcs:
				bc_Fixed_general(omb, bcVals, imin, imax);
#endif
			}

		} else {

			int imin[3] = {gdata.mx[0]+1, -rim, -rim};
			int imax[3] = {gdata.mx[0]+rim, gdata.mx[1]+rim, gdata.mx[2]+rim};

			if(generic_vars) {
				NumMatrix<REAL,3> & bcVals = bcVals_xe[q];
				// Doing bcs:
				bc_Fixed_general(omb, bcVals, imin, imax);
#if (OMS_USER == TRUE)
			} else {
				NumMatrix<REAL,3> & bcVals = bcVals_User_xe[q];
				// Doing bcs:
				bc_Fixed_general(omb, bcVals, imin, imax);
#endif
			}

		}
	} else if(dir == 1) { // y-direction
		if(above == 0) {

			int imin[3] = {-rim, -rim, -rim};
			int imax[3] = {gdata.mx[0]+rim, -1, gdata.mx[2]+rim};

			// Do not overwrite B on face between computing domain and
			// ghost cells.
			if(q==pr.q_Bx) {
				imax[1] = -2;
			}

			if(generic_vars) {
				NumMatrix<REAL,3> & bcVals = bcVals_yb[q];
				// Doing bcs:
				bc_Fixed_general(omb, bcVals, imin, imax);
#if (OMS_USER == TRUE)
			} else {
				NumMatrix<REAL,3> & bcVals = bcVals_User_yb[q];
				// Doing bcs:
				bc_Fixed_general(omb, bcVals, imin, imax);
#endif
			}

		} else {

			int imin[3] = {-rim, gdata.mx[1]+1, -rim};
			int imax[3] = {gdata.mx[0]+rim, gdata.mx[1]+rim, gdata.mx[2]+rim};

			if(generic_vars) {
				NumMatrix<REAL,3> & bcVals = bcVals_ye[q];
				// Doing bcs:
				bc_Fixed_general(omb, bcVals, imin, imax);
#if (OMS_USER == TRUE)
			} else {
				NumMatrix<REAL,3> & bcVals = bcVals_User_ye[q];
				// Doing bcs:
				bc_Fixed_general(omb, bcVals, imin, imax);
#endif
			}

		}
	} else { // z-direction
		if(above == 0) {

			int imin[3] = {-rim, -rim, -rim};
			int imax[3] = {gdata.mx[0]+rim, gdata.mx[1]+rim, -1};

			// Do not overwrite B on face between computing domain and
			// ghost cells.
			if(q==pr.q_Bz) {
				imax[2] = -2;
			}

			if(generic_vars) {
				NumMatrix<REAL,3> & bcVals = bcVals_zb[q];
				// Doing bcs:
				bc_Fixed_general(omb, bcVals, imin, imax);
#if (OMS_USER == TRUE)
			} else {
				NumMatrix<REAL,3> & bcVals = bcVals_User_zb[q];
				// Doing bcs:
				bc_Fixed_general(omb, bcVals, imin, imax);
#endif
			}

		} else {

			int imin[3] = {-rim, -rim, gdata.mx[2]+1};
			int imax[3] = {gdata.mx[0]+rim, gdata.mx[1]+rim, gdata.mx[2]+rim};

			if(generic_vars) {
				NumMatrix<REAL,3> & bcVals = bcVals_ze[q];
				// Doing bcs:
				bc_Fixed_general(omb, bcVals, imin, imax);
#if (OMS_USER == TRUE)
			} else {
				NumMatrix<REAL,3> & bcVals = bcVals_User_ze[q];
				// Doing bcs:
				bc_Fixed_general(omb, bcVals, imin, imax);
#endif
			}

		}
	}


}

void gridFunc::bc_Fixed_general(NumMatrix<REAL,3> &omb,
                                NumMatrix<REAL,3> &bcVals,
                                int imin[3], int imax[3]) {

	for (int iz = imin[2]; iz <= imax[2]; iz++){
		for (int iy = imin[1]; iy <= imax[1]; ++iy){
			for (int ix = imin[0]; ix <= imax[0]; ++ix){

				omb(ix,iy,iz) = bcVals(ix,iy,iz);

			}
		}
	}
}



// void gridFunc::dataout(Data &gdata, string &filename, ProblemType & Problem,
//                        int numout, bool isfloat)
void gridFunc::dataout(Data &gdata,  Hdf5Stream &h5out, ProblemType & Problem,
                       int numout, bool isfloat, bool is_collective)
{

	// Hdf5Stream h5out(filename, numout, gdata.rank);
  
	Problem.WriteToH5(h5out);

	int rim(0);
	if(isfloat) {
		rim = BOUT_FLT;
	} else {
		rim = BOUT_DBL;
	}

	// Correction for boundary points if such are taken into account
	double xmin[3];
	// xmin[0] = gdata.xb[0]-rim*gdata.dx[0];
	// xmin[1] = gdata.xb[1]-rim*gdata.dx[1];
	// xmin[2] = gdata.xb[2]-rim*gdata.dx[2];
	xmin[0] = gdata.getCen_x(-rim);
	xmin[1] = gdata.getCen_y(-rim);
	xmin[2] = gdata.getCen_z(-rim);
  
	int itime = static_cast<int>(gdata.time/gdata.dt);
    
#ifdef parallel
	h5out.AddGlobalAttr("procs",gdata.nproc,3);
#if(HDF_PARALLEL_IO == CRONOS_OFF)
	h5out.AddGlobalAttr("coords",gdata.coords,3);
	h5out.AddGlobalAttr("rank",gdata.rank);
#else
	if(is_collective) {
		int coords[3] = {0,0,0};
		int rank_global = 0;
		h5out.AddGlobalAttr("coords",coords,3);
		h5out.AddGlobalAttr("rank",rank_global);
	} else {
		h5out.AddGlobalAttr("coords",gdata.coords,3);
		h5out.AddGlobalAttr("rank",gdata.rank);
	}
#endif
#else
	int nproc[3] = {1,1,1};
	int coords[3] = {0,0,0};
	h5out.AddGlobalAttr("procs",nproc,3);
	h5out.AddGlobalAttr("coords",coords,3);
	h5out.AddGlobalAttr("rank",gdata.rank);
#endif

#if (FLUID_TYPE == CRONOS_HYDRO)
	h5out.AddGlobalAttr("fluid_type", "hydro");
#elif (FLUID_TYPE == CRONOS_MHD)
	h5out.AddGlobalAttr("fluid_type", "mhd");
#elif (FLUID_TYPE == CRONOS_RHD)
	h5out.AddGlobalAttr("fluid_type", "rhd");
#else
	h5out.AddGlobalAttr("fluid_type", "N/A");
#endif

	h5out.AddGlobalAttr("geometry",GEOM);
	h5out.AddGlobalAttr("itime",itime);
	h5out.AddGlobalAttr("timestep",gdata.tstep);
	h5out.AddGlobalAttr("time",gdata.time);
	h5out.AddGlobalAttr("rim",rim);
	h5out.AddGlobalAttr("dt",gdata.dt);
	h5out.AddGlobalAttr("xmin",xmin,3);
#if (NON_LINEAR_GRID == _CRONOS_ON)
	h5out.AddGlobalAttr("dx",gdata.dx,3);
#endif
	h5out.AddGlobalAttr("CflNumber",gdata.cfl);

	// Prepare mx and so on for parallel output if necessary
#ifdef parallel
#if(HDF_PARALLEL_IO == CRONOS_ON)
	NumArray<int> mx_global(3), mx_local(3), rank_shift(3);
	//	NumArray<int> rank_pos(3);
	NumArray<float> xmin_global(3), delx(3);
	for(int dir=0; dir<3; ++dir) {
		mx_global(dir) = gdata.global_mx[dir];
		mx_local(dir) = gdata.mx[dir];
		rank_shift(dir) = gdata.get_RankShift(dir);
		// rank_pos(dir) = gdata.coords[dir];
		xmin_global(dir) = gdata.get_pos_global(dir, -rim, 0);
		delx(dir) = gdata.dx[dir];
	}
#endif
#endif

	// In case of multifluid simulation do separate output for individual fluids
#if(FLUID_TYPE == CRONOS_MULTIFLUID)
	int numFluids = gdata.fluids->get_numFluids();
	for(int iFluid=0; iFluid<numFluids; ++iFluid) { // Loop over all fluids
		string groupName = gdata.fluids->fluids[iFluid].get_Name();
		hid_t group = h5out.AddGroup(groupName);

		// Get min and max output index
		int qmin = gdata.fluids->fluids[iFluid].get_IndexGlobal(0);
		int qmax = gdata.fluids->fluids[iFluid].get_N_OMINT() + qmin;

		// Check whether magnetized fluid -- if so vector potential needs to be stored
		bool has_vecPot;
		if(gdata.fluids->get_fluidType(iFluid) == CRONOS_MHD) {
			has_vecPot = true;
		} else {
			has_vecPot = false;
		}
		if(has_vecPot) {
			qmax += 3; // Add 3 fields for vector potential
		}

//		if(iFluid==numFluids-1) {
//			qmax = numout-n_omIntUser;
//		}

		int qmin_user = gdata.fluids->fluids[iFluid].get_IndexGlobalUser(0);
		int qmax_user = gdata.fluids->fluids[iFluid].get_N_OMINT_USER() + qmin;

#else
#if(CRONOS_OUTPUT_COMPATIBILITY == CRONOS_ON)
		hid_t group = h5out.get_defaultGroup();
#else
		string groupName = gdata.fluid.get_Name();
		hid_t group = h5out.AddGroup(groupName);
#endif
		int qmin = 0;
		int qmax = numout-n_omIntUser;

		int qmin_user = 0;
		int qmax_user = n_omIntUser;
#endif

		int q_index = -1;
 	for (int q = qmin; q < qmax; ++q) {
    
		int qout = q;
#if(FLUID_TYPE == CRONOS_MULTIFLUID)
		if(has_vecPot && (q > qmax-4)) {
			// Use index of vector potential:
			qout = q - qmax + 3 + n_omIntAll + N_ADD;
		}
		q_index = qout;
//		cout << " Using " << q << " " << qout << " " << q_index << endl;
#else
		if(q >= n_omInt) {
			qout = q+N_ADD;
		}
#endif

		/*
		if(gdata.mag || (qout <4 || qout>6)) {
		*/

		string dsetName = gdata.om[qout].getName();
#if (USE_COROTATION == CRONOS_ON)
		if(dsetName=="v_x_Corot") dsetName = "v_x";
		if(dsetName=="v_y_Corot") dsetName = "v_y";
		if(dsetName=="v_z_Corot") dsetName = "v_z";
#endif

		if(isfloat) {
			NumMatrix<float,3> data(Index::set(-rim,-rim,-rim),
			                        Index::set(gdata.mx[0]+rim,gdata.mx[1]+rim,
			                                   gdata.mx[2]+rim));
			data = gdata.float_data(qout,rim,true);
			// Possibility of parallel IO into one file so far only
			// for float data (because only those are analysed
			// directly)


#ifdef parallel
#if(HDF_PARALLEL_IO == CRONOS_ON)
			h5out.Write3DMatrix_withMPI_IO(dsetName, data,
			                               mx_global, mx_local,
			                               rank_shift,
			                               //			                               rank_pos,
			                               xmin_global, delx, group);
#else
			h5out.Write3DMatrix(dsetName, data, xmin, gdata.dx, group);
#endif
#else
			h5out.Write3DMatrix(dsetName, data, xmin, gdata.dx, group);
#endif
			// if possible add unit to field
			if(Problem.TrafoNorm != NULL) {

#if(FLUID_TYPE == CRONOS_MULTIFLUID)
				string fieldName;
				if(qout < gdata.fluids->fluids[iFluid].get_N_OMINT()) {
					fieldName = gdata.fluids->fluids[iFluid].get_fieldName(qout);
				} else {
					fieldName = dsetName;
				}
#else
				string fieldName;
				if(qout < gdata.fluid.get_N_OMINT()) {
					fieldName = gdata.fluid.get_fieldName(qout);
					fieldName = gdata.om[qout].getName();
				} else {
					fieldName = dsetName;
				}
#endif
				h5out.AddAttributeToArray(dsetName, "NormFactor",
						Problem.TrafoNorm->get_normConst(fieldName), group);
				h5out.AddAttributeToArray(dsetName, "NormUnit",
						Problem.TrafoNorm->get_unitNormConst(fieldName), group);
				if(gdata.rank==0 && false) {
					cout << " Setting Norm: " << Problem.TrafoNorm->get_normConst(fieldName) << " ";
					cout << Problem.TrafoNorm->get_unitNormConst(fieldName) << " for field: ";
					cout << fieldName << endl;
				}
			}

		} else {
			NumMatrix<double,3> data(Index::set(-rim,-rim,-rim),
			                         Index::set(gdata.mx[0]+rim, gdata.mx[1]+rim,
			                                    gdata.mx[2]+rim));
			
			for (int k = -rim; k <= gdata.mx[2]+rim; k++) {
				for (int j = -rim; j <= gdata.mx[1]+rim; j++) {
					for (int i = -rim; i <= gdata.mx[0]+rim; i++) {
						data(i,j,k) = gdata.om[qout](i,j,k);
					}
				}
			}
#ifdef parallel
#if(HDF_PARALLEL_IO == CRONOS_ON)
			h5out.Write3DMatrix_withMPI_IO(dsetName, data,
			                               mx_global, mx_local,
			                               rank_shift,
			                               //			                               rank_pos,
			                               xmin_global, delx, group, q_index);
#else
			h5out.Write3DMatrix(dsetName, data, xmin, gdata.dx, group, q_index);
#endif
#else
			h5out.Write3DMatrix(dsetName, data, xmin, gdata.dx, group, q_index);
#endif
		}

		/*
		  }
		*/

	}

 	if(isfloat) {
 		// If desired write output of flag for carbuncle problem
 		if(gdata.use_carbuncleFlag) {
 			if(value_exists("write_carbuncleFlag")) {
 				NumMatrix<float,3> data(Index::set(-rim,-rim,-rim),
 						Index::set(gdata.mx[0]+rim, gdata.mx[1]+rim, gdata.mx[2]+rim));
 				for (int k = -rim; k <= gdata.mx[2]+rim; k++) {
 					for (int j = -rim; j <= gdata.mx[1]+rim; j++) {
 						for (int i = -rim; i <= gdata.mx[0]+rim; i++) {
 							data(i,j,k) = static_cast<float>(gdata.carbuncleFlag(i,j,k));
 						}
 					}
 				}
#ifdef parallel
#if(HDF_PARALLEL_IO == CRONOS_ON)
 				h5out.Write3DMatrix_withMPI_IO("carbuncle_flag", data,
 						mx_global, mx_local,
 						rank_shift,
 						xmin_global, delx, group);
#else
	h5out.Write3DMatrix("carbuncle_flag", data, xmin, gdata.dx, group);
#endif
#else
	h5out.Write3DMatrix("carbuncle_flag", data, xmin, gdata.dx, group);
#endif
 			}
 		}
 	}




#if (OMS_USER == TRUE)
	for (int q = qmin_user; q < qmax_user; ++q) {
    
		if(isfloat) {
			NumMatrix<float,3> data(Index::set(-rim,-rim,-rim),
			                        Index::set(gdata.mx[0]+rim,gdata.mx[1]+rim,
			                                   gdata.mx[2]+rim));
			data = gdata.float_data(q,rim,false);
#ifdef parallel
#if(HDF_PARALLEL_IO == CRONOS_ON)
			h5out.Write3DMatrix_withMPI_IO(gdata.om_user[q].getName(), data,
			                               mx_global, mx_local,
			                               rank_shift,
			                               xmin_global, delx, group);
#else
			h5out.Write3DMatrix(gdata.om_user[q].getName(), data,
			                    xmin, gdata.dx, group);
#endif
#else
			h5out.Write3DMatrix(gdata.om_user[q].getName(), data,
			                    xmin, gdata.dx, group);
#endif
		} else {
			NumMatrix<double,3> data(Index::set(-rim,-rim,-rim),
			                         Index::set(gdata.mx[0]+rim,
			                                    gdata.mx[1]+rim,
			                                    gdata.mx[2]+rim));
			
			for (int k = -rim; k <= gdata.mx[2]+rim; k++) {
				for (int j = -rim; j <= gdata.mx[1]+rim; j++) {
					for (int i = -rim; i <= gdata.mx[0]+rim; i++) {
						data(i,j,k) = gdata.om_user[q](i,j,k);
					}
				}
			}
#ifdef parallel
#if(HDF_PARALLEL_IO == CRONOS_ON)
			h5out.Write3DMatrix_withMPI_IO(gdata.om_user[q].getName(), data,
			                               mx_global, mx_local,
			                               rank_shift,
			                               xmin_global, delx, group);
#else
			h5out.Write3DMatrix(gdata.om_user[q].getName(), data,
					xmin, gdata.dx, group);
#endif
#else
			h5out.Write3DMatrix(gdata.om_user[q].getName(), data,
					xmin, gdata.dx, group);
#endif

		}

	}
#endif
#if(FLUID_TYPE == CRONOS_MULTIFLUID)
	// Now close the new group:
	h5out.CloseGroup(group);
	}
#else
#if(CRONOS_OUTPUT_COMPATIBILITY != CRONOS_ON)
	h5out.CloseGroup(group);
#endif
#endif

	// If we are using a non-linear grid -- write specific grid:
#if(NON_LINEAR_GRID == CRONOS_ON)
	gdata.add_nonLinGridToHdf(h5out, rim);
#endif

}





// void gridFunc::datain(Data &gdata, string &filename, ProblemType &Problem)
// {
void gridFunc::datain_old(Data &gdata, Hdf5iStream & h5in, string &filename,
                      ProblemType &Problem)
{

	// Hdf5iStream h5in(filename);


	Problem.ReadFromH5(h5in);

#ifdef parallel
	int nprocOld[DIM];
	int facProc[DIM];
	int coordsOld[DIM];
	h5in.ReadGlobalAttr("procs",*nprocOld);
	h5in.ReadGlobalAttr("coords",*coordsOld);

	for(int i=0; i<DIM; ++i) {
		facProc[i] = 1;
		if(nprocOld[i] != gdata.nproc[i]) {
			if(gdata.nproc[i] > nprocOld[i]) {
				facProc[i] = gdata.nproc[i]/nprocOld[i];
			} else {
				facProc[i] = -nprocOld[i]/gdata.nproc[i];
			}
		}
	}
#endif
	h5in.ReadGlobalAttr("time",gdata.time);
	h5in.ReadGlobalAttr("dt",gdata.dt);
	int rim(0); // Boundary points to be taken into account
	h5in.ReadGlobalAttr("rim",rim);
	h5in.ReadGlobalAttr("CflNumber",gdata.cfl);
	//   h5in.ReadGlobalAttr("xmin",*xb);
  
	int mxsm[DIM];

	// get group name:
#if(FLUID_TYPE == CRONOS_MULTIFLUID)
	string groupName = gdata.fluids->fluids[0].get_Name();
#else
	string groupName = gdata.fluid.get_Name();
#endif

	h5in.getSize(groupName, gdata.om[0].getName(), mxsm, DIM);

	for(int i=0; i<DIM; ++i) {
		mxsm[i] -= 2*rim+1;
#ifdef parallel
		if(facProc[i] != 1) {
			if(facProc[i] > 1) {
				mxsm[i] /= facProc[i];
			} else {
				mxsm[i] *= -facProc[i];
			}
		}
#endif
	}

	NumMatrix<double,3> Data;
	// Check if grid ratio can be devided by zero:
	for(int dir=0; dir<3; ++dir) {
		if(mxsm[dir] > 0) {
			if(gdata.mx[dir]%mxsm[dir] != 0) {
				throw CException(" New grid-extent no devisor of old one ");
			}
		}
	}

	int fac(1);
	int facloc[3];
	for(int dir=0; dir<3; ++dir) {
		if(mxsm[dir] > 0) {
			facloc[dir] = gdata.mx[dir]/mxsm[dir];
		} else {
			facloc[dir] = 1;
		}
	}
	if(facloc[0] != facloc[1] || facloc[1] != facloc[2]) {
		if(gdata.rank == 0) {
			cerr << " Change of resolution not uniform " << endl;
			cerr << " x-direction:    " << facloc[0] << endl;
			cerr << " y-direction:    " << facloc[1] << endl;
			cerr << " z-direction:    " << facloc[2] << endl;
		}
		exit(-102);
	} else {
		fac = facloc[0];
	}

	if(fac < 1) {
		throw CException(" Grid can only be extended ");
	}

	// Adjust timestep:
	gdata.dt /= fac;


	int numin(0);
	h5in.ReadGlobalAttr("Entries",numin);
	int diffnum = 0;
	if(!gdata.mag) {
		diffnum = 3;
	}
	if(numin < N_OMINT_ALL + N_SUBS - diffnum) {
		cerr << " Not enough entries in input file " << endl;
		cerr << numin << " " << n_omInt << " " << N_SUBS << endl;
		exit(-54);
	}


#ifdef parallel
	for(int i=0; i<DIM; ++i) {
		if(facProc[i] != 1) {
			if(facProc[i] > 1) {
				mxsm[i] *=  facProc[i];
			} else {
				mxsm[i] /= -facProc[i];
			}
		}	
	}
#endif

	int shift[DIM];
	for(int i=0; i<DIM; ++i) {
		shift[i] = 0;
	}
	for (int q = 0; q < n_omInt+n_omIntUser+N_SUBS; ++q) {

//		bool is_generic = true;

		int qin = q;
		if(q < n_omInt) {
			qin = q;
		} else if(q >= n_omInt && q < n_omInt+N_SUBS) {
			qin = q+N_ADD;
		} else {
			qin = q-(n_omInt+N_SUBS);
//			is_generic = false;
		}
    
		/*
		  if(gdata.mag || (qin <4 || qin>6)) {
		*/

		Data.resize(Index::set(-rim, -rim, -rim),
		            Index::set(mxsm[0]+rim, mxsm[1]+rim, mxsm[2]+rim));
		Data.clear();
		
#if (OMS_USER == TRUE)
		if(q < n_omInt+N_SUBS) {
#endif
			gdata.om[qin].rename("dump");
#if (OMS_USER == TRUE)
		} else {
			gdata.om_user[qin].rename("dump");
		}
#endif

		int qshift = 0;
		// if(!gdata.mag && qin > 6) {
		// This should only have been applicable to the non-user fields?!
		// I do not see how this was ever sensible???
		// if(!gdata.mag && qin > 6) {
		// 	qshift = 3;
		// }

#if (OMS_USER == TRUE)
		if(q < n_omInt+N_SUBS) {
#endif
			gdata.om[qin].rename(h5in.GetDatasetName(q-qshift));
			h5in.Read3DMatrix(gdata.om[qin].getName(), Data);
#if (OMS_USER == TRUE)
		} else {
			gdata.om_user[qin].rename(h5in.GetDatasetName(q-qshift));
			h5in.Read3DMatrix(gdata.om_user[qin].getName(), Data);
		}
#endif

#ifdef parallel
		Assign(gdata, Data, mxsm, facProc, coordsOld, shift, rim);
#endif

		/*
		  }
		*/

		REAL facred = fac;
#ifdef parallel
		int mxloc[3] = {mxsm[0]/facProc[0], mxsm[1]/facProc[1],
		                mxsm[2]/facProc[2]};
#else
		int mxloc[3] = {mxsm[0], mxsm[1], mxsm[2]};
#endif
		if(fac != 1 && gdata.rank == 0) {
			cout << " Prolongating: ";
		}
		while(facred != 1) {
			if(gdata.rank == 0) {
				cout << facred << " ";
			}
			Prolongate(gdata, Data, qin, mxloc, rim);
			facred /= 2;
		}
		if(fac != 1) {
			cout << endl;
		}

		int imin[DIM];
		int imax[DIM];
		for(int i=0; i<DIM; ++i) {


#ifdef parallel
			if(facProc[i] > 0) {
				imin[i] = -rim;
				imax[i] = gdata.mx[i]+rim;
			} else {
				imin[i] = 0;
				imax[i] = -gdata.mx[i]/facProc[i];
			}
#else
			imin[i] = -rim;
			imax[i] = gdata.mx[i]+rim;
#endif
		}
  
#ifdef parallel
		if(q < n_omInt+N_SUBS) {

			if(facProc[0] < 0) {
				imin[0] = 0;
				imax[0] = -gdata.mx[0]/facProc[0];
				if (gdata.om[qin].getName() == "B_x") {
					imin[0] = -1;
				} else if (gdata.om[qin].getName() == "A_y") {
					imax[0] += 1;
				} else if (gdata.om[qin].getName() == "A_z") {
					imax[0] += 1;
				}
			}
			
			if(facProc[1] < 0) {
				imin[1] = 0;
				imax[1] = -gdata.mx[1]/facProc[1];
				if (gdata.om[qin].getName() == "B_y") {
					imin[1] = 1;
				} else if (gdata.om[qin].getName() == "A_x") {
					imax[1] += 1;
				} else if (gdata.om[qin].getName() == "A_z") {
					imax[1] += 1;
				}
			}
			
			if(facProc[2] < 0) {
				imin[2] = 0;
				imax[2] = -gdata.mx[2]/facProc[2];
				if (gdata.om[qin].getName() == "B_z") {
					imin[2] = 1;
				} else if (gdata.om[qin].getName() == "A_x") {
					imax[2] += 1;
				} else if (gdata.om[qin].getName() == "A_y") {
					imax[2] += 1;
				}
			}
		}
			
#endif

		
		/*
		  if(gdata.mag || (qin <4 || qin>6)) {
		*/

		for(int k=imin[2]; k<=imax[2]; ++k) {
			for(int j=imin[1]; j<=imax[1]; ++j) {
				for(int i=imin[0]; i<=imax[0]; ++i) {
					
#if (OMS_USER == TRUE)
					if(q < n_omInt+N_SUBS) {
#endif
						gdata.om[qin](i+shift[0],j+shift[1],
						              k+shift[2]) = Data(i+shift[0],j+shift[1],k+shift[2]);
#if (OMS_USER == TRUE)
					} else {
						gdata.om_user[qin](i+shift[0],j+shift[1],
						                   k+shift[2]) = Data(i+shift[0],j+shift[1],k+shift[2]);
					}
#endif
					
				}
			}
		}
		/*
		  }
		*/
	}

}



void gridFunc::datain(Data &gdata, Hdf5iStream & h5in, string &filename,
                      ProblemType &Problem)
{

	// Hdf5iStream h5in(filename);


	Problem.ReadFromH5(h5in);

#ifdef parallel
	int nprocOld[DIM];
	int facProc[DIM];
	int coordsOld[DIM];
	h5in.ReadGlobalAttr("procs",*nprocOld);
	h5in.ReadGlobalAttr("coords",*coordsOld);

	for(int i=0; i<DIM; ++i) {
		facProc[i] = 1;
		if(nprocOld[i] != gdata.nproc[i]) {
			if(gdata.nproc[i] > nprocOld[i]) {
				facProc[i] = gdata.nproc[i]/nprocOld[i];
			} else {
				facProc[i] = -nprocOld[i]/gdata.nproc[i];
			}
		}
	}
#endif
	h5in.ReadGlobalAttr("time",gdata.time);
	h5in.ReadGlobalAttr("dt",gdata.dt);
	int rim(0); // Boundary points to be taken into account
	h5in.ReadGlobalAttr("rim",rim);
	h5in.ReadGlobalAttr("CflNumber",gdata.cfl);
	//   h5in.ReadGlobalAttr("xmin",*xb);

	int mxsm[DIM];

	// get group name:
#if(FLUID_TYPE == CRONOS_MULTIFLUID)
	string groupName = gdata.fluids->fluids[0].get_Name();
#else

#if(CRONOS_OUTPUT_COMPATIBILITY == CRONOS_ON)
		string groupName = "Data";
#else
		string groupName = gdata.fluid.get_Name();
#endif

#endif

	h5in.getSize(groupName, gdata.om[0].getName(), mxsm, DIM);

	for(int i=0; i<DIM; ++i) {
		mxsm[i] -= 2*rim+1;
#ifdef parallel
		if(facProc[i] != 1) {
			if(gdata.get_EdgeGridding()) {
				if(facProc[i] > 1) {
					mxsm[i] = (mxsm[i] + 1)/facProc[i] - 1;
				} else {
					mxsm[i] = (mxsm[i] + 1)*(-facProc[i]) - 1;
				}
			} else {
				if(facProc[i] > 1) {
					mxsm[i] /= facProc[i];
				} else {
					mxsm[i] *= -facProc[i];
				}
			}
		}
#endif
	}

	NumMatrix<double,3> Data;
	// Check if grid ratio can be devided by zero:
	for(int dir=0; dir<3; ++dir) {
		if(mxsm[dir] > 0) {
			if(gdata.get_EdgeGridding()) {
				if((gdata.mx[dir]-1)%(mxsm[dir]-1) != 0) {
					throw CException(" New grid-extent no devisor of old one ");
				}
			} else {
				if(gdata.mx[dir]%mxsm[dir] != 0) {
					throw CException(" New grid-extent no devisor of old one ");
				}
			}
		}
	}

	int fac(1);
	int facloc[3];
	for(int dir=0; dir<3; ++dir) {
		if(mxsm[dir] > 0) {
			if(gdata.get_EdgeGridding()) {
				facloc[dir] = (gdata.mx[dir]+1)/(mxsm[dir]+1);
			} else {
				facloc[dir] = gdata.mx[dir]/mxsm[dir];
			}
		} else {
			facloc[dir] = 1;
		}
	}
	if(facloc[0] != facloc[1] || facloc[1] != facloc[2]) {
		if(gdata.rank == 0) {
			cerr << " Change of resolution not uniform " << endl;
			cerr << " x-direction:    " << facloc[0] << endl;
			cerr << " y-direction:    " << facloc[1] << endl;
			cerr << " z-direction:    " << facloc[2] << endl;
		}
		exit(-102);
	} else {
		fac = facloc[0];
	}

	if(fac < 1) {
		throw CException(" Grid can only be extended ");
	}

	// Adjust timestep:
	gdata.dt /= fac;


	int numin(0);
	h5in.ReadGlobalAttr("Entries",numin);
	int diffnum = 0;
	if(!gdata.mag) {
		diffnum = 3;
	}
	if(numin < N_OMINT_ALL + N_SUBS - diffnum) {
		cerr << " Not enough entries in input file " << endl;
		cerr << numin << " " << n_omInt << " " << N_SUBS << endl;
		exit(-54);
	}


#ifdef parallel
	for(int i=0; i<DIM; ++i) {
		if(facProc[i] != 1) {
			if(gdata.get_EdgeGridding()) {
				if(facProc[i] > 1) {
					mxsm[i] = (mxsm[i] + 1)*facProc[i] - 1;
				} else {
					mxsm[i] = (mxsm[i] + 1)/(-facProc[i]) - 1;
				}
			} else {
				if(facProc[i] > 1) {
					mxsm[i] *=  facProc[i];
				} else {
					mxsm[i] /= -facProc[i];
				}
			}
		}
	}
#endif

	int shift[DIM];
	for(int i=0; i<DIM; ++i) {
		shift[i] = 0;
	}

	// In case of multifluid simulation separate output for individual fluids is done
#if(FLUID_TYPE == CRONOS_MULTIFLUID)
	int numFluids = gdata.fluids->get_numFluids();
	for(int iFluid=0; iFluid<numFluids; ++iFluid) { // Loop over all fluids
		groupName = gdata.fluids->fluids[iFluid].get_Name();

		// Get min and max output index
		int qmin = gdata.fluids->fluids[iFluid].get_IndexGlobal(0);
		int qmax = gdata.fluids->fluids[iFluid].get_N_OMINT() + qmin;
//		if(iFluid==numFluids-1) {
//			qmax = numout-n_omIntUser;
//		}

		// Check whether magnetized fluid -- if so vector potential needs to be stored
		bool has_vecPot;
		if(gdata.fluids->get_fluidType(iFluid) == CRONOS_MHD) {
			has_vecPot = true;
		} else {
			has_vecPot = false;
		}
		if(has_vecPot) {
			qmax += 3; // Add 3 fields for vector potential
		}


		int qmin_user = gdata.fluids->fluids[iFluid].get_IndexGlobalUser(0);
//		int qmax_user = gdata.fluids->fluids[iFluid].get_N_OMINT_USER() + qmin_user;
		int qmax_user = gdata.fluids->fluids[iFluid].get_N_OMINT_USER();

#else
#if(CRONOS_OUTPUT_COMPATIBILITY == CRONOS_ON)
	groupName = "Data";
#else
	groupName = gdata.fluid.get_Name();
#endif
	int qmin = 0;
	int qmax = n_omIntAll + N_SUBS - n_omIntUser;

	int qmin_user = 0;
	int qmax_user = n_omIntUser;
#endif

//	for (int q = 0; q < qmax + qmax_user; ++q) {
	for (int q = qmin; q < qmax + qmax_user; ++q) {

//		bool is_generic = true;

		int qin = q;

#if(FLUID_TYPE == CRONOS_MULTIFLUID)
		qin = q;
		if(has_vecPot && (q > qmax-4)) {
			qin = q - qmax + 3 + n_omIntAll + N_ADD;
		}
#else

		if(q < n_omInt) {
			qin = q;
		} else if(q >= n_omInt && q < n_omInt+N_SUBS) {
			qin = q+N_ADD;
		} else {
			qin = q-(n_omInt+N_SUBS);
//			is_generic = false;
		}

#endif
//		cerr << " Lese: " << q << " " << qin << " " << gdata.om[qin].getName() << endl;
//		cerr << " von " << qmax << " " << qmax_user << endl;

		/*
		  if(gdata.mag || (qin <4 || qin>6)) {
		*/

		Data.resize(Index::set(-rim, -rim, -rim),
		            Index::set(mxsm[0]+rim, mxsm[1]+rim, mxsm[2]+rim));
		Data.clear();

#if (OMS_USER == TRUE)
		if(q < n_omInt+N_SUBS) {
#endif
			gdata.om[qin].rename("dump");
#if (OMS_USER == TRUE)
		} else {
			gdata.om_user[qin].rename("dump");
		}
#endif

		int qshift = 0;
		// if(!gdata.mag && qin > 6) {
		// This should only have been applicable to the non-user fields?!
		// I do not see how this was ever sensible???
		// if(!gdata.mag && qin > 6) {
		// 	qshift = 3;
		// }

#if (OMS_USER == TRUE)
		if(q < n_omInt+N_SUBS) {
#endif

#if(FLUID_TYPE == CRONOS_MULTIFLUID)
			gdata.om[qin].rename(h5in.GetDatasetName(groupName, q-qmin));
#else
			gdata.om[qin].rename(h5in.GetDatasetName(q-qshift));
#endif
//			cout << " Lese2 " << q << " " << qin << " " << gdata.om[qin].getName() << " " << OMS_USER << endl;
			h5in.Read3DMatrix(groupName, gdata.om[qin].getName(), Data);
#if (OMS_USER == TRUE)
		} else {
			gdata.om_user[qin].rename(h5in.GetDatasetName(q-qshift));
			h5in.Read3DMatrix(groupName, gdata.om_user[qin].getName(), Data);
		}
#endif

#ifdef parallel
		Assign(gdata, Data, mxsm, facProc, coordsOld, shift, rim);
#endif

		/*
		  }
		*/

		REAL facred = fac;
#ifdef parallel
		int mxloc[3] = {mxsm[0]/facProc[0], mxsm[1]/facProc[1],
		                mxsm[2]/facProc[2]};
#else
		int mxloc[3] = {mxsm[0], mxsm[1], mxsm[2]};
#endif
		if(fac != 1 && gdata.rank == 0) {
			cout << " Prolongating: ";
		}
		while(facred != 1) {
			if(gdata.rank == 0) {
				cout << facred << " ";
			}
			Prolongate(gdata, Data, qin, mxloc, rim);
			facred /= 2;
		}
		if(fac != 1) {
			cout << endl;
		}

		int imin[DIM];
		int imax[DIM];
		for(int i=0; i<DIM; ++i) {


#ifdef parallel
			if(facProc[i] > 0) {
				imin[i] = -rim;
				imax[i] = gdata.mx[i]+rim;
			} else {
				imin[i] = 0;
				imax[i] = -gdata.mx[i]/facProc[i];
			}
#else
			imin[i] = -rim;
			imax[i] = gdata.mx[i]+rim;
#endif
		}

#ifdef parallel
		if(q < n_omInt+N_SUBS) {

			if(facProc[0] < 0) {
				imin[0] = 0;
				imax[0] = -gdata.mx[0]/facProc[0];
				if (gdata.om[qin].getName() == "B_x") {
					imin[0] = -1;
				} else if (gdata.om[qin].getName() == "A_y") {
					imax[0] += 1;
				} else if (gdata.om[qin].getName() == "A_z") {
					imax[0] += 1;
				}
			}

			if(facProc[1] < 0) {
				imin[1] = 0;
				imax[1] = -gdata.mx[1]/facProc[1];
				if (gdata.om[qin].getName() == "B_y") {
					imin[1] = 1;
				} else if (gdata.om[qin].getName() == "A_x") {
					imax[1] += 1;
				} else if (gdata.om[qin].getName() == "A_z") {
					imax[1] += 1;
				}
			}

			if(facProc[2] < 0) {
				imin[2] = 0;
				imax[2] = -gdata.mx[2]/facProc[2];
				if (gdata.om[qin].getName() == "B_z") {
					imin[2] = 1;
				} else if (gdata.om[qin].getName() == "A_x") {
					imax[2] += 1;
				} else if (gdata.om[qin].getName() == "A_y") {
					imax[2] += 1;
				}
			}
		}

#endif


		/*
		  if(gdata.mag || (qin <4 || qin>6)) {
		*/
		for(int k=imin[2]; k<=imax[2]; ++k) {
			for(int j=imin[1]; j<=imax[1]; ++j) {
				for(int i=imin[0]; i<=imax[0]; ++i) {

#if (OMS_USER == TRUE)
					if(q < n_omInt+N_SUBS) {
#endif
						gdata.om[qin](i+shift[0],j+shift[1],
						              k+shift[2]) = Data(i+shift[0],j+shift[1],k+shift[2]);
#if (OMS_USER == TRUE)
					} else {
						gdata.om_user[qin](i+shift[0],j+shift[1],
						                   k+shift[2]) = Data(i+shift[0],j+shift[1],k+shift[2]);
					}
#endif

				}
			}
		}
		/*
		  }
		*/
	}
#if(FLUID_TYPE == CRONOS_MULTIFLUID)
	}
#endif

}


int gridFunc::get_powerTwo(int val) {
	int powerTwo=0;
	while(val >= 2) {
		val /= 2;
		powerTwo++;
	}
	return powerTwo;
}

double gridFunc::datain_collective(Data &gdata, Hdf5iStream & h5in, string &filename,
                      ProblemType &Problem)
{

	// Hdf5iStream h5in(filename);


	Problem.ReadFromH5(h5in);

	h5in.ReadGlobalAttr("time",gdata.time);
	h5in.ReadGlobalAttr("dt",gdata.dt);
	int rim_data(0); // Boundary points to be taken into account
	h5in.ReadGlobalAttr("rim",rim_data);
	h5in.ReadGlobalAttr("CflNumber",gdata.cfl);
	//   h5in.ReadGlobalAttr("xmin",*xb);

	// check if the vector potential is supposed to be used
	int use_vecpot;
	if(h5in.doesAttrExist("HyperbolicSolver_IntegrateA")) {
		h5in.ReadGlobalAttr("HyperbolicSolver_IntegrateA", use_vecpot);
	} else {
		use_vecpot = 0;
	}

	// get group name:
#if(FLUID_TYPE == CRONOS_MULTIFLUID)
	string groupName = gdata.fluids->fluids[0].get_Name();
#else

#if(CRONOS_OUTPUT_COMPATIBILITY == CRONOS_ON)
	string groupName = "Data";
#else
	string groupName = gdata.fluid.get_Name();
#endif

#endif

	// Get grid size in hdf5 file
	int mx_data[DIM];
	int Nx_data[DIM];
	h5in.getSize(groupName, gdata.om[0].getName(), Nx_data, DIM);

	// Get number of ghost cells used for writing:

	// Take boundary cells into account:
	for(int dir=0; dir<DIM; ++dir) {
		Nx_data[dir] -= 2*rim_data;
		mx_data[dir] = Nx_data[dir]-1;
//		mx_data[dir] -= 2*rim_data+1;
//		Nx_data[dir] = mx_data[dir]-1;
	}

	// Now, determine the respective grid ratios
	NumArray<int> ratio(DIM);
	REAL eps = 1.e-42;
	for(int dir=0; dir<DIM; ++dir) {
		if(gdata.get_EdgeGridding()) {
			ratio(dir) = gdata.get_numCells_global(dir)/Nx_data[dir];
//			ratio(dir) = Nx_data[dir]/(gdata.get_numCells_global(dir));
			// Test if we have an integer number
			REAL d_ratio = gdata.get_numCells_global(dir)/(1.*Nx_data[dir]-eps);
			if(std::abs(ratio(dir)-d_ratio) > 1.e-12) {
				cerr << " Error in datain_collective - ratio of grids is not an integer " << endl;
				exit(3);
			}
//			cout << " ratio: " << ratio(dir) << " " << Nx_data[dir] << " " << gdata.get_numCells_global(dir) << endl;
		} else {
			ratio(dir) = gdata.global_mx[dir]/(mx_data[dir]);
//			ratio(dir) = mx_data[dir]/(gdata.global_mx[dir]);
			// Test if we have an integer number
			REAL d_ratio = gdata.global_mx[dir]/(1.*mx_data[dir]-eps);
//			REAL d_ratio = mx_data[dir]/(1.*gdata.global_mx[dir]-eps);
			if(ratio(dir)-d_ratio > 1.e-12) {
				cerr << " Error in datain_collective - ratio of grids is not an integer " << endl;
				exit(3);
			}
		}
	}

	// Test whether ratio is acceptable (can be divided by 2)
	NumArray<int> levels(DIM); // number of levels for resolution change
	for(int dir=0; dir<DIM; ++dir) {
		levels(dir) =  get_powerTwo(ratio(dir));
		// value as expected from level
		int value = cronos::power(2,levels(dir));
		// comparison to real value
		cout << " level " << dir << " " << levels(dir) << endl;

		cout << " val " << value << " " << ratio(dir) << " " << dir << endl;

		if(value != ratio(dir)) {
			cerr << " Error in datain_collective - ratio of grids is not power of 2 " << endl;
			exit(3);
		}
		cout << " levels " << levels(dir) << endl;
	}

	NumArray<int> mx_data_local(DIM);
#ifdef parallel
#if(HDF_PARALLEL_IO == CRONOS_ON)
	// in MPI-parallel case: set local grid for reading part of hdf5 file
	NumArray<int> mx_global;
	NumArray<int> rank_shift;
	mx_data_local.resize(DIM);
	rank_shift.resize(DIM);

	for(int dir=0; dir<DIM; ++dir) {
		mx_data_local(dir) = gdata.mx[dir];
	    rank_shift(dir) = gdata.get_RankShift(dir);
	}

	// Adapt this for smaller grid size in file as compared to what is
	// given in cat file
	if(gdata.get_EdgeGridding()) {
		for(int dir=0; dir<DIM; ++dir) {
			mx_data_local(dir) = (mx_data_local(dir)+1)/ratio(dir)-1;
			rank_shift(dir) /= ratio(dir);
		}
	} else {
		for(int dir=0; dir<DIM; ++dir) {
			mx_data_local(dir) /= ratio(dir);
			rank_shift(dir) /= ratio(dir);
		}
	}

#else
	for(int dir=0; dir<DIM; ++dir) {
		mx_data_local(dir) = mx_data[dir];
	}
#endif
#else
	for(int dir=0; dir<DIM; ++dir) {
		mx_data_local(dir) = mx_data[dir];
	}
#endif


	NumMatrix<double,3> Data;


//	int numin(0);
//	h5in.ReadGlobalAttr("Entries",numin);
//	int diffnum = 0;
//	if(!gdata.mag) {
//		diffnum = 3;
//	}
//	if(numin < N_OMINT_ALL + N_SUBS - diffnum) {
//		cerr << " Not enough entries in input file " << endl;
//		cerr << numin << " " << n_omInt << " " << N_SUBS << endl;
//		exit(-54);
//	}


	// This might be completely unnnecessary, now
	int shift[DIM];
	for(int i=0; i<DIM; ++i) {
		shift[i] = 0;
	}

	// In case of multifluid simulation separate output for individual fluids is done
#if(FLUID_TYPE == CRONOS_MULTIFLUID)
	int numFluids = gdata.fluids->get_numFluids();
	for(int iFluid=0; iFluid<numFluids; ++iFluid) { // Loop over all fluids
		groupName = gdata.fluids->fluids[iFluid].get_Name();

		// Get min and max output index
		int qmin = gdata.fluids->fluids[iFluid].get_IndexGlobal(0);
		int qmax = gdata.fluids->fluids[iFluid].get_N_OMINT() + qmin;
//		if(iFluid==numFluids-1) {
//			qmax = numout-n_omIntUser;
//		}

		// Check whether magnetized fluid -- if so vector potential needs to be stored
		bool has_vecPot;
		if(gdata.fluids->get_fluidType(iFluid) == CRONOS_MHD) {
			has_vecPot = true;
		} else {
			has_vecPot = false;
		}
		if(has_vecPot) {
			qmax += 3; // Add 3 fields for vector potential
		}


		int qmin_user = gdata.fluids->fluids[iFluid].get_IndexGlobalUser(0);
//		int qmax_user = gdata.fluids->fluids[iFluid].get_N_OMINT_USER() + qmin_user;
		int qmax_user = gdata.fluids->fluids[iFluid].get_N_OMINT_USER();

#else
#if(CRONOS_OUTPUT_COMPATIBILITY == CRONOS_ON)
	groupName = "Data";
#else
	groupName = gdata.fluid.get_Name();
#endif
	int qmin = 0;
	int qmax = n_omIntAll + N_SUBS - n_omIntUser;

	int qmin_user = 0;
	int qmax_user = n_omIntUser;
#endif




	// Show all available fields in file

	for (int q = qmin; q < qmax + qmax_user; ++q) {

//		bool is_generic = true;

		int qin = q;

#if(FLUID_TYPE == CRONOS_MULTIFLUID)
		qin = q;
		if(has_vecPot && (q > qmax-4)) {
			qin = q - qmax + 3 + n_omIntAll + N_ADD;
		}
#else

		if(q < n_omInt) {
			qin = q;
		} else if(q >= n_omInt && q < n_omInt+N_SUBS) {
			qin = q+N_ADD;
		} else {
			qin = q-(n_omInt+N_SUBS);
//			is_generic = false;
		}

#endif
		std::string dataName;


		int qshift = 0;
		// if(!gdata.mag && qin > 6) {
		// This should only have been applicable to the non-user fields?!
		// I do not see how this was ever sensible???
		// if(!gdata.mag && qin > 6) {
		// 	qshift = 3;
		// }

#if (OMS_USER == TRUE)
		if(q < n_omInt+N_SUBS) {
#endif


			string fieldName;

#if(FLUID_TYPE == CRONOS_MULTIFLUID)
			fieldName = h5in.GetDatasetName(groupName, q-qmin);
#elif(CRONOS_OUTPUT_COMPATIBILITY == CRONOS_ON)
			fieldName = h5in.GetDatasetName(q-qshift);
#else
			fieldName = h5in.GetDatasetName(groupName, q-qshift);
#endif
			dataName = fieldName;
			if(gdata.rank==0) {
				cout << " Available: " << dataName << " as " << q << endl;
			}
#if (OMS_USER == TRUE)
		}
#endif
	}






//	for (int q = 0; q < qmax + qmax_user; ++q) {
	for (int q = qmin; q < qmax + qmax_user; ++q) {
//		bool is_generic = true;

		int qin = q;

#if(FLUID_TYPE == CRONOS_MULTIFLUID)
		qin = q;
		if(has_vecPot && (q > qmax-4)) {
			qin = q - qmax + 3 + n_omIntAll + N_ADD;
		}
#else

		if(q < n_omInt) {
			qin = q;
		} else if(q >= n_omInt && q < n_omInt+N_SUBS) {
			qin = q+N_ADD;
		} else {
			qin = q-(n_omInt+N_SUBS);
//			is_generic = false;
		}

#endif
//		cerr << " Lese: " << q << " " << qin << " " << gdata.om[qin].getName() << endl;
//		cerr << " von " << qmax << " " << qmax_user << endl;

		/*
		  if(gdata.mag || (qin <4 || qin>6)) {
		*/
//		cout << " resize " << endl;
		Data.resize(Index::set(-rim_data, -rim_data, -rim_data),
				Index::set(mx_data_local[0]+rim_data, mx_data_local[1]+rim_data, mx_data_local[2]+rim_data));
//		Data.resize(Index::set(-rim, -rim, -rim),
//		            Index::set(mx_data[0]+rim, mx_data[1]+rim, mx_data[2]+rim));
		Data.clear();
		std::string dataName;

#if (OMS_USER == TRUE)
		if(q < n_omInt+N_SUBS) {
#endif
			gdata.om[qin].rename("dump");
#if (OMS_USER == TRUE)
		} else {
			gdata.om_user[qin].rename("dump");
		}
#endif

		int qshift = 0;
		// if(!gdata.mag && qin > 6) {
		// This should only have been applicable to the non-user fields?!
		// I do not see how this was ever sensible???
		// if(!gdata.mag && qin > 6) {
		// 	qshift = 3;
		// }

#if (OMS_USER == TRUE)
		if(q < n_omInt+N_SUBS) {
#endif

#if(FLUID_TYPE == CRONOS_MULTIFLUID)
			gdata.om[qin].rename(h5in.GetDatasetName(groupName, q-qmin));
#elif(CRONOS_OUTPUT_COMPATIBILITY == CRONOS_ON)
			gdata.om[qin].rename(h5in.GetDatasetName(q-qshift));
#else
//			cerr << " Going for name of " << q-qshift << " " << qin << endl;
			gdata.om[qin].rename(h5in.GetDatasetName(groupName, q-qshift));
#endif
			dataName = gdata.om[qin].getName();

#ifdef parallel
#if(HDF_PARALLEL_IO == CRONOS_ON)
			h5in.Read3DMatrix_parallel(groupName, gdata.om[qin].getName(),
					mx_data_local, rank_shift, Data);
#else
			h5in.Read3DMatrix(groupName, gdata.om[qin].getName(), Data);
#endif
#else
			h5in.Read3DMatrix(groupName, gdata.om[qin].getName(), Data);
#endif

#if (OMS_USER == TRUE)
		} else {
			gdata.om_user[qin].rename(h5in.GetDatasetName(groupName, q-qshift));
			dataName = gdata.om_user[qin].getName();
#ifdef parallel
#if(HDF_PARALLEL_IO == CRONOS_ON)
			h5in.Read3DMatrix_parallel(groupName, gdata.om_user[qin].getName(),
					mx_data_local, rank_shift, Data);
#else
			h5in.Read3DMatrix(groupName, gdata.om_user[qin].getName(), Data);
#endif
#else
			h5in.Read3DMatrix(groupName, gdata.om_user[qin].getName(), Data);
#endif
		}
#endif

		NumArray<int> staggered(DIM);
		for(int dir=0; dir<DIM; ++dir) {
			staggered(dir) = 0;
		}
		// Now, destinguish different fields
		if(dataName=="A_x") {
			staggered(1) = 1;
			staggered(2) = 1;
		} else if (dataName=="A_y") {
			staggered(0) = 1;
			staggered(2) = 1;
		} else if (dataName=="A_z") {
			staggered(0) = 1;
			staggered(1) = 1;
		}

		// if we do not have a vector potential resize of the magnetic field is not implemented, yet.
		if(levels(0)>0 || levels(1)>0 || levels(2)>0) {
			if((dataName=="B_x" || dataName=="B_y" || dataName=="B_z") && !use_vecpot) {
				cerr << " Error in datain_collective" << endl;
				cerr << " Refinement for magnetic field has NOT been implemented " << endl;
				exit(3);
			}
		}

		// Project data to finer grid
//		cout << " refinement " << Data.getHigh(0) << " " << Data.getHigh(1) << " " << Data.getLow(1)<< endl;
		RefineData(gdata, Data, mx_data_local,levels, staggered, rim_data);

		int imin[DIM];
		int imax[DIM];
		for(int i=0; i<DIM; ++i) {

			imin[i] = -rim_data;
			imax[i] = gdata.mx[i]+rim_data;

		}

		/*
		  if(gdata.mag || (qin <4 || qin>6)) {
		*/
		for(int k=imin[2]; k<=imax[2]; ++k) {
			for(int j=imin[1]; j<=imax[1]; ++j) {
				for(int i=imin[0]; i<=imax[0]; ++i) {

#if (OMS_USER == TRUE)
					if(q < n_omInt+N_SUBS) {
#endif
						gdata.om[qin](i+shift[0],j+shift[1],
						              k+shift[2]) = Data(i+shift[0],j+shift[1],k+shift[2]);
#if (OMS_USER == TRUE)
									} else {
										gdata.om_user[qin](i+shift[0],j+shift[1],
														   k+shift[2]) = Data(i+shift[0],j+shift[1],k+shift[2]);
									}
#endif

				}
			}
		}
//	cout << " rho ";
//	cout << gdata.om[0](0,0,0) << " " << gdata.om[0](1,0,0) << endl;
//	exit(3);
		/*
		  }
		*/
	}


#if(FLUID_TYPE == CRONOS_MULTIFLUID)
	}
#endif

	// Compute factor for timestep adaption
	int max_level = std::max(levels(0), std::max(levels(1), levels(2)));
	double factor = cronos::power(2.,max_level);
//	gdata.dt /= factor;
//	cout << " fac " << factor << " " << gdata.dt << endl;
//	exit(3);

	return factor;

}


void gridFunc::RefineData(Data &gdata, NumMatrix<double,3> &dataOld, int mx[DIM],
		NumArray<int> &factor, NumArray<int> &staggered, int rim) {
	// start with x-dimesion
	NumArray<int> mxOld(3);
	for(int dir=0; dir<DIM; ++dir) {
		mxOld(dir) = mx[dir];
	}

	NumArray<int> mxNew(3);
	mxNew = mxOld;
	NumMatrix<double,3> dataNew;

//	cout << " factors " << factor(0) << " " << factor(1) << endl;

	for(int dir=0; dir<DIM; ++dir) {
		if(factor(dir)>0) {
			for(int level=0; level<factor(dir); ++level) {
				// Compute new x-resolution:
				if(gdata.get_EdgeGridding()) {

					mxNew(dir) = 2*mxOld(dir) + 1;
					dataNew.resize(Index::set(-rim,-rim,-rim),
							Index::set(mxNew[0]+rim,mxNew[1]+rim,mxNew[2]+rim));

//					cout << " Refining " << dataNew.getHigh(0) << " " << dataNew.getHigh(1) << endl;
//					cout << " staggered? " << staggered(dir) << " " << dir << endl;

					// Need to take into account additional grid point in direction of staggering
					NumArray<int> ix_min(3), ix_max(3);
					for(int dir=0; dir<DIM; ++dir) {
						ix_max(dir) = mxOld(dir);
						ix_min(dir) = 0;
						if(staggered(dir)) {
							ix_max(dir) += 1;
						}
					}


//					cout << " Refinement of " << ix_max(0) << " " << ix_max(1) << " " << ix_max(2) << endl;
//					cout << "         onto  " << 2*ix_max(0)+1 << " "<< ix_max(1) << " " << ix_max(2) << " " << endl;

					// interpolate data
					for(int iz=ix_min(2); iz<=ix_max(2); iz++) {
						for(int iy=ix_min(1); iy<=ix_max(1); iy++) {
							for(int ix=ix_min(0); ix<=ix_max(0); ix++) {
								if(dir==0) { // x-dimension
									if(staggered(dir)) {
										// even cells
										dataNew(2*ix  ,iy,iz) = dataOld(ix,iy,iz);
										// odd cells
										dataNew(2*ix+1,iy,iz) = 0.5*dataOld(ix,iy,iz) + 0.5*dataOld(ix+1,iy,iz);
									} else {
										// even cells
										dataNew(2*ix  ,iy,iz) = dataOld(ix,iy,iz);
										// odd cells
										dataNew(2*ix+1,iy,iz) = dataOld(ix,iy,iz);
									}
								} else if (dir==1) { // y-dimension
									if(staggered(dir)) {
										// even cells
										dataNew(ix,2*iy  ,iz) = dataOld(ix,iy,iz);
										// odd cells
										dataNew(ix,2*iy+1,iz) = 0.5*dataOld(ix,iy,iz) + 0.5*dataOld(ix,iy+1,iz);
									} else {
										// even cells
										dataNew(ix,2*iy  ,iz) = dataOld(ix,iy,iz);
										// odd cells
										dataNew(ix,2*iy+1,iz) = dataOld(ix,iy,iz);
									}
								} else { // z-dimension
									if(staggered(dir)) {
										// even cells
										dataNew(ix,iy,2*iz  ) = dataOld(ix,iy,iz);
										// odd cells
										dataNew(ix,iy,2*iz+1) = 0.5*dataOld(ix,iy,iz) + 0.5*dataOld(ix,iy,iz+1);
									} else {
										// even cells
										dataNew(ix,iy,2*iz  ) = dataOld(ix,iy,iz);
										// odd cells
										dataNew(ix,iy,2*iz+1) = dataOld(ix,iy,iz);
									}
								}
							}
						}
					}
//					cout << " after ref " << dataOld.getName() << " " << dataNew(2*mxOld(0),0,0) << endl;
//					cout << " after ref ";
//					cout << mxOld(0) << " " << mxOld(1) << endl;
//					cout << dataOld(-1,0,0) << " " << dataOld(0,0,0) << " " << dataOld(1,0,0) << endl;
//					cout << dataNew(0,0,0) << " " << dataNew(1,0,0) << " " << dataNew(2,0,0) << endl;
					// save re-assigned data
					dataOld = dataNew;
					// use resized grid as new starting point
					mxOld = mxNew;


				} else {
					mxNew(dir) = 2*mxOld(dir);
					cerr << " Well, I have not implemeted this, yet " << endl;
					exit(3);
				}
			}
		}
	}
	// Resized data is stored in dataOld.


}


#if(FLUID_TYPE != CRONOS_MULTIFLUID)
void gridFunc::load_flt(Data &gdata, Hdf5iStream & h5in, string &filename,
                      ProblemType &Problem)
{
	//
	//	This routine is used to load om and om_user data from a hdf5 file
	//
	//	It is assumed that every relevant field is already initialized and ready to go
	// 	The fields will be loaded using the names as set in gdata.
	//
	//	most parts are borrowed from the dataout routine

	int rim = BOUT_FLT;

	// Prepare mx and so on for parallel input if necessary
#ifdef parallel
#if(HDF_PARALLEL_IO == CRONOS_ON)
	NumArray<int> mx_local(3), rank_shift(3);
	for(int dir=0; dir<3; ++dir) {
		mx_local(dir) = gdata.mx[dir];
		rank_shift(dir) = gdata.get_RankShift(dir);
	}
#endif
#endif

	// Buffer
	NumMatrix<float,3> data(Index::set(-rim,-rim,-rim),
	                        Index::set(gdata.mx[0]+rim,gdata.mx[1]+rim,
	                                   gdata.mx[2]+rim));

#if(CRONOS_OUTPUT_COMPATIBILITY == CRONOS_ON)
	string groupName = h5in.get_defaultGroupName();
#else
	string groupName = gdata.fluid.get_Name();
#endif

	int qmin = 0;
	int qmax = n_omInt;

	int qmin_user = 0;
	int qmax_user = n_omIntUser;

	for (int q = qmin; q < qmax; ++q) {
		int qin = q;
		string dsetName = gdata.om[qin].getName();

#ifdef parallel
#if(HDF_PARALLEL_IO == CRONOS_ON)
		// Parallel version
		h5in.Read3DMatrix_withMPI_IO(dsetName, data, mx_local, rank_shift, groupName);
#endif
#else
		// Serial version
		h5in.Read3DMatrix(groupName, dsetName, data);
#endif

		// write it to gdata
		for (int k = -rim; k <= gdata.mx[2]+rim; k++) {
			for (int j = -rim; j <= gdata.mx[1]+rim; j++) {
				for (int i = -rim; i <= gdata.mx[0]+rim; i++) {
					gdata.om[qin](i,j,k) = data(i,j,k);
				}
			}
		}
	}


	// Same for user fields
#if (OMS_USER == TRUE)
	for (int q = qmin_user; q < qmax_user; ++q) {
		int qin = q;
		string dsetName = gdata.om_user[qin].getName();

#ifdef parallel
#if(HDF_PARALLEL_IO == CRONOS_ON)
		// Parallel version
		h5in.Read3DMatrix_withMPI_IO(dsetName, data, mx_local, rank_shift, groupName);
#endif
#else
		// Serial version
		h5in.Read3DMatrix(groupName, dsetName, data);
#endif

		// write it to gdata
		for (int k = -rim; k <= gdata.mx[2]+rim; k++) {
			for (int j = -rim; j <= gdata.mx[1]+rim; j++) {
				for (int i = -rim; i <= gdata.mx[0]+rim; i++) {
					gdata.om_user[qin](i,j,k) = data(i,j,k);
				}
			}
		}
	}
#endif

}
#endif

void gridFunc::Prolongate(Data &gdata, NumMatrix<double,3> &data, const int &q,
                          int mxloc[3], int rim) 
{
	/* 
	   Expanding data to a finer grid. The idea is to keep the center of
	   cells 0 and mx constant.
	*/

	if(data.getLow(0) > -1 || data.getLow(1) > -1 || data.getLow(2) > -1) {
		cerr << " Field has to contain first gridpoint " << endl;
		exit(-22);
	}

	if(data.getHigh(0) < mxloc[0]+1 || data.getHigh(1) < mxloc[1]+1 ||
	   data.getHigh(2) < mxloc[2]+1) {
		cerr << " Field has to contain last gridpoint " << endl;
		exit(-22);
	}

	rim = std::max(rim, 2);

	mxloc[0] *= 2;
	mxloc[1] *= 2;
	mxloc[2] *= 2;
	NumMatrix<double,3> dataloc(Index::set(-rim,-rim,-rim),
	                            Index::set(mxloc[0]+rim,mxloc[1]+rim,mxloc[2]+rim));

	// Setting direct neighbours:
	for(int k=-2; k<=mxloc[2]+2; k+=2) {
		for(int j=-2; j<=mxloc[1]+2; j+=2) {
			for(int i=-2; i<=mxloc[0]+2; i+=2) {		
	
				dataloc(i,j,k) = data(i/2,j/2,k/2);

			}
		}
	}

	data.resize(Index::set(-rim,-rim,-rim),
	            Index::set(mxloc[0]+rim,mxloc[1]+rim,mxloc[2]+rim));

	data = dataloc;

	// All other gridpoints are computed from the above:

	for(int k=-2; k<=mxloc[2]+2; k+=2) {
		for(int j=-2; j<=mxloc[1]+2; j+=2) {
			for(int i=-1; i<=mxloc[0]+1; i+=2) {
	
				data(i,j,k) = 0.5*(data(i+1,j,k) + data(i-1,j,k));

			}
		}
	}

	for(int k=-2; k<=mxloc[2]+2; k+=2) {
		for(int j=-1; j<=mxloc[1]+1; j+=2) {
			for(int i=-2; i<=mxloc[0]+2; i+=2) {

				data(i,j,k) = 0.5*(data(i,j+1,k) + data(i,j-1,k));

			}
		}
	}

	for(int k=-1; k<=mxloc[2]+2; k+=2) {
		for(int j=-2; j<=mxloc[1]+2; j+=2) {
			for(int i=-2; i<=mxloc[0]+2; i+=2) {
	
				data(i,j,k) = 0.5*(data(i,j,k+1) + data(i,j,k-1));

			}
		}
	}



	for(int k=-2; k<=mxloc[2]+2; k+=2) {
		for(int j=-1; j<=mxloc[1]+2; j+=2) {
			for(int i=-1; i<=mxloc[0]+2; i+=2) {
	
				data(i,j,k) = 0.25*(data(i+1,j+1,k) + data(i-1,j+1,k) +
				                    data(i+1,j-1,k) + data(i-1,j-1,k));

			}
		}
	}

	for(int k=-1; k<=mxloc[2]+2; k+=2) {
		for(int j=-2; j<=mxloc[1]+2; j+=2) {
			for(int i=-1; i<=mxloc[0]+2; i+=2) {
	
				data(i,j,k) = 0.25*(data(i+1,j,k+1) + data(i-1,j,k+1) +
				                    data(i+1,j,k-1) + data(i-1,j,k-1));

			}
		}
	}

	for(int k=-1; k<=mxloc[2]+1; k+=2) {
		for(int j=-1; j<=mxloc[1]+1; j+=2) {
			for(int i=-2; i<=mxloc[0]+2; i+=2) {

				data(i,j,k) = 0.25*(data(i,j+1,k+1) + data(i,j-1,k+1) +
				                    data(i,j+1,k-1) + data(i,j-1,k-1));

			}
		}
	}


	for(int k=-1; k<=mxloc[2]+1; k+=2) {
		for(int j=-1; j<=mxloc[1]+1; j+=2) {
			for(int i=-1; i<=mxloc[0]+1; i+=2) {
	
				data(i,j,k) = 0.125*(data(i+1,j+1,k+1) + data(i-1,j+1,k+1) +
				                     data(i+1,j-1,k+1) + data(i-1,j-1,k+1) +
				                     data(i+1,j+1,k-1) + data(i-1,j+1,k-1) +
				                     data(i+1,j-1,k-1) + data(i-1,j-1,k-1));
	
			}
		}
	}
}



#ifdef parallel
void gridFunc::Assign(Data &gdata, NumMatrix<double,3> &data,
                      int mxFile[DIM], int facProc[DIM], 
                      int coordsOld[DIM], int shiftCode[DIM], int rim) 
{
	int fac = 1;
	/*
	  Parameters:
	  mxFile -> Extent of data array from file
	  mxData -> Extent of data array in the code
	*/
  
	int mxData[DIM];
	for(int i=0; i<DIM; ++i) {
		shiftCode[i] = 0;
		fac *= facProc[i];
		if(facProc[i] != 1) {
			if(gdata.get_EdgeGridding()) {
				if(facProc[i] > 1) {
					mxData[i] =  (mxFile[i] + 1)/facProc[i] - 1;
				} else {
					mxData[i] = -(mxFile[i] + 1)*facProc[i] - 1;
				}
			} else {
				if(facProc[i] > 1) {
					mxData[i] =  mxFile[i] / facProc[i];
				} else {
					mxData[i] = -mxFile[i] * facProc[i];
				}
			}
		} else {
			mxData[i] = mxFile[i];
		}
	}

	if(abs(fac) == 1) {
		return;
	} else {

		NumMatrix<double,3> dataLoc;
		dataLoc.resize(Index::set(-rim, -rim, -rim),
		               Index::set(mxData[0]+rim, mxData[1]+rim, mxData[2]+rim));

		// Getting the shift:
		int shiftFile[DIM];
		int imax[DIM];
		for(int i=0; i<DIM; ++i) {
			shiftFile[i] = 0;
			if(facProc[i] > 1) {
				if(gdata.get_EdgeGridding()) {
					shiftFile[i] = (gdata.coords[i]%facProc[i])*(mxData[i]+1);
				} else {
					shiftFile[i] = (gdata.coords[i]%facProc[i])*mxData[i];
				}
				imax[i] = mxData[i];
			} else {
				if(gdata.get_EdgeGridding()) {
					shiftCode[i] = (coordsOld[i]%(-facProc[i]))*(mxFile[i]+1);
				} else {
					shiftCode[i] = (coordsOld[i]%(-facProc[i]))*mxFile[i];
				}
				imax[i] = mxFile[i];
			}
		}

		for(int i=-rim; i<=imax[0]+rim; ++i) {
			for(int j=-rim; j<=imax[1]+rim; ++j) {
				for(int k=-rim; k<=imax[2]+rim; ++k) {
	  
					dataLoc(i+shiftCode[0],j+shiftCode[1],
					        k+shiftCode[2]) = data(i+shiftFile[0],j+shiftFile[1],
					                               k+shiftFile[2]);

				}
			}
		}
	      
		data.resize(Index::set(-rim, -rim, -rim),
		            Index::set(mxData[0]+rim, mxData[1]+rim, mxData[2]+rim));

		data = dataLoc;
	}
	return;
  
}
#endif


