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

	// Use a value that is not used (-1 used by user fields -> so, I am using -2)
	q_Bx = -2;
	q_By = -2;
	q_Bz = -2;

	n_om = N_OM;
	n_omInt = N_OMINT;
#if (OMS_USER == TRUE)
	n_omUser = N_OM_USER;
	n_omIntUser = N_OMINT_USER;
#else
	n_omIntUser = 0;
#endif
	n_omIntAll = N_OMINT_ALL;

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

	// if(atAxis) {
#if (GEOM == CYLINDRICAL)
	double eps(1.e-5*(gdata.global_xe[0]-gdata.global_xb[0]));
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
		double eps(1.e-5*(gdata.global_xe[1]-gdata.global_xb[1]));
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
		double eps(0.);
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
	CronosFluid fluid = gdata.fluid;
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
		q = Problem.get_q(omb);
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

void gridFunc::boundary(Queue &queue, Data &gdata, ProblemType &Problem,
                        NumMatrix<double,3> &omb, int rim, int q, int iFluid)
{

	//-----------------------------------------------------------
	// x-direction
	//-----------------------------------------------------------

	cout << "q_bctype: " << q << " " << bc_Type[q] << endl;

	bc_Periodic(queue, gdata, Problem, omb, 0, q, rim);

	if(bc_Type[0] > 1) {
		if (bc_Type[0] == 2) {
			bc_Extrapolate(queue, gdata, Problem, omb,0,0,q,rim);
		} else if (bc_Type[0] == 3) {
			bc_Outflow(queue, gdata, Problem, omb,0,0,q,rim);
		} else if (bc_Type[0] == 4) {
			Problem.bc_User(queue, gdata, omb,0,0,q,rim);
		} else if (bc_Type[0] == 5) {
			cerr << " Not implemented" << endl;
		} else if (bc_Type[0] == 6) {
			cerr << " Not implemented" << endl;
		} else if (bc_Type[0] == 7) {
			cerr << " Not implemented" << endl;
		} else {
			cerr << " No such boundary condition 1 " << endl;
		}
	}

	// Right:
	if(bc_Type[1] > 1) {
		if (bc_Type[1] == 2){
			bc_Extrapolate(queue, gdata, Problem, omb,0,1,q,rim);
		} else if (bc_Type[1] == 3) {
			bc_Outflow(queue, gdata, Problem, omb,0,1,q,rim);
		} else if (bc_Type[1] == 4) {
			Problem.bc_User(queue, gdata, omb,0,1,q,rim);
		} else if (bc_Type[1] == 5) {
			cerr << " Not implemented" << endl;
		} else if (bc_Type[1] == 7) {
			cerr << " Not implemented" << endl;
		} else {
			cerr << " No such boundary condition 2 " << endl;
		}
	}

	//-----------------------------------------------------------
	// y-direction
	//-----------------------------------------------------------

	bc_Periodic(queue, gdata, Problem, omb, 1, q, rim);

	if(bc_Type[2] > 1) {
		if(bc_Type[2] == 2) {
			bc_Extrapolate(queue, gdata, Problem, omb,1,0,q,rim);
		} else if(bc_Type[2] == 3) {
			bc_Outflow(queue, gdata, Problem, omb,1,0,q,rim);
		} else if(bc_Type[2] == 4) {
			Problem.bc_User(queue, gdata, omb,1,0,q,rim);
		} else if(bc_Type[2] == 5) {
			cerr << " Not implemented" << endl;
		} else if (bc_Type[2] == 6) {
			cerr << " Not implemented" << endl;
		} else if (bc_Type[2] == 7) {
			cerr << " Not implemented" << endl;
		} else {
			cerr << " No such boundary condition 3 " << " " << bc_Type[1]<< endl;
		}
	}

	if(bc_Type[3] > 1) {
		if(bc_Type[3] == 2) {
			bc_Extrapolate(queue, gdata, Problem, omb,1,1,q,rim);
		} else if(bc_Type[3] == 3) {
			bc_Outflow(queue, gdata, Problem, omb,1,1,q,rim);
		} else if(bc_Type[3] == 4) {
			Problem.bc_User(queue, gdata, omb,1,1,q,rim);
		} else if(bc_Type[3] == 5) {
			cerr << " Not implemented" << endl;
		} else if (bc_Type[3] == 6) {
			cerr << " Not implemented" << endl;
		} else if (bc_Type[3] == 7) {
			cerr << " Not implemented" << endl;
		} else {
			cerr << " No such boundary condition 4 " << endl;
		}
	}

	//-----------------------------------------------------------
	// z-direction
	//-----------------------------------------------------------

	bc_Periodic(queue, gdata, Problem, omb, 2, q, rim);

	if(bc_Type[4] > 1) {
		if (bc_Type[4] == 2) {
			bc_Extrapolate(queue, gdata, Problem, omb,2,0,q,rim);
		} else if (bc_Type[4] == 3) {
			bc_Outflow(queue, gdata, Problem, omb,2,0,q,rim);
		} else if (bc_Type[4] == 4) {
			Problem.bc_User(queue, gdata, omb,2,0,q,rim);
		} else if (bc_Type[4] == 5) {
			cerr << " Not implemented" << endl;
		} else if (bc_Type[4] == 7) {
			cerr << " Not implemented" << endl;
		} else {
			cerr << " No such boundary condition 5 " << endl;
		}
	}

	if(bc_Type[5] > 1) {
		if (bc_Type[5] == 2) {
			bc_Extrapolate(queue, gdata, Problem, omb,2,1,q,rim);
		} else if (bc_Type[5] == 3) {
			bc_Outflow(queue, gdata, Problem, omb,2,1,q,rim);
		} else if (bc_Type[5] == 4) {
			Problem.bc_User(queue, gdata, omb,2,1,q,rim);
		} else if (bc_Type[5] == 5) {
			cerr << " Not implemented" << endl;
		} else if (bc_Type[5] == 7) {
			cerr << " Not implemented" << endl;
		} else {
			cerr << " No such boundary condition 6 " << endl;
		}
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
		q = Problem.get_q(omb);
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
	bc_Periodic_serial(gdata, Problem, omb, 0, q, rim);

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

	// Non-MPI periodic boundary conditions
	bc_Periodic_serial(gdata, Problem, omb, 1, q, rim);

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

	// Non-MPI periodic boundary conditions
	bc_Periodic_serial(gdata, Problem, omb, 2, q, rim);

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

void gridFunc::bc_SendRecv(Data &gdata, NumMatrix<double,3> &Send,
                           NumMatrix<double,3> &Recv, int dir, int above,
                           int q, bool do_Send, bool do_Receive)
{
	// (Nearly) Nothing to be done
	Recv = Send;
}


void gridFunc::bc_SendRecv(Data &gdata, NumArray<double> &Send,
		NumArray<double> &Recv, int dir, int size_send, int size_recv, int above, int q, bool do_Send, bool do_Receive)
{
	// (Nearly) Nothing to be done
	for(int iter=0; iter<size_send; iter++) {
		Recv(iter) = Send(iter);
	}
}



void gridFunc::bc_ISendRecv(Data &gdata, NumArray<double> &Send, NumArray<double> &Recv,
		int dir, int size_send, int size_recv, int above, int q, bool do_Send, bool do_Receive)
{
	// (Nearly) Nothing to be done
	for(int iter=0; iter<size_send; iter++) {
		Recv(iter) = Send(iter);
	}
}



void gridFunc::bc_Periodic_old(Data &gdata, ProblemType &Problem,
                           NumMatrix<double,3> &omb,
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
		SendLeft = false;
		RecvRight = false;
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

void gridFunc::bc_Periodic(Queue &queue, Data &gdata, ProblemType &Problem,
		NumMatrix<double,3> &omb, int dir, int q, int rim) {

	int shift(0);
	if(gdata.get_EdgeGridding()) {
		shift = 1;
	}
	auto omRange = gdata.omSYCL[q].get_range();
	int i_min, i_max;

	if((bc_Type[2*dir] < 2) && (bc_Type[2*dir+1] < 2)) {

		if (dir == 0) {

			get_bcBuffRange(gdata, i_min, i_max, 0, 1, q, rim);
			auto left_to_right = Range<3>(i_max + 1 - i_min, omRange.get(1), omRange.get(2));
			queue.submit(celerity::allow_by_ref, [=, &gdata](celerity::handler& cgh) {
				celerity::accessor om_acc{gdata.omSYCL[q], cgh, celerity::access::neighborhood{(size_t) gdata.mx[0] + 2*rim + 1,1,1}, celerity::read_write};
				cgh.parallel_for<class Periodic_Dir1_Left_to_Right>(left_to_right, [=, &gdata](celerity::item<3> item){
					size_t ix = item.get_id(0) + i_min;
					size_t iy = item.get_id(1);
					size_t iz = item.get_id(2);

					om_acc[gdata.mx[0]+rim+ix][iy][iz] = om_acc[ix+rim-shift][iy][iz];
				});
			});

			get_bcBuffRange(gdata, i_min, i_max, 0, 0, q, rim);
			auto right_to_left = Range<3>(-(i_min + 1) - i_max, omRange.get(1), omRange.get(2));
			queue.submit(celerity::allow_by_ref, [=, &gdata](celerity::handler& cgh) {
				celerity::accessor om_acc{gdata.omSYCL[q], cgh, celerity::access::neighborhood{(size_t) gdata.mx[0] + 2*rim + 1,1,1}, celerity::read_write};
				cgh.parallel_for<class Periodic_Dir1_Right_to_Left>(right_to_left, [=, &gdata](celerity::item<3> item){
					size_t ix = item.get_id(0);
					size_t iy = item.get_id(1);
					size_t iz = item.get_id(2);

					om_acc[ix][iy][iz] = om_acc[gdata.mx[0]+rim+ix+shift][iy][iz];
				});
			});

		} else if (dir == 2) {

			get_bcBuffRange(gdata, i_min, i_max, 1, 1, q, rim);
			auto left_to_right = Range<3>(omRange.get(0), i_max + 1 - i_min, omRange.get(2));
			queue.submit(celerity::allow_by_ref, [=, &gdata](celerity::handler& cgh) {
				celerity::accessor om_acc{gdata.omSYCL[q], cgh, celerity::access::neighborhood{1,(size_t) gdata.mx[1] + 2*rim + 1,1}, celerity::read_write};
				cgh.parallel_for<class Periodic_Dir1_Left_to_Right>(left_to_right, [=, &gdata](celerity::item<3> item){
					size_t ix = item.get_id(0);
					size_t iy = item.get_id(1) + i_min;
					size_t iz = item.get_id(2);

					om_acc[ix][gdata.mx[1]+rim+iy][iz] = om_acc[ix][iy+rim-shift][iz];
				});
			});

			get_bcBuffRange(gdata, i_min, i_max, 1, 0, q, rim);
			auto right_to_left = Range<3>(omRange.get(0), -(i_min + 1) - i_max, omRange.get(2));
			queue.submit(celerity::allow_by_ref, [=, &gdata](celerity::handler& cgh) {
				celerity::accessor om_acc{gdata.omSYCL[q], cgh, celerity::access::neighborhood{1,(size_t) gdata.mx[1] + 2*rim + 1,1}, celerity::read_write};
				cgh.parallel_for<class Periodic_Dir1_Right_to_Left>(right_to_left, [=, &gdata](celerity::item<3> item){
					size_t ix = item.get_id(0);
					size_t iy = item.get_id(1);
					size_t iz = item.get_id(2);

					om_acc[ix][iy][iz] = om_acc[ix][gdata.mx[1]+rim+iy+shift][iz];
				});
			});

		} else if (dir == 3) {

			get_bcBuffRange(gdata, i_min, i_max, 2, 1, q, rim);
			auto left_to_right = Range<3>(omRange.get(0), omRange.get(1), i_max + 1 - i_min);
			queue.submit(celerity::allow_by_ref, [=, &gdata](celerity::handler& cgh) {
				celerity::accessor om_acc{gdata.omSYCL[q], cgh, celerity::access::neighborhood{1,1,(size_t) gdata.mx[2] + 2*rim + 1}, celerity::read_write};
				cgh.parallel_for<class Periodic_Dir1_Left_to_Right>(left_to_right, [=, &gdata](celerity::item<3> item){
					size_t ix = item.get_id(0);
					size_t iy = item.get_id(1);
					size_t iz = item.get_id(2) + i_min;

					om_acc[ix][iy][gdata.mx[2]+rim+iz] = om_acc[ix][iy][iz+rim-shift];
				});
			});

			get_bcBuffRange(gdata, i_min, i_max, 2, 0, q, rim);
			auto right_to_left = Range<3>(omRange.get(0), omRange.get(1), -(i_min + 1) - i_max);
			queue.submit(celerity::allow_by_ref, [=, &gdata](celerity::handler& cgh) {
				celerity::accessor om_acc{gdata.omSYCL[q], cgh, celerity::access::neighborhood{1,1,(size_t) gdata.mx[2] + 2*rim + 1}, celerity::read_write};
				cgh.parallel_for<class Periodic_Dir1_Right_to_Left>(right_to_left, [=, &gdata](celerity::item<3> item){
					size_t ix = item.get_id(0);
					size_t iy = item.get_id(1);
					size_t iz = item.get_id(2);

					om_acc[ix][iy][iz] = om_acc[ix][iy][gdata.mx[2]+rim+iz+shift];
				});
			});
		}
	} 
}


void gridFunc::bc_Periodic_serial(Data &gdata, ProblemType &Problem,
		NumMatrix<double,3> &omb, int dir, int q, int rim) {

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


void gridFunc::bc_Periodic(Data &gdata, ProblemType &Problem,
                           NumMatrix<double,3> &omb,
                           int dir, int q, int rim, bool internal_only) {



	bool SendLeft(false), RecvRight(false);

	if(bc_Type[2*dir] < 2) SendLeft = true;
	if(bc_Type[2*dir+1] < 2) RecvRight = true;

	if(internal_only) {
		SendLeft = false;
		RecvRight = false;
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


void gridFunc::ExchangePositions(Data &gdata, int top) {
	// First we need to find if we have a corresponding opposite rank
	// (should moslty be the case - only when deviding by odd number
	// would we expect problems)

	// cerr << " Exchanging " << gdata.rank << " " << gdata.coords[0] << endl;

#if (GEOM != CARTESIAN)
	PhiPos[top].resize(Index::set(-B),
	                   Index::set(gdata.global_mx[gdata.phi_dir]+B));

	// Compute Phi-positions also for serial run
	for(int ip = -B; ip <= gdata.mx[gdata.phi_dir]+B; ++ip) { 
		double phi = gdata.getCen(gdata.phi_dir, ip);
		
		PhiPos[top](ip) = phi;
		
	}

#endif // not cartesian

}


void gridFunc::compute_AxisPartners(Data &gdata, int top) {
	// Computing the partners for axis boundary conditions - so far
	// only for cylindrical coordinates

#if (GEOM != CARTESIAN)
	ExchangePositions(gdata, top);

	// Array over full phi-domain
	AxisPartnersPhi[top].resize(Index::set(-B),
			Index::set(gdata.mx[gdata.phi_dir]+B));

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
	int i_phi_other = AxisPartnersPhi[top](iz+shift);
	return i_phi_other;
}




void gridFunc::bc_Axis(Data &gdata, ProblemType &Problem,
		NumMatrix<double,3> &omb, int dir, int top, int q, int rim) {
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
		NumMatrix<double,3> &omb, int q, int rim) {

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
	NumMatrix<double,3> & omb_global = omb;

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
		NumMatrix<double,3> &omb, int top, int q, int rim) {
#if (GEOM == 3)

	NumMatrix<double,3> & omb_global = omb;

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
		double Bcx_sum = 0.;  // new values for Cartesian Bx, By
		double Bcy_sum = 0.;  //  on the singular axis
		int N_terms(0);
			//		for (int j = 0; j <= gdata.mx[1]; j++) {
			// Do not use phi=0 twice!
		if(gdata.get_singularity_treatment(0) == 1) { // at lower r-bound
			for (int iy = 0; iy < jnum; iy++) {
				// br0_j, bp0_j: projected to (r,phi)=(0,phi_j) via avg.
				double br0_j = (+gdata.om[q_Bx]( 0,iy,  iz)
				              +gdata.om[q_Bx](-2,iy,  iz))/2.;
				double bp0_j = (+gdata.om[q_By]( 0,iy-1,iz)
				              +gdata.om[q_By]( 0,iy  ,iz)
				              +gdata.om[q_By](-1,iy-1,iz)
				              +gdata.om[q_By](-1,iy  ,iz))/4.;
				double phi = gdata.getCen(1, iy);
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

		double Bcx_ave(Bcx_sum/N_terms);
		double Bcy_ave(Bcy_sum/N_terms);

		// if((gdata.rank==0 && k==gdata.mx[2]) ||
		//    (gdata.rank==1 && k==-1)) {
		// 	cout << setiosflags(ios::scientific) << setprecision(12);
		// 	cout << " corr val: " << gdata.rank << " " << Bcx_ave << " " <<	Bcy_ave << endl;
		// 	cout << resetiosflags(ios::scientific) << setprecision(6);
		// }

		if(gdata.get_singularity_treatment(0) == 1) { // at lower r-bound
			for (int iy = -B; iy <= gdata.mx[1]+B; iy++) { // set Br(0)s
				double phi = gdata.getCen(1, iy);
				double Br0_ave = ( Bcx_ave * cos(phi) +
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
		double cost0 = 1.-2.*top;      // cos({0|Pi}) = {+1|-1}
		int jS = top*(gdata.mx[1]+1)-1; // {-1|mxtet} (sing. index)
		for (int ix = -B; ix <= gdata.mx[0]+B; ++ix) {
			double Bcx_sum = 0.;  // new values for Cartesian Bx, By
			double Bcy_sum = 0.;  //  on the singular axis
			for (int iz = 0; iz <= gdata.mx[2]; ++iz) {
				// bt0_k, bp0_k: projected to (tet,phi)=(0,phi_k) via avg.
				double bt0_k = (+gdata.om[q_By](ix,jS-1, iz  )
						+gdata.om[q_By](ix,jS+1, iz  ))*0.5;
				double bp0_k = (+gdata.om[q_Bz](ix,jS  , iz-1)
						+gdata.om[q_Bz](ix,jS  , iz  )
						+gdata.om[q_Bz](ix,jS+1, iz-1)
						+gdata.om[q_Bz](ix,jS+1, iz  ))*0.25;

				Bcx_sum += bt0_k * cost0 * cosPos[2](iz) - bp0_k * sinPos[2](iz);
				Bcy_sum += bt0_k * cost0 * sinPos[2](iz) + bp0_k * cosPos[2](iz);
			}
			int N_terms = gdata.mx[2]+1;

			double Bcx_ave(Bcx_sum/N_terms);
			double Bcy_ave(Bcy_sum/N_terms);

			for (int iz = -B; iz <= gdata.mx[2]+B; iz++) { // set Bt(0)s
				double Bt0_ave = ( Bcx_ave * cost0 * cosPos[2](iz) +
						Bcy_ave * cost0 * sinPos[2](iz) );
				//         NB:  -Bcz_sum * sint0 ~ sin({0|Pi}) = 0
				gdata.om[q_By](ix,jS,iz) = Bt0_ave;
				// gdata.om[q_By](i,jS,k) = (gdata.om[q_By](i,jS-1,k)+
				//                           gdata.om[q_By](i,jS+1,k))/2.;
			} // k
		} // i
	}

//#endif
}

void gridFunc::bc_Extrapolate(Queue &queue, Data &gdata, ProblemType &Problem, 
                              NumMatrix<double,3> &omb,
                              int dir, int above, int q, int rim)
{

	int mx[3] = {gdata.mx[0], gdata.mx[1], gdata.mx[2]};

	if(dir == 0) {
		if(above == 0) {
			int ixmin = 0;

			auto omRange = gdata.omSYCL[q].get_range();
			auto range = Range<3>(rim + ixmin, omRange.get(1), omRange.get(2));

			queue.submit(celerity::allow_by_ref, [=, &gdata](celerity::handler& cgh) {

				celerity::accessor om_acc{gdata.omSYCL[q], cgh, celerity::access::neighborhood{(size_t) rim+1,1,1}, celerity::read_write};

				cgh.parallel_for<class IntegrationKernel>(range, [=](celerity::item<3> item){

					size_t ix = item.get_id(0);
					size_t iy = item.get_id(1);
					size_t iz = item.get_id(2);

					om_acc[ix][iy][iz] = - om_acc[ix + 1][iy][iz];

				});
			});

		} else if (above == 1) {
			int ixmin = 1;

			auto omRange = gdata.omSYCL[q].get_range();
			auto range = Range<3>(rim + ixmin, omRange.get(1), omRange.get(2));

			queue.submit(celerity::allow_by_ref, [=, &gdata](celerity::handler& cgh) {

				celerity::accessor om_acc{gdata.omSYCL[q], cgh, celerity::access::neighborhood{(size_t) rim+1,1,1}, celerity::read_write};

				cgh.parallel_for<class IntegrationKernel>(range, [=](celerity::item<3> item){

					size_t ix = item.get_id(0) + ixmin;
					size_t iy = item.get_id(1);
					size_t iz = item.get_id(2);

					om_acc[mx[0]+rim+ix][iy][iz] = - om_acc[mx[0]+rim+ix-1][iy][iz];

				});
			});
		}
	}

	if(dir == 1) {
		if(above == 0) {
			int iymin = 0;

			auto omRange = gdata.omSYCL[q].get_range();
			auto range = Range<3>(omRange.get(0), rim + iymin, omRange.get(2));

			queue.submit(celerity::allow_by_ref, [=, &gdata](celerity::handler& cgh) {

				celerity::accessor om_acc{gdata.omSYCL[q], cgh, celerity::access::neighborhood{1,(size_t) rim+1,1}, celerity::read_write};

				cgh.parallel_for<class IntegrationKernel>(range, [=](celerity::item<3> item){

					size_t ix = item.get_id(0);
					size_t iy = item.get_id(1);
					size_t iz = item.get_id(2);

					om_acc[ix][iy][iz] = - om_acc[ix][iy + 1][iz];

				});
			});

		} else if (above == 1) {
			int iymin = 1;

			auto omRange = gdata.omSYCL[q].get_range();
			auto range = Range<3>(omRange.get(0), rim + iymin, omRange.get(2));

			queue.submit(celerity::allow_by_ref, [=, &gdata](celerity::handler& cgh) {

				celerity::accessor om_acc{gdata.omSYCL[q], cgh, celerity::access::neighborhood{1,(size_t) rim+1,1}, celerity::read_write};

				cgh.parallel_for<class IntegrationKernel>(range, [=](celerity::item<3> item){

					size_t ix = item.get_id(0);
					size_t iy = item.get_id(1) + iymin;
					size_t iz = item.get_id(2);

					om_acc[ix][mx[1]+rim+iy][iz] = - om_acc[ix][mx[1]+rim+iy-1][iz];

				});
			});
		}
	}

	if(dir == 1) {
		if(above == 0) {
			int izmin = 0;

			auto omRange = gdata.omSYCL[q].get_range();
			auto range = Range<3>(omRange.get(0), omRange.get(1), rim + izmin);

			queue.submit(celerity::allow_by_ref, [=, &gdata](celerity::handler& cgh) {

				celerity::accessor om_acc{gdata.omSYCL[q], cgh, celerity::access::neighborhood{1,1,(size_t) rim+1}, celerity::read_write};

				cgh.parallel_for<class IntegrationKernel>(range, [=](celerity::item<3> item){

					size_t ix = item.get_id(0);
					size_t iy = item.get_id(1);
					size_t iz = item.get_id(2);

					om_acc[ix][iy][iz] = - om_acc[ix][iy][iz + 1];

				});
			});

		} else if (above == 1) {
			int izmin = 1;

			auto omRange = gdata.omSYCL[q].get_range();
			auto range = Range<3>(omRange.get(0), omRange.get(1), rim + izmin);

			queue.submit(celerity::allow_by_ref, [=, &gdata](celerity::handler& cgh) {

				celerity::accessor om_acc{gdata.omSYCL[q], cgh, celerity::access::neighborhood{1,1,(size_t) rim+1}, celerity::read_write};

				cgh.parallel_for<class IntegrationKernel>(range, [=](celerity::item<3> item){

					size_t ix = item.get_id(0);
					size_t iy = item.get_id(1);
					size_t iz = item.get_id(2) + izmin;

					om_acc[ix][iy][mx[2]+rim+iz] = - om_acc[ix][iy][mx[2]+rim+iz-1];

				});
			});
		}
	}
}

void gridFunc::bc_Extrapolate(Data &gdata, ProblemType &Problem, 
                              NumMatrix<double,3> &omb,
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


void gridFunc::bc_Outflow(Queue &queue, Data &gdata, ProblemType &pr,
                          NumMatrix<double,3> &omb,
                          int dir, int above, int q, int rim)

{
	if(dir == 0) {
		if(above == 0) {
			
			auto omRange = gdata.omSYCL[q].get_range();
			auto range = Range<3>(rim, omRange.get(1), omRange.get(2));

			queue.submit(celerity::allow_by_ref, [=, &gdata](celerity::handler& cgh) {

				celerity::accessor om_acc{gdata.omSYCL[q], cgh, celerity::access::neighborhood{(size_t) rim,1,1}, celerity::read_write};
				celerity::accessor om_sx_acc{gdata.omSYCL[1], cgh, celerity::access::one_to_one{}, celerity::read_only};

				cgh.parallel_for<class BCOutflowKernel>(range, [=](celerity::item<3> item){

					size_t ix = item.get_id(0);
					size_t iy = item.get_id(1);
					size_t iz = item.get_id(2);

					if (om_sx_acc[rim][iy][iz] < 0.) {
						om_acc[ix][iy][iz] = om_acc[ix + 1][iy][iz];
					} else {
						if (q == 1) {
							om_acc[ix][iy][iz] = - om_acc[rim - ix - 1][iy][iz];
						} else {
							om_acc[ix][iy][iz] = om_acc[rim - ix - 1][iy][iz];
						}
					}

				});
			});
		} else {
			auto omRange = gdata.omSYCL[q].get_range();
			auto range = Range<3>(rim, omRange.get(1), omRange.get(2));

			queue.submit(celerity::allow_by_ref, [=, &gdata](celerity::handler& cgh) {

				celerity::accessor om_acc{gdata.omSYCL[q], cgh, celerity::access::neighborhood{(size_t) rim,1,1}, celerity::read_write};
				celerity::accessor om_sx_acc{gdata.omSYCL[1], cgh, celerity::access::one_to_one{}, celerity::read_only};

				cgh.parallel_for<class BCOutflowKernel>(range, [=, &gdata](celerity::item<3> item){

					size_t ix = item.get_id(0) + 1;
					size_t iy = item.get_id(1);
					size_t iz = item.get_id(2);

					if (om_sx_acc[gdata.mx[0] + rim][iy][iz] > 0.) {
						om_acc[gdata.mx[0] + rim + ix][iy][iz] = om_acc[gdata.mx[0] + rim + ix - 1][iy][iz];
					} else {
						if (q == 1) {
							om_acc[gdata.mx[0] + rim + ix][iy][iz] = - om_acc[gdata.mx[0] + rim - ix + 1][iy][iz];
						} else {
							om_acc[gdata.mx[0] + rim + ix][iy][iz] = om_acc[gdata.mx[0] + rim - ix + 1][iy][iz];
						}
					}

				});
			});
		}
	}

	if(dir == 1) {
		if(above == 0) {
			
			auto omRange = gdata.omSYCL[q].get_range();
			auto range = Range<3>(omRange.get(0), rim, omRange.get(2));

			queue.submit(celerity::allow_by_ref, [=, &gdata](celerity::handler& cgh) {

				celerity::accessor om_acc{gdata.omSYCL[q], cgh, celerity::access::neighborhood{1,(size_t) rim,1}, celerity::read_write};
				celerity::accessor om_sy_acc{gdata.omSYCL[2], cgh, celerity::access::one_to_one{}, celerity::read_only};

				cgh.parallel_for<class BCOutflowKernel>(range, [=](celerity::item<3> item){

					size_t ix = item.get_id(0);
					size_t iy = item.get_id(1);
					size_t iz = item.get_id(2);

					if (om_sy_acc[ix][rim][iz] < 0.) {
						om_acc[ix][iy][iz] = om_acc[ix][iy + 1][iz];
					} else {
						if (q == 2) {
							om_acc[ix][iy][iz] = - om_acc[ix][rim - iy - 1][iz];
						} else {
							om_acc[ix][iy][iz] = om_acc[ix][rim - iy - 1][iz];
						}
					}

				});
			});
		} else {
			auto omRange = gdata.omSYCL[q].get_range();
			auto range = Range<3>(omRange.get(0), rim, omRange.get(2));

			queue.submit(celerity::allow_by_ref, [=, &gdata](celerity::handler& cgh) {

				celerity::accessor om_acc{gdata.omSYCL[q], cgh, celerity::access::neighborhood{1,(size_t) rim,1}, celerity::read_write};
				celerity::accessor om_sy_acc{gdata.omSYCL[2], cgh, celerity::access::one_to_one{}, celerity::read_only};

				cgh.parallel_for<class BCOutflowKernel>(range, [=, &gdata](celerity::item<3> item){

					size_t ix = item.get_id(0);
					size_t iy = item.get_id(1) + 1;
					size_t iz = item.get_id(2);

					if (om_sy_acc[ix][gdata.mx[1] + rim][iz] > 0.) {
						om_acc[ix][gdata.mx[1] + rim + iy][iz] = om_acc[ix][gdata.mx[1] + rim + iy - 1][iz];
					} else {
						if (q == 2) {
							om_acc[ix][gdata.mx[1] + rim + iy][iz] = - om_acc[ix][gdata.mx[1] + rim - iy + 1][iz];
						} else {
							om_acc[ix][gdata.mx[1] + rim + iy][iz] = om_acc[ix][gdata.mx[1] + rim - iy + 1][iz];
						}
					}

				});
			});
		}
	}

	if(dir == 2) {
		if(above == 0) {
			
			auto omRange = gdata.omSYCL[q].get_range();
			auto range = Range<3>(omRange.get(0), omRange.get(1), rim);

			queue.submit(celerity::allow_by_ref, [=, &gdata](celerity::handler& cgh) {

				celerity::accessor om_acc{gdata.omSYCL[q], cgh, celerity::access::neighborhood{1,1,(size_t) rim}, celerity::read_write};
				celerity::accessor om_sz_acc{gdata.omSYCL[3], cgh, celerity::access::one_to_one{}, celerity::read_only};

				cgh.parallel_for<class BCOutflowKernel>(range, [=](celerity::item<3> item){

					size_t ix = item.get_id(0);
					size_t iy = item.get_id(1);
					size_t iz = item.get_id(2);

					if (om_sz_acc[ix][iy][rim] < 0.) {
						om_acc[ix][iy][iz] = om_acc[ix][iy][iz + 1];
					} else {
						if (q == 2) {
							om_acc[ix][iy][iz] = - om_acc[ix][iy][rim - iz - 1];
						} else {
							om_acc[ix][iy][iz] = om_acc[ix][iy][rim - iz - 1];
						}
					}

				});
			});
		} else {
			auto omRange = gdata.omSYCL[q].get_range();
			auto range = Range<3>(omRange.get(0), omRange.get(1), rim);

			queue.submit(celerity::allow_by_ref, [=, &gdata](celerity::handler& cgh) {

				celerity::accessor om_acc{gdata.omSYCL[q], cgh, celerity::access::neighborhood{1,1,(size_t) rim}, celerity::read_write};
				celerity::accessor om_sz_acc{gdata.omSYCL[3], cgh, celerity::access::one_to_one{}, celerity::read_only};

				cgh.parallel_for<class BCOutflowKernel>(range, [=, &gdata](celerity::item<3> item){

					size_t ix = item.get_id(0);
					size_t iy = item.get_id(1);
					size_t iz = item.get_id(2) + 1;

					if (om_sz_acc[ix][iy][gdata.mx[2] + rim] > 0.) {
						om_acc[ix][iy][gdata.mx[2] + rim + iz] = om_acc[ix][iy][gdata.mx[2] + rim + iz - 1];
					} else {
						if (q == 3) {
							om_acc[ix][iy][gdata.mx[2] + rim + iz] = - om_acc[ix][iy][gdata.mx[2] + rim - iz + 1];
						} else {
							om_acc[ix][iy][gdata.mx[2] + rim + iz] = om_acc[ix][iy][gdata.mx[2] + rim - iz + 1];
						}
					}

				});
			});
		}
	}


}



void gridFunc::bc_Outflow(Data &gdata, ProblemType &pr,
                          NumMatrix<double,3> &omb,
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
                             NumMatrix<double,3> &omb,
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
                        NumMatrix<double,3> &omb,
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
				NumMatrix<double,3> & bcVals = bcVals_xb[q];
				// Doing bcs:
				bc_Fixed_general(omb, bcVals, imin, imax);
#if (OMS_USER == TRUE)
			} else {
				NumMatrix<double,3> & bcVals = bcVals_User_xb[q];
				// Doing bcs:
				bc_Fixed_general(omb, bcVals, imin, imax);
#endif
			}

		} else {

			int imin[3] = {gdata.mx[0]+1, -rim, -rim};
			int imax[3] = {gdata.mx[0]+rim, gdata.mx[1]+rim, gdata.mx[2]+rim};

			if(generic_vars) {
				NumMatrix<double,3> & bcVals = bcVals_xe[q];
				// Doing bcs:
				bc_Fixed_general(omb, bcVals, imin, imax);
#if (OMS_USER == TRUE)
			} else {
				NumMatrix<double,3> & bcVals = bcVals_User_xe[q];
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
				NumMatrix<double,3> & bcVals = bcVals_yb[q];
				// Doing bcs:
				bc_Fixed_general(omb, bcVals, imin, imax);
#if (OMS_USER == TRUE)
			} else {
				NumMatrix<double,3> & bcVals = bcVals_User_yb[q];
				// Doing bcs:
				bc_Fixed_general(omb, bcVals, imin, imax);
#endif
			}

		} else {

			int imin[3] = {-rim, gdata.mx[1]+1, -rim};
			int imax[3] = {gdata.mx[0]+rim, gdata.mx[1]+rim, gdata.mx[2]+rim};

			if(generic_vars) {
				NumMatrix<double,3> & bcVals = bcVals_ye[q];
				// Doing bcs:
				bc_Fixed_general(omb, bcVals, imin, imax);
#if (OMS_USER == TRUE)
			} else {
				NumMatrix<double,3> & bcVals = bcVals_User_ye[q];
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
				NumMatrix<double,3> & bcVals = bcVals_zb[q];
				// Doing bcs:
				bc_Fixed_general(omb, bcVals, imin, imax);
#if (OMS_USER == TRUE)
			} else {
				NumMatrix<double,3> & bcVals = bcVals_User_zb[q];
				// Doing bcs:
				bc_Fixed_general(omb, bcVals, imin, imax);
#endif
			}

		} else {

			int imin[3] = {-rim, -rim, gdata.mx[2]+1};
			int imax[3] = {gdata.mx[0]+rim, gdata.mx[1]+rim, gdata.mx[2]+rim};

			if(generic_vars) {
				NumMatrix<double,3> & bcVals = bcVals_ze[q];
				// Doing bcs:
				bc_Fixed_general(omb, bcVals, imin, imax);
#if (OMS_USER == TRUE)
			} else {
				NumMatrix<double,3> & bcVals = bcVals_User_ze[q];
				// Doing bcs:
				bc_Fixed_general(omb, bcVals, imin, imax);
#endif
			}

		}
	}


}

void gridFunc::bc_Fixed_general(NumMatrix<double,3> &omb,
                                NumMatrix<double,3> &bcVals,
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
    
	int nproc[3] = {1,1,1};
	int coords[3] = {0,0,0};
	h5out.AddGlobalAttr("procs",nproc,3);
	h5out.AddGlobalAttr("coords",coords,3);
	h5out.AddGlobalAttr("rank",gdata.rank);

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

	// In case of multifluid simulation do separate output for individual fluids
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

		int q_index = -1;
 	for (int q = qmin; q < qmax; ++q) {
    
		int qout = q;
		if(q >= n_omInt) {
			qout = q+N_ADD;
		}

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

			h5out.Write3DMatrix(dsetName, data, xmin, gdata.dx, group);

			// if possible add unit to field
			if(Problem.TrafoNorm != NULL) {

				string fieldName;
				if(qout < gdata.fluid.get_N_OMINT()) {
					fieldName = gdata.fluid.get_fieldName(qout);
					fieldName = gdata.om[qout].getName();
				} else {
					fieldName = dsetName;
				}

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
			h5out.Write3DMatrix(dsetName, data, xmin, gdata.dx, group, q_index);
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
	h5out.Write3DMatrix("carbuncle_flag", data, xmin, gdata.dx, group);

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

			h5out.Write3DMatrix(gdata.om_user[q].getName(), data,
			                    xmin, gdata.dx, group);

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

			h5out.Write3DMatrix(gdata.om_user[q].getName(), data,
					xmin, gdata.dx, group);


		}

	}
#endif

#if(CRONOS_OUTPUT_COMPATIBILITY != CRONOS_ON)
	h5out.CloseGroup(group);
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

	h5in.ReadGlobalAttr("time",gdata.time);
	h5in.ReadGlobalAttr("dt",gdata.dt);
	int rim(0); // Boundary points to be taken into account
	h5in.ReadGlobalAttr("rim",rim);
	h5in.ReadGlobalAttr("CflNumber",gdata.cfl);
	//   h5in.ReadGlobalAttr("xmin",*xb);
  
	int mxsm[DIM];

	// get group name:
	string groupName = gdata.fluid.get_Name();


	h5in.getSize(groupName, gdata.om[0].getName(), mxsm, DIM);

	for(int i=0; i<DIM; ++i) {
		mxsm[i] -= 2*rim+1;
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

		/*
		  }
		*/

		double facred = fac;
		int mxloc[3] = {mxsm[0], mxsm[1], mxsm[2]};
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
			imin[i] = -rim;
			imax[i] = gdata.mx[i]+rim;
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

#if(CRONOS_OUTPUT_COMPATIBILITY == CRONOS_ON)
		string groupName = "Data";
#else
		string groupName = gdata.fluid.get_Name();
#endif

	h5in.getSize(groupName, gdata.om[0].getName(), mxsm, DIM);

	for(int i=0; i<DIM; ++i) {
		mxsm[i] -= 2*rim+1;
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

	int shift[DIM];
	for(int i=0; i<DIM; ++i) {
		shift[i] = 0;
	}

	// In case of multifluid simulation separate output for individual fluids is done
#if(CRONOS_OUTPUT_COMPATIBILITY == CRONOS_ON)
	groupName = "Data";
#else
	groupName = gdata.fluid.get_Name();
#endif
	int qmin = 0;
	int qmax = n_omIntAll + N_SUBS - n_omIntUser;

	int qmin_user = 0;
	int qmax_user = n_omIntUser;

//	for (int q = 0; q < qmax + qmax_user; ++q) {
	for (int q = qmin; q < qmax + qmax_user; ++q) {

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

			gdata.om[qin].rename(h5in.GetDatasetName(q-qshift));

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

		double facred = fac;
		int mxloc[3] = {mxsm[0], mxsm[1], mxsm[2]};

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
			imin[i] = -rim;
			imax[i] = gdata.mx[i]+rim;
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
		/*
		  }
		*/
	}

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

#if(CRONOS_OUTPUT_COMPATIBILITY == CRONOS_ON)
	string groupName = "Data";
#else
	string groupName = gdata.fluid.get_Name();
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
	double eps = 1.e-42;
	for(int dir=0; dir<DIM; ++dir) {
		if(gdata.get_EdgeGridding()) {
			ratio(dir) = gdata.get_numCells_global(dir)/Nx_data[dir];
//			ratio(dir) = Nx_data[dir]/(gdata.get_numCells_global(dir));
			// Test if we have an integer number
			double d_ratio = gdata.get_numCells_global(dir)/(1.*Nx_data[dir]-eps);
			if(std::abs(ratio(dir)-d_ratio) > 1.e-12) {
				cerr << " Error in datain_collective - ratio of grids is not an integer " << endl;
				exit(3);
			}
//			cout << " ratio: " << ratio(dir) << " " << Nx_data[dir] << " " << gdata.get_numCells_global(dir) << endl;
		} else {
			ratio(dir) = gdata.global_mx[dir]/(mx_data[dir]);
//			ratio(dir) = mx_data[dir]/(gdata.global_mx[dir]);
			// Test if we have an integer number
			double d_ratio = gdata.global_mx[dir]/(1.*mx_data[dir]-eps);
//			double d_ratio = mx_data[dir]/(1.*gdata.global_mx[dir]-eps);
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

#if(CRONOS_OUTPUT_COMPATIBILITY == CRONOS_ON)
	groupName = "Data";
#else
	groupName = gdata.fluid.get_Name();
#endif
	int qmin = 0;
	int qmax = n_omIntAll + N_SUBS - n_omIntUser;

	int qmin_user = 0;
	int qmax_user = n_omIntUser;

	// Show all available fields in file

	for (int q = qmin; q < qmax + qmax_user; ++q) {

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

#if(CRONOS_OUTPUT_COMPATIBILITY == CRONOS_ON)
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
		if(q < n_omInt) {
			qin = q;
		} else if(q >= n_omInt && q < n_omInt+N_SUBS) {
			qin = q+N_ADD;
		} else {
			qin = q-(n_omInt+N_SUBS);
//			is_generic = false;
		}

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

#if(CRONOS_OUTPUT_COMPATIBILITY == CRONOS_ON)
			gdata.om[qin].rename(h5in.GetDatasetName(q-qshift));
#else
//			cerr << " Going for name of " << q-qshift << " " << qin << endl;
			gdata.om[qin].rename(h5in.GetDatasetName(groupName, q-qshift));
#endif
			dataName = gdata.om[qin].getName();
			h5in.Read3DMatrix(groupName, gdata.om[qin].getName(), Data);

#if (OMS_USER == TRUE)
		} else {
			gdata.om_user[qin].rename(h5in.GetDatasetName(groupName, q-qshift));
			dataName = gdata.om_user[qin].getName();
			h5in.Read3DMatrix(groupName, gdata.om_user[qin].getName(), Data);
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
		// Serial version
		h5in.Read3DMatrix(groupName, dsetName, data);

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
		// Serial version
		h5in.Read3DMatrix(groupName, dsetName, data);

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
