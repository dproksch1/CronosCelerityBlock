#ifndef CRONOS_GRIDFUNC_H
#define CRONOS_GRIDFUNC_H

#include "data.H"
#include "problem.H"
#include "Hdf5File_cbase_typemap.H"

namespace distr_io {

	template <typename T>
	std::pair<hid_t, hid_t> allocation_window_to_dataspace(const celerity::buffer_allocation_window<T, 3>& layout)
	{
		hsize_t file_size[3], file_start[3];
		hsize_t allocation_size[3], allocation_start[3];
		hsize_t count[3];
		for(int d = 0; d < 3; ++d) {
			file_size[d] = layout.get_buffer_range()[d];
			file_start[d] = layout.get_window_offset_in_buffer()[d];
			allocation_size[d] = layout.get_allocation_range()[d];
			allocation_start[d] = layout.get_window_offset_in_allocation()[d];
			count[d] = layout.get_window_range()[d];
		}

		auto file_space = H5Screate_simple(3, file_size, nullptr);
		H5Sselect_hyperslab(file_space, H5S_SELECT_SET, file_start, nullptr, count, nullptr);
		auto allocation_space = H5Screate_simple(3, allocation_size, nullptr);
		H5Sselect_hyperslab(allocation_space, H5S_SELECT_SET, allocation_start, nullptr, count, nullptr);

		return {file_space, allocation_space};
	}

	template<typename DataT>
	bool Write3DMatrix(Hdf5Stream &h5out, const std::string &ArrayName, celerity::buffer_allocation_window<DataT, 3> allocation_window, 
											const NumArray<int> &mx_global, int rim, const double *xb, const double *dx, hid_t my_group, int q_index,
											bool with_opendxinfo, bool is_float) 
	{

		h5out.AddDatasetName(ArrayName, my_group);
		h5out.AddDatasetName(ArrayName);
		h5out.increase_num();

		int dim = 3;
		hid_t datatype;
		if (is_float) {
		  	datatype = get_hdf5_data_type<float>();
		} else {
			datatype = get_hdf5_data_type<double>();
		}
		H5Tset_order(datatype, H5T_ORDER_LE);

		hsize_t DimsData[dim];
		DimsData[0] = mx_global[2] + 1 + 2 * rim;
		DimsData[1] = mx_global[1] + 1 + 2 * rim;
		DimsData[2] = mx_global[0] + 1 + 2 * rim;

		hid_t dataspace = H5Screate_simple(dim, DimsData, NULL);

		// Supplying additional attributes for opendx input

		hid_t datatypefloat = get_hdf5_data_type<double>();
		H5Tset_order(datatypefloat, H5T_ORDER_LE);

		// Create dataspace for attribute
		hsize_t DimsAttr = 3;
		hid_t attrspace = H5Screate_simple(1, &DimsAttr, NULL);

		double Origin[3];
		double Delta[3];
		if (with_opendxinfo) {
			for (int q = 0; q < 3; ++q) {
			Origin[q] = float(xb[q]);
			Delta[q] = float(dx[q]);
			}
		}

		// Create dataset
		hid_t dataset_id = H5Dcreate2(my_group, ArrayName.c_str(), datatype, dataspace, H5P_DEFAULT,
										H5P_DEFAULT, H5P_DEFAULT);

		if (with_opendxinfo) {
			// Create attributes
			hid_t origin =
				H5Acreate2(dataset_id, "origin", datatypefloat, attrspace, H5P_DEFAULT, H5P_DEFAULT);
			hid_t delta =
				H5Acreate2(dataset_id, "delta", datatypefloat, attrspace, H5P_DEFAULT, H5P_DEFAULT);

			// hid_t test_attr =
			// 	H5Acreate2(dataset_id, "test_attr", datatypefloat, attrspace, H5P_DEFAULT, H5P_DEFAULT);

			// Write attributes
			H5Awrite(origin, datatypefloat, &Origin);
			H5Awrite(delta, datatypefloat, &Delta);

			// Close attributes
			H5Aclose(origin);
			H5Aclose(delta);
			H5Sclose(attrspace);
		}

		// if (q_index > 0) {
		// 	// Create dataspace
		// 	hid_t info_id = H5Screate(H5S_SCALAR);
		// 	// Create Attribute
		// 	hid_t info =
		// 		H5Acreate2(dataset_id, "q_index", H5T_NATIVE_INT, info_id, H5P_DEFAULT, H5P_DEFAULT);
		// 	H5Awrite(info, H5T_NATIVE_INT, &q_index);
		// 	H5Aclose(info);
		// 	H5Sclose(info_id);
		// }

		auto plist = H5Pcreate(H5P_DATASET_XFER);
			H5Pset_dxpl_mpio(plist, H5FD_MPIO_COLLECTIVE);

		// auto allocation_window = data_acc.get_allocation_window(part);
		auto [file_space, allocation_space] = allocation_window_to_dataspace(allocation_window);
		H5Dwrite(dataset_id, datatype, allocation_space, file_space, plist, allocation_window.get_allocation());

		H5Sclose(dataspace);

		// Close the dataset:
		H5Dclose(dataset_id);
		H5Sclose(allocation_space);
		H5Sclose(file_space);
		H5Pclose(plist);


		return true;
	}

	template<typename DataT>
	void dataout(Data &gdata,  Hdf5Stream &h5out,
						ProblemType & Problem, int numout, bool isfloat, celerity::buffer_allocation_window<DataT, 3> om_rho_acc,
					   celerity::buffer_allocation_window<DataT, 3> om_sx_acc, celerity::buffer_allocation_window<DataT, 3> om_sy_acc, 
					   celerity::buffer_allocation_window<DataT, 3> om_sz_acc, celerity::buffer_allocation_window<DataT, 3> om_Eges_acc)
	{
		Problem.WriteToH5(h5out);

		int rim(0);
		if(isfloat) {
			rim = BOUT_FLT;
		} else {
			rim = BOUT_DBL;
		}

		NumArray<int> mx_global(3);
		NumArray<float> xmin_global(3), delx(3);
		for (int dir = 0; dir < 3; dir++) {
			mx_global(dir) = gdata.global_mx[dir];
			xmin_global(dir) = gdata.get_pos_global(dir, rim, 0);
			delx(dir) = gdata.dx[dir];
		}


		double xmin[3];
		xmin[0] = gdata.getCen_x(-rim);
		xmin[1] = gdata.getCen_y(-rim);
		xmin[2] = gdata.getCen_z(-rim);

		int itime = static_cast<int>(gdata.time/gdata.dt);

		int nproc[3] = {1,1,1};
		int coords[3] = {0,0,0};
		h5out.AddGlobalAttr("procs",nproc,3);
		h5out.AddGlobalAttr("coords",coords,3);
		h5out.AddGlobalAttr("rank",gdata.rank);
		h5out.AddGlobalAttr("fluid_type", "hydro");

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

	#if(CRONOS_OUTPUT_COMPATIBILITY == CRONOS_ON)
		hid_t group = h5out.get_defaultGroup();
	#else
		string groupName = gdata.fluid.get_Name();
		hid_t group = h5out.AddGroup(groupName);
	#endif

		// if (isfloat) {

		// } else {
			Write3DMatrix(h5out, gdata.om[0].getName(), om_rho_acc, mx_global, rim,
						xmin, gdata.dx, group, 0, true, isfloat);
			Write3DMatrix(h5out, gdata.om[1].getName(), om_sz_acc, mx_global, rim,
						xmin, gdata.dx, group, 1, true, isfloat);
			Write3DMatrix(h5out, gdata.om[2].getName(), om_sy_acc, mx_global, rim,
						xmin, gdata.dx, group, 2, true, isfloat);
			Write3DMatrix(h5out, gdata.om[3].getName(), om_sx_acc, mx_global, rim,
						xmin, gdata.dx, group, 3, true, isfloat);
			Write3DMatrix(h5out, gdata.om[4].getName(), om_Eges_acc, mx_global, rim,
						xmin, gdata.dx, group, 4, true, isfloat);
		// }

	#if(CRONOS_OUTPUT_COMPATIBILITY != CRONOS_ON)
		h5out.CloseGroup(group);
	#endif

	}
	
}




class gridFunc {
public:
	gridFunc(Data &);
	void boundary(Data &, ProblemType &);
	void boundary(Data &, ProblemType &, NumMatrix<double,3> &omb, int);
	// void boundary_old(Data &, ProblemType &, NumMatrix<double,3> &omb, int, int, int iFluid=0);

	/*!
	 * Optimised boundary condition routine for MPI parallel case - trying to avoid stalls by
	 * MPI communication
	 * @param gdata Grid and Array data class
	 * @param problem Class containing the user setup
	 * @param omb NumMatrix array holding the fields, for which bcs are to be computed
	 * @param rim number of ghost cells to be taken into account
	 * @param q index of omb field in global data array
	 * @param iFluid index of fluid considered (only applicable in multifluid simulations)
	 * */
	void boundary(Data &gdata, ProblemType &problem, NumMatrix<double,3> &omb,
			int rim, int q, int iFluid=0);

	void boundary(Queue &, Data &gdata, ProblemType &problem, NumMatrix<double,3> &omb,
			int rim, int q, int iFluid=0);

	// void dataout(Data &, string &, ProblemType &, int, bool);

	void dataout(Data &, Hdf5Stream &, ProblemType &, int, bool, bool is_collective);
	// void dataout2(Queue &, Data &, Hdf5Stream &, ProblemType &, int, bool, bool is_collective);
	void datain(Data &, Hdf5iStream &, string &, ProblemType &);
	// void datain_old(Data &, Hdf5iStream &, string &, ProblemType &);
	/*!
	 * Function to load restart data from runs, which wrote all data into a single file
	 * the maximum refinement factor is returned
	 * */
	double datain_collective(Data &, Hdf5iStream &, string &, ProblemType &);
#if(FLUID_TYPE != CRONOS_MULTIFLUID)
	void load_flt(Data &, Hdf5iStream &, string &, ProblemType &);
#endif
	void prep_boundaries(Data &, ProblemType &);
	int get_bc_Type(int);
private:
	// void bc_SendRecv(Data &, NumMatrix<double,3> &, NumMatrix<double,3> &,
	//                  int, int, int, bool, bool);
	// void bc_SendRecv(Data &, NumArray<double> &, NumArray<double> &,
	//                  int dir, int size_send, int size_recv, int, int, bool, bool);
// #ifdef parallel
// 	void bc_ISendRecv(Data &, NumArray<double> &, NumArray<double> &, MPI_Request *,
// 	                 int dir, int size_send, int size_recv, int, int, bool, bool);
// #else
// 	void bc_ISendRecv(Data &, NumArray<double> &, NumArray<double> &,
// 	                 int dir, int size_send, int size_recv, int, int, bool, bool);
// #endif
	// void bc_Periodic_old(Data &gdata, ProblemType &, NumMatrix<double,3> &omb,
	//                  int dir, int q, int rim, bool internal_only=false);

	// void bc_Periodic(Data &gdata, ProblemType &, NumMatrix<double,3> &omb,
	//                  int dir, int q, int rim, bool internal_only=false);
	/*!
	 * Periodic boundary conditions for non-mpi version of the code, only
	 * */
	void bc_Periodic_serial(Data &gdata, ProblemType &, NumMatrix<double,3> &omb,
			int dir, int q, int rim);
		
	void bc_Periodic(Queue &, Data &gdata, ProblemType &, NumMatrix<double,3> &omb,
			int dir, int q, int rim);

	/*!
	 * Determine size of send buffer in send direction
	 * @param gdata Grid and Array data class
	 * @param range_min begin of buffer
	 * @param end of buffer
	 * @param dir Cartesian direction for boundary communication
	 * @param leftToRight data transfer direction (1 for left to right, 0 else)
	 * @param q index of omb field in global data array
	 * @param rim number of ghost cells to be taken into account
	 * */
	void get_bcBuffRange(Data &gdata, int &range_min, int & range_max,
			int dir, int leftToRight, int q, int rim );

	void bc_AxisCyl(Data &gdata, ProblemType &Problem,
            NumMatrix<double,3> &omb, int q, int rim);
	void bc_AxisSph(Data &gdata, ProblemType &Problem,
			NumMatrix<double,3> &omb, int top, int q, int rim);
	void bc_Axis(Data &gdata, ProblemType &Problem,
	             NumMatrix<double,3> &omb, int dir, int top, int q, int rim);
	/*!
	 * Compute location of related ghost cells at coordinate axes
	 * @param top 0 for bottom and 1 for top end of perp axis
	 * */
	void compute_AxisPartners(Data &gdata, int top=0);
	/*!
	 * Get index of related ghost cell at coordinate axis
	 * -> routine takes care about correct mapping for MPI-parallel case
	 * @param top 0 for bottom and 1 for top end of perp axis
	 * @param iz index of coordinate surrounding axis
	 * */
	int get_AxisPartners(Data &gdata, int top, int iz);
	void do_AxisValCorrectionSph(Data &gdata, ProblemType &Problem, bool top);
	void bc_Extrapolate(Data &, ProblemType &, NumMatrix<double,3> &,
	                    int, int, int, int);
	void bc_Outflow(Data &, ProblemType &, NumMatrix<double,3> &,
	                int, int, int, int);
	void bc_Reflecting(Data &, ProblemType &, NumMatrix<double,3> &,
	                   int, int, int, int);
	void bc_Fixed(Data &, ProblemType &, NumMatrix<double,3> &,
	              int, int, int, int);
	void bc_Fixed_general(NumMatrix<double,3> &, NumMatrix<double,3> &,
	                      int imin[3], int imax[3]);
	void Prolongate(Data &, NumMatrix<double,3> &, const int &, int mxloc[], int);
	int get_powerTwo(int val);

	void bc_Outflow(Queue &queue, Data &, ProblemType &, NumMatrix<double,3> &,
	                int, int, int, int);
	void bc_Extrapolate(Queue &queue, Data &, ProblemType &, NumMatrix<double,3> &,
	                    int, int, int, int);
	void bc_Reflecting(Queue &queue, Data &, ProblemType &, NumMatrix<double,3> &,
	                   int, int, int, int);

	/**
	 * Adapt resolution of input data to version of cat file
	 * */
	void RefineData(Data &gdata, NumMatrix<double,3> &, int mx[DIM],
			NumArray<int> &factor, NumArray<int> &staggered, int rim);
	int bc_Type[6];
	NumMatrix<int,1> AxisPartnersPhi[2];  // array of mapping partners
	NumMatrix<double,1> PhiPos[2];
	NumMatrix<double,1> cosPos[3], sinPos[3];
	NumMatrix<double,3> bcVals_xb[N_OMINT], bcVals_xe[N_OMINT];
	NumMatrix<double,3> bcVals_yb[N_OMINT], bcVals_ye[N_OMINT];
	NumMatrix<double,3> bcVals_zb[N_OMINT], bcVals_ze[N_OMINT];

	NumMatrix<double,3> Recv_x, Recv_y, Recv_z, Recv_empty;
	NumMatrix<double,3> Send_x, Send_y, Send_z, Send_empty;

	NumArray<double> SendArr_buff, RecvArr_buff;
	NumArray<double> RecvArrSmall_x, RecvArrSmall_y, RecvArrSmall_z;
	NumArray<double> SendArrSmall_x, SendArrSmall_y, SendArrSmall_z;
	NumArray<double> RecvArr_xLR, RecvArr_yLR, RecvArr_zLR;
	NumArray<double> RecvArr_xRL, RecvArr_yRL, RecvArr_zRL;
	NumArray<double> SendArr_xLR, SendArr_yLR, SendArr_zLR;
	NumArray<double> SendArr_xRL, SendArr_yRL, SendArr_zRL;

	NumArray<double> SendArrRL, SendArrLR, RecvArrRL, RecvArrLR;

#ifdef parallel
	MPI_Request requests_xRL[2], requests_xLR[2];
	MPI_Request requests_yRL[2], requests_yLR[2];
	MPI_Request requests_zRL[2], requests_zLR[2];
	MPI_Request requests_RL[2], requests_LR[2];
#endif

#if (OMS_USER == TRUE)
	NumMatrix<double,3> bcVals_User_xb[N_OMINT], bcVals_User_xe[N_OMINT];
	NumMatrix<double,3> bcVals_User_yb[N_OMINT], bcVals_User_ye[N_OMINT];
	NumMatrix<double,3> bcVals_User_zb[N_OMINT], bcVals_User_ze[N_OMINT];
#endif
	int n_om, n_omInt, n_omIntAll;
	int n_omUser, n_omIntUser;
	int q_Bx, q_By, q_Bz;
};

#endif
