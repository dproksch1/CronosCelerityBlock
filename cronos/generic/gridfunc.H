#ifndef CRONOS_GRIDFUNC_H
#define CRONOS_GRIDFUNC_H

#include "data.H"
#include "problem.H"


class gridFunc {
public:
	gridFunc(Data &);
	void boundary(Data &, ProblemType &);
	void boundary(Data &, ProblemType &, NumMatrix<double,3> &omb, int);
	void boundary_old(Data &, ProblemType &, NumMatrix<double,3> &omb, int, int, int iFluid=0);

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

	/*! internal MPI-boundary conditions only
	 * BCs will NOT be applied to the real boundaries of the domain
	 * routine needed, e.g., for the magnetic vector potential
	 * @param gdata Grid and Array data class
	 * @param problem Class containing the user setup
	 * @param omb NumMatrix array holding the fields, for which bcs are to be computed
	 * */
	void boundary_MPI(Data &gdata, ProblemType &problem, NumMatrix<double,3> &omb, int, int, int iFluid=0);
	// void dataout(Data &, string &, ProblemType &, int, bool);
	void dataout(Data &, Hdf5Stream &, ProblemType &, int, bool, bool is_collective);
	void datain(Data &, Hdf5iStream &, string &, ProblemType &);
	void datain_old(Data &, Hdf5iStream &, string &, ProblemType &);
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
	void bc_SendRecv(Data &, NumMatrix<double,3> &, NumMatrix<double,3> &,
	                 int, int, int, bool, bool);
	void bc_SendRecv(Data &, NumArray<double> &, NumArray<double> &,
	                 int dir, int size_send, int size_recv, int, int, bool, bool);
#ifdef parallel
	void bc_ISendRecv(Data &, NumArray<double> &, NumArray<double> &, MPI_Request *,
	                 int dir, int size_send, int size_recv, int, int, bool, bool);
#else
	void bc_ISendRecv(Data &, NumArray<double> &, NumArray<double> &,
	                 int dir, int size_send, int size_recv, int, int, bool, bool);
#endif
	void bc_Periodic_old(Data &gdata, ProblemType &, NumMatrix<REAL,3> &omb,
	                 int dir, int q, int rim, bool internal_only=false);

	void bc_Periodic(Data &gdata, ProblemType &, NumMatrix<REAL,3> &omb,
	                 int dir, int q, int rim, bool internal_only=false);
	/*!
	 * Periodic boundary conditions for non-mpi version of the code, only
	 * */
#ifndef parallel
	void bc_Periodic_serial(Data &gdata, ProblemType &, NumMatrix<REAL,3> &omb,
			int dir, int q, int rim);
#endif

#ifdef parallel
	void bc_MPISend(Data &gdata, ProblemType &problem, NumMatrix<REAL,3> &omb,
			int dir, int q, int rim);
	void bc_MPIRecv(Data &gdata, ProblemType &, NumMatrix<REAL,3> &omb,
			int dir, int q, int rim);
#endif

#ifdef parallel
	/*!
	 * First routine for MPI boundary conditions. This part contains filling of
	 * send buffers and non-blocking send and receive routines.
	 * @param gdata Grid and Array data class
	 * @param problem Class containing the user setup
	 * @param omb NumMatrix array holding the fields, for which bcs are to be computed
	 * @param q index of omb field in global data array
	 * @param rim number of ghost cells to be taken into account
	 * */
	void bc_MPISend(Data &gdata, ProblemType &problem, NumMatrix<REAL,3> &omb,
			int q, int rim);

	/*!
	 * First routine for MPI boundary conditions.
	 * In this part MpiWait is executed until corresponding recv-buffers are filled
	 * Then recv-buffers are written back to computational fields
	 * @param gdata Grid and Array data class
	 * @param problem Class containing the user setup
	 * @param omb NumMatrix array holding the fields, for which bcs are to be computed
	 * @param q index of omb field in global data array
	 * @param rim number of ghost cells to be taken into account
	 * */
	void bc_MPIRecv(Data &gdata, ProblemType &, NumMatrix<REAL,3> &omb,
			int q, int rim);
#endif

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
            NumMatrix<REAL,3> &omb, int q, int rim);
	void bc_AxisSph(Data &gdata, ProblemType &Problem,
			NumMatrix<REAL,3> &omb, int top, int q, int rim);
	void bc_Axis(Data &gdata, ProblemType &Problem,
	             NumMatrix<REAL,3> &omb, int dir, int top, int q, int rim);
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
	void do_AxisValCorrectionCyl(Data &gdata, ProblemType &Problem);
	void do_AxisValCorrectionSph(Data &gdata, ProblemType &Problem, bool top);
	void bc_Extrapolate(Data &, ProblemType &, NumMatrix<double,3> &,
	                    int, int, int, int);
	void bc_Outflow(Data &, ProblemType &, NumMatrix<double,3> &,
	                int, int, int, int);
	void bc_Reflecting(Data &, ProblemType &, NumMatrix<double,3> &,
	                   int, int, int, int);
	void bc_Fixed(Data &, ProblemType &, NumMatrix<double,3> &,
	              int, int, int, int);
	void bc_Fixed_general(NumMatrix<REAL,3> &, NumMatrix<REAL,3> &,
	                      int imin[3], int imax[3]);
	void Prolongate(Data &, NumMatrix<double,3> &, const int &, int mxloc[], int);
	void ExchangePositions(Data &, int top=0);
	int get_powerTwo(int val);
#ifdef parallel
	void Assign(Data &, NumMatrix<double,3> &, int [DIM], int [DIM], int [DIM],
	            int [DIM],int );
	void get_AxisData(Data &gdata, 
			NumMatrix<REAL,3> &myBound,
			NumMatrix<REAL,3> &global_bcData, int top=0);
	void store_AxisData(Data &gdata, int lbound[], int ubound[],
	                    int rank_local,
	                    NumMatrix<REAL,3> &local_bcData,
	                    NumMatrix<REAL,3> &global_bcData);
	void ExchangeData(Data &);
#endif
	/**
	 * Adapt resolution of input data to version of cat file
	 * */
	void RefineData(Data &gdata, NumMatrix<double,3> &, int mx[DIM],
			NumArray<int> &factor, NumArray<int> &staggered, int rim);
	int bc_Type[6];
	NumMatrix<int,1> AxisPartnersPhi[2];  // array of mapping partners
	NumMatrix<REAL,1> PhiPos[2];
	NumMatrix<REAL,1> cosPos[3], sinPos[3];
	// NumMatrix<REAL,3> BoundVals[2*DIM];
	// NumMatrix<REAL,3> BoundVals_User[2*DIM];
	NumMatrix<REAL,3> bcVals_xb[N_OMINT], bcVals_xe[N_OMINT];
	NumMatrix<REAL,3> bcVals_yb[N_OMINT], bcVals_ye[N_OMINT];
	NumMatrix<REAL,3> bcVals_zb[N_OMINT], bcVals_ze[N_OMINT];

	NumMatrix<double,3> Recv_x, Recv_y, Recv_z, Recv_empty;
	NumMatrix<double,3> Send_x, Send_y, Send_z, Send_empty;

//	NumArray<REAL> RecvArr_x, RecvArr_y, RecvArr_z, RecvArr_empty;
//	NumArray<REAL> SendArr_x, SendArr_y, SendArr_z, SendArr_empty;
	NumArray<REAL> SendArr_buff, RecvArr_buff;
	NumArray<REAL> RecvArrSmall_x, RecvArrSmall_y, RecvArrSmall_z;
	NumArray<REAL> SendArrSmall_x, SendArrSmall_y, SendArrSmall_z;
	NumArray<REAL> RecvArr_xLR, RecvArr_yLR, RecvArr_zLR;
	NumArray<REAL> RecvArr_xRL, RecvArr_yRL, RecvArr_zRL;
	NumArray<REAL> SendArr_xLR, SendArr_yLR, SendArr_zLR;
	NumArray<REAL> SendArr_xRL, SendArr_yRL, SendArr_zRL;

	NumArray<REAL> SendArrRL, SendArrLR, RecvArrRL, RecvArrLR;

#ifdef parallel
	MPI_Request requests_xRL[2], requests_xLR[2];
	MPI_Request requests_yRL[2], requests_yLR[2];
	MPI_Request requests_zRL[2], requests_zLR[2];
	MPI_Request requests_RL[2], requests_LR[2];
#endif

#if (OMS_USER == TRUE)
	NumMatrix<REAL,3> bcVals_User_xb[N_OMINT], bcVals_User_xe[N_OMINT];
	NumMatrix<REAL,3> bcVals_User_yb[N_OMINT], bcVals_User_ye[N_OMINT];
	NumMatrix<REAL,3> bcVals_User_zb[N_OMINT], bcVals_User_ze[N_OMINT];
#endif
	int n_om, n_omInt, n_omIntAll;
	int n_omUser, n_omIntUser;
	int q_Bx, q_By, q_Bz;
};

#endif
