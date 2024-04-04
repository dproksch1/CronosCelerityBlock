#ifndef CRONOS_GRID_H
#define CRONOS_GRID_H = 1

#include <matrix.H>
#include "constants.H"
#include "Hdf5File_cbase.H"

#include "gridParam.H"

//#define REAL double

// Class holding the position of the grid-cells and cell-interfaces
// Here we use the matrix class in order to be able to use negative
// indices.
class Grid {
	void include_coordinateAxis();

public:
  /** Construct a grid according to the new gridding layout (no centers on the boundaries) */
	static Grid newInDomainGrid(const NumArray<double>& xb, const  NumArray<double>& xe, const NumArray<int>& Nx);
	static Grid newLegacyGrid(const NumArray<double>& xb, const  NumArray<double>& xe, const NumArray<int>& mx);
	static Grid newFromCat();

	Grid(bool constructFromCat = true);
	virtual ~Grid(){}
	// get position of cell centre - direction dependet versions
	double getCen_x(int) const, getCen_y(int) const, getCen_z(int) const;
	// get position of cell centre - direction independet version
	double getCen(int, int) const;
	// get position of cell eges - direction dependet versions
	// double getEdge_x(int), getEdge_y(int), getEdge_z(int);
	// get position of cell eges - direction independet version
	// double getEdge(int, int);
	// get position of left cell edge - direction dependent versions
	double getEdgL_x(int) const;
	double getEdgL_y(int) const;
	double getEdgL_z(int) const;
	// Same for right edge
	double getEdgR_x(int) const;
	double getEdgR_y(int) const;
	double getEdgR_z(int) const;
	// get position of cell eges - direction independet version
	double getEdgL(int, int) const;
	double getEdgR(int, int) const;
	double get_x(int, int) const;
	double get_y(int, int) const;
	double get_z(int, int) const;

	// Get positions globally:
	double get_x_global(int, int) const;
	double get_y_global(int, int) const;
	double get_z_global(int, int) const;
	double get_pos_global(int, int, int) const;

	double getCen_dx_global(int dir, int ii) const;

	bool get_EdgeGridding() const;
	int get_RankWidth(int dir) const;
	int get_RankShift(int dir) const;
	int get_numCells_global(int dir) const;
	
	double get_pos(int, int, int) const;
	double getCen_dx(int dir, int ii) const;
	double getCen_hx(int dir, int ii) const;
	double getCen_idx(int dir, int ii) const;
	int get_CellIndex(int dir, double position) const;
	double get_dx(int dir, int ii, int shift);
	double get_idx(int, int, int);
	double get_hx(int, int, int);
	double get_CellGeomTrafo(int, int, int);
	double get_AreaGeomTrafo(int, int, int, int);
	double get_AreaGeomTrafo_x(int, int, int);
	double get_AreaGeomTrafo_y(int, int, int);
	double get_AreaGeomTrafo_z(int, int, int);
	double get_CellVolume(int, int, int) const;
	double get_CellArea_x(int, int, int);
	double get_CellArea_y(int, int, int);
	double get_CellArea_z(int, int, int);
	double get_CellArea();
	double get_Volume();

	/** Maps from local (mpi) indices to global ones */
	int get_idx_global(int dir, int ii) const;

	void generate_Grid();
	int axis_treatment(int dir) const;
	int get_singularity_treatment(int dir, int top) const;
	int get_singularity_treatment(int bound) const;
	int singularity_treatment(int dir, int top) const;
#ifdef GEOM
	double h0(double, double, double) const;
	double h1(double, double, double) const;
	double h2(double, double, double) const;
	double h0(int, int, int, int, int, int) const;
	double h1(int, int, int, int, int, int) const;
	double h2(int, int, int, int, int, int) const;
	double getCen_h0(int, int, int) const;
	double getCen_h1(int, int, int) const;
	double getCen_h2(int, int, int) const;
	double getEdgL_h0(int, int, int) const;
	double getEdgL_h1(int, int, int) const;
	double getEdgL_h2(int, int, int) const;
#endif
	int mx[DIM] = { 0 }, global_mx[DIM] = { 0 }, coords[DIM] = { 0 };
	double xb[DIM] = { 0 }, xe[DIM] = { 0 }, dx[DIM] = { 0 }, Lx[DIM] = { 0 }, hx[DIM] = { 0 }, idx[DIM] = { 0 };
	double grid_ref[DIM] = { 0.0 };
	double global_xb[DIM] = { 0.0 }, global_xe[DIM] = { 0.0 };
	int rank = 0, rank_constz = 0, rank_constx = 0;
	int rank_constxy = 0;
	static constexpr short phi_dir = GEOM-1; // phi dir: 1 for cyl, 2 for sph

	// int ntasks = 0;
	// int MPI_periods[DIM] = { 0 };
	// int nproc[DIM] = { 0 };
	// int left = 0, right = 0;
	// int front = 0,back = 0;
	// int top = 0,bottom = 0;
	// MPI_Comm comm3d;

protected:
	int numCellsEff[DIM] = { 0 };
private:
	double get_x(double) const;
	double get_y(double) const;
	double get_z(double) const;
	// double get_x(int, double);

	static void ScaleGrid(NumArray<double>& xb, NumArray<double>& xe, const int &, const double &);
	static void CheckScale(NumArray<double>& xb, NumArray<double>& xe);
	void WriteGridInfo();
	bool WriteGridData();
	bool WriteGridStructure();
	int Shift_Rank(int, int) const;
	void MakeGridDir();
	// void GetMpi();

#if (NON_LINEAR_GRID == CRONOS_ON)
	NumMatrix<double,1> cellCentres[DIM], global_cellCentres[DIM];
	NumMatrix<double,1> cellEdges[DIM], global_cellEdges[DIM];
	NumMatrix<double,1> cellsDouble[DIM], global_cellsDouble[DIM];
	NumMatrix<double,1> dxCentres[DIM], global_dxCentres[DIM];
	NumMatrix<double,1> idxCentres[DIM];
	NumMatrix<double,1> dxDouble[DIM];
	NumMatrix<double,1> idxDouble[DIM];
	int numDouble[DIM], global_numDouble[DIM];
#endif
	int numCells[DIM] = { 0 }, global_numCells[DIM] = { 0 };
	int parametrisation[DIM] = { 0 };
	double x0[3] = { 0.0 };
	bool EdgeGridding = false;
	int include_coordAxisPhi = 0;
	int include_Singularity[2*DIM] = { 0 };
	char griddirname[255] = { '\0' };
	Interface_GridFunction *myNonlinGrid[DIM];

};

#endif