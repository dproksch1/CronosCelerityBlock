#ifndef CRONOS_GRID_H
#define CRONOS_GRID_H = 1

#include <matrix.H>
#include "constants.H"
#include "Hdf5File_cbase.H"

#include "gridParam.H"

#define REAL double

// Class holding the position of the grid-cells and cell-interfaces
// Here we use the matrix class in order to be able to use negative
// indices.
class Grid {
	void include_coordinateAxis();
	void init_nonlinearGrid();

public:
  /** Construct a grid according to the new gridding layout (no centers on the boundaries) */
	static Grid newInDomainGrid(const NumArray<double>& xb, const  NumArray<double>& xe, const NumArray<int>& Nx);
	static Grid newLegacyGrid(const NumArray<double>& xb, const  NumArray<double>& xe, const NumArray<int>& mx);
	static Grid newFromCat();

	Grid(bool constructFromCat = true);
	virtual ~Grid(){}
	// get position of cell centre - direction dependet versions
	REAL getCen_x(int) const, getCen_y(int) const, getCen_z(int) const;
	// get position of cell centre - direction independet version
	REAL getCen(int, int) const;
	// get position of cell eges - direction dependet versions
	// REAL getEdge_x(int), getEdge_y(int), getEdge_z(int);
	// get position of cell eges - direction independet version
	// REAL getEdge(int, int);
	// get position of left cell edge - direction dependent versions
	REAL getEdgL_x(int) const;
	REAL getEdgL_y(int) const;
	REAL getEdgL_z(int) const;
	// Same for right edge
	REAL getEdgR_x(int) const;
	REAL getEdgR_y(int) const;
	REAL getEdgR_z(int) const;
	// get position of cell eges - direction independet version
	REAL getEdgL(int, int) const;
	REAL getEdgR(int, int) const;
	REAL get_x(int, int) const;
	REAL get_y(int, int) const;
	REAL get_z(int, int) const;

	// Get positions globally:
	REAL get_x_global(int, int) const;
	REAL get_y_global(int, int) const;
	REAL get_z_global(int, int) const;
	REAL get_pos_global(int, int, int) const;

	REAL getCen_dx_global(int dir, int ii) const;

	bool get_EdgeGridding() const;
	int get_RankWidth(int dir) const;
	int get_RankShift(int dir) const;
	int get_numCells_global(int dir) const;
	
	REAL get_pos(int, int, int) const;
	REAL getCen_dx(int dir, int ii) const;
	REAL getCen_hx(int dir, int ii) const;
	REAL getCen_idx(int dir, int ii) const;
	int get_CellIndex(int dir, REAL position) const;
	REAL get_dx(int dir, int ii, int shift);
	REAL get_idx(int, int, int);
	REAL get_hx(int, int, int);
	REAL get_CellGeomTrafo(int, int, int);
	REAL get_AreaGeomTrafo(int, int, int, int);
	REAL get_AreaGeomTrafo_x(int, int, int);
	REAL get_AreaGeomTrafo_y(int, int, int);
	REAL get_AreaGeomTrafo_z(int, int, int);
	REAL get_CellVolume(int, int, int) const;
	REAL get_CellArea_x(int, int, int);
	REAL get_CellArea_y(int, int, int);
	REAL get_CellArea_z(int, int, int);
	REAL get_CellArea();
	REAL get_Volume();

	/** Maps from local (mpi) indices to global ones */
	int get_idx_global(int dir, int ii) const;

#if (NON_LINEAR_GRID == CRONOS_ON)
	// REAL get_GridFunction(ProblemType &Problem, int dir, REAL ratio);
	//	REAL get_GridFunction(int dir, REAL ratio);
	REAL get_GridFunction(int dir, REAL ratio, REAL (*func)(REAL));
	REAL get_ExpGrid(int dir, REAL ratio);
	void set_NonLinGridFunction();
	void set_UserGridFunction(int dir, REAL beg, REAL end, REAL len);
	void add_nonLinGridToHdf(Hdf5Stream &, int num_ghost);
#endif
	// void generate_Grid(ProblemType &Problem);
	//	void generate_Grid();
	//	void generate_Grid(REAL (*func)(int, REAL));
	void generate_Grid();
	int axis_treatment(int dir) const;
	int get_singularity_treatment(int dir, int top) const;
	int get_singularity_treatment(int bound) const;
	int singularity_treatment(int dir, int top) const;
#ifdef GEOM
	double h0(REAL, REAL, REAL) const;
	double h1(REAL, REAL, REAL) const;
	double h2(REAL, REAL, REAL) const;
	REAL h0(int, int, int, int, int, int) const;
	REAL h1(int, int, int, int, int, int) const;
	REAL h2(int, int, int, int, int, int) const;
	REAL getCen_h0(int, int, int) const;
	REAL getCen_h1(int, int, int) const;
	REAL getCen_h2(int, int, int) const;
	REAL getEdgL_h0(int, int, int) const;
	REAL getEdgL_h1(int, int, int) const;
	REAL getEdgL_h2(int, int, int) const;
#endif
	int mx[DIM] = { 0 }, global_mx[DIM] = { 0 }, coords[DIM] = { 0 };
	REAL xb[DIM] = { 0 }, xe[DIM] = { 0 }, dx[DIM] = { 0 }, Lx[DIM] = { 0 }, hx[DIM] = { 0 }, idx[DIM] = { 0 };
	REAL grid_ref[DIM] = { 0.0 };
	REAL global_xb[DIM] = { 0.0 }, global_xe[DIM] = { 0.0 };
	int rank = 0, rank_constz = 0, rank_constx = 0;
	int rank_constxy = 0;
	static constexpr short phi_dir = GEOM-1; // phi dir: 1 for cyl, 2 for sph

protected:
	int numCellsEff[DIM] = { 0 };
private:
	REAL get_x(REAL) const;
	REAL get_y(REAL) const;
	REAL get_z(REAL) const;
	// REAL get_x(int, REAL);

	static void ScaleGrid(NumArray<double>& xb, NumArray<double>& xe, const int &, const REAL &);
	static void CheckScale(NumArray<double>& xb, NumArray<double>& xe);
	void WriteGridInfo();
	bool WriteGridData();
	bool WriteGridStructure();
	int Shift_Rank(int, int) const;
	void MakeGridDir();

#if (NON_LINEAR_GRID == CRONOS_ON)
	NumMatrix<REAL,1> cellCentres[DIM], global_cellCentres[DIM];
	NumMatrix<REAL,1> cellEdges[DIM], global_cellEdges[DIM];
	NumMatrix<REAL,1> cellsDouble[DIM], global_cellsDouble[DIM];
	NumMatrix<REAL,1> dxCentres[DIM], global_dxCentres[DIM];
	NumMatrix<REAL,1> idxCentres[DIM];
	NumMatrix<REAL,1> dxDouble[DIM];
	NumMatrix<REAL,1> idxDouble[DIM];
	int numDouble[DIM], global_numDouble[DIM];
#endif
	int numCells[DIM] = { 0 }, global_numCells[DIM] = { 0 };
	int parametrisation[DIM] = { 0 };
	REAL x0[3] = { 0.0 };
	bool EdgeGridding = false;
	int include_coordAxisPhi = 0;
	int include_Singularity[2*DIM] = { 0 };
	char griddirname[255] = { '\0' };
	Interface_GridFunction *myNonlinGrid[DIM];

};

#endif