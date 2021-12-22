#include "grid.H"

#include <stdio.h>
#include <string.h>
#include <sys/stat.h>
#ifdef _MSC_VER
#include "timewrapper.H"
#else
#include <utime.h>
#endif
#include <stdlib.h>
#include <iostream>
#include <fstream>
#include <filesystem>

// Included for eclipse
#include <assert.h>
#include <math.h>
using namespace std;

// Constructor for grid class - set size of matrix holding the grid
Grid::Grid(bool constructFromCat) {
  if (constructFromCat) {
    *this = newFromCat();
  }
}

Grid Grid::newInDomainGrid(const NumArray<double>& _xb, const  NumArray<double>& _xe, const NumArray<int>& _Nx) {
  Grid grid(false);

  // New gridding scheme - Nx = number of in-domain cells
  grid.EdgeGridding = true;

  for (int dir = 0; dir < DIM; ++dir) {
    grid.xb[dir] = _xb[dir];
    grid.xe[dir] = _xe[dir];
    grid.numCells[dir] = _Nx[dir];

    grid.numCellsEff[dir] = grid.numCells[dir];
    grid.global_numCells[dir] = grid.numCells[dir];
    grid.mx[dir] = grid.numCells[dir] - 1;
    grid.global_mx[dir] = grid.mx[dir];
  }

#if (NON_LINEAR_GRID == CRONOS_OFF)
  // For a linear grid we need dx and similar stuff:
  for (int dir = 0; dir < DIM; ++dir) {
    grid.Lx[dir] = grid.xe[dir] - grid.xb[dir];
    grid.dx[dir] = grid.Lx[dir] / grid.numCells[dir];
    grid.idx[dir] = 1. / grid.dx[dir];
    grid.hx[dir] = 0.5 / grid.dx[dir];
    // Set centre of first cell (at i=0)
    grid.x0[dir] = grid.xb[dir] + 0.5 * grid.dx[dir];
  }
#endif
  for (int dir = 0; dir < DIM; ++dir) {
    grid.global_xb[dir] = grid.xb[dir];
    grid.global_xe[dir] = grid.xe[dir];
  }

  grid.rank = 0;
  for (int i = 0; i < DIM; ++i) {
    grid.coords[i] = 0;
  }

  return grid;
}

Grid Grid::newLegacyGrid(const NumArray<double>& _xb, const  NumArray<double>& _xe, const NumArray<int>& _mx) {
  Grid grid(false);

  // Old gridding scheme - mx = number of lass in-domain cell
  grid.EdgeGridding = false;

  for (int dir = 0; dir < DIM; ++dir) {
    grid.xb[dir] = _xb[dir];
    grid.xe[dir] = _xe[dir];
    grid.mx[dir] = _mx[dir];

    grid.numCellsEff[dir] = grid.mx[dir];
    grid.global_mx[dir] = grid.mx[dir];
    grid.numCells[dir] = grid.mx[dir] + 1;
    grid.global_numCells[dir] = grid.numCells[dir];
  }

#if (NON_LINEAR_GRID == CRONOS_OFF)
  // For a linear grid we need dx and similar stuff:
  for (int dir = 0; dir < DIM; ++dir) {
    grid.Lx[dir] = grid.xe[dir] - grid.xb[dir];
    grid.dx[dir] = grid.Lx[dir] / grid.mx[dir];
    grid.idx[dir] = 1. / grid.dx[dir];
    grid.hx[dir] = 0.5 / grid.dx[dir];
    // Set centre of first cell (at i=0)
    grid.x0[dir] = grid.xb[dir];
  }
#endif
  for (int dir = 0; dir < DIM; ++dir) {
    grid.global_xb[dir] = grid.xb[dir];
    grid.global_xe[dir] = grid.xe[dir];
  }

  grid.rank = 0;
  for (int i = 0; i < DIM; ++i) {
    grid.coords[i] = 0;
  }

  return grid;
}

Grid Grid::newFromCat() {
	NumArray<double> xb = NumArray<double>::set(value("xb"), value("yb"), value("zb"));
	NumArray<double> xe = NumArray<double>::set(value("xe"), value("ye"), value("ze"));

  bool newGridding = value_exists("Nx") && value_exists("Ny") && value_exists("Nz");

  Grid grid(false);

	// Scaling grid by some constant factor
  bool Scale = static_cast<bool>(value("ScaleGrid"));

  if (Scale) {
    CheckScale(xb, xe);
  }

  if (newGridding) {
    NumArray<int> Nx = NumArray<int>::set(value("Nx"), value("Ny"), value("Nz"));
    grid = Grid::newInDomainGrid(xb, xe, Nx);
  } else {
    NumArray<int> mx = NumArray<int>::set(value("mx"), value("my"), value("mz"));
    grid = Grid::newLegacyGrid(xb, xe, mx);
  }

  grid.init_nonlinearGrid();

  if (value_exists("Include_CoordinateAxis")) {
    grid.include_coordinateAxis();
  }

  grid.WriteGridInfo();

  // Writing the grid structure to be used by paraview for visualisation
  if (value_exists((char *)"WriteGrid")) {
    // WriteGridData();
    grid.WriteGridStructure();
  }

  return grid;
}

void Grid::include_coordinateAxis() {
  include_Singularity[0] = 0;  // lower x
  include_Singularity[1] = 0;  // upper x
  include_Singularity[2] = 0;  // lower y
  include_Singularity[3] = 0;  // upper y
  include_Singularity[4] = 0;  // lower z
  include_Singularity[5] = 0;  // upper z

#if (GEOM == CYLINDRICAL)

  // dx at lower global boundary
  double dxl_glob = getCen_dx_global(0, 0);
  // Check if grid on axis:
  if ((get_EdgeGridding() && (global_xb[0] < 0.1 * dxl_glob)) ||
      (!get_EdgeGridding() && (global_xb[0] < 0.6 * dxl_glob))) {
    //			include_coordAxisPhi = 2;
    include_Singularity[0] = 2;

    // Check if grid on axis for current rank:
    // dx at lower boundary
    double dxl_loc = getCen_dx(0, 0);
    if ((get_EdgeGridding() && (xb[0] < 0.1 * dxl_loc)) || (!get_EdgeGridding() && (xb[0] < 0.6 * dxl_loc))) {
      //				include_coordAxisPhi = 1;
      include_Singularity[0] = 1;
    }
  } else {
#ifndef parallel
    cerr << " ERROR: Grid does not end on-axis. " << endl;
    cerr << "   switch off /Include_CoordinateAxis/ for code to be allowed to run " << endl;
    exit(3);
#endif
  }
#elif (GEOM == SPHERICAL)

  // Check if theta reaches the boundary
  // cout << " gridding " << get_EdgeGridding() << " ";
  // cout << global_xb[1] << " " << dx[1] << endl;

  // cout << getCen_dx_global(1,0) << endl;

  // exit(3);

  if (get_EdgeGridding()) {
    // dy at lower global y-boundary
    double dyl_glob = getCen_dx_global(1, 0);

    // Check if axis is reached on some rank
    if ((global_xb[1] < 0.1 * dyl_glob) || (global_xb[1] < 0.6 * dyl_glob)) {
      include_Singularity[2] = 2;

      // dy at lower local y-boundary
      double dyl_loc = getCen_dx(1, 0);

      // Check if boundary included in current rank
      if ((xb[1] < 0.1 * dyl_loc) || (xb[1] < 0.6 * dyl_loc)) {
        include_Singularity[2] = 1;
      }
    }

    // same at upper boundary

    // dy at upper global y-boundary
    double dyr_glob = getCen_dx_global(1, global_mx[1]);

    if ((global_xe[1] > pi - 0.1 * dyr_glob) || (global_xe[1] > pi - 0.6 * dyr_glob)) {
      include_Singularity[3] = 2;

      // dy at upper global y-boundary
      double dyr_loc = getCen_dx_global(1, mx[1]);

      // Check if boundary included in current rank
      if ((xe[1] > pi - 0.1 * dyr_loc) || (xe[1] > pi - 0.6 * dyr_loc)) {
        include_Singularity[3] = 1;
      }
    }
  }
#endif
}

void Grid::init_nonlinearGrid() {
#if (NON_LINEAR_GRID == CRONOS_ON)
  // Set size of 1D grids
  for (int dir = 0; dir < DIM; ++dir) {
    // cellCentres[q].resize(Index::set(-B, -B, -B),
    //                       Index::set(mx[0]+B, mx[1]+B, mx[2]+B));
    // cellCentres[q].resize(Index::set(-B, -B, -B),
    //                       Index::set(numCells[0]+B, numCells[1]+B,
    //                                  numCells[2]+B));
    cellCentres[dir].resize(Index::set(-B), Index::set(numCells[dir] + B));
    global_cellCentres[dir].resize(Index::set(-B), Index::set(global_numCells[dir] + B));
    // One more edge than centers...
    // cellEdges[q].resize(Index::set(-B, -B, -B),
    //                     Index::set(mx[0]+B+1, mx[1]+B+1, mx[2]+B+1));
    // cellEdges[q].resize(Index::set(-B, -B, -B),
    //                     Index::set(numCells[0]+B+1, numCells[1]+B+1,
    //                                numCells[2]+B+1));
    cellEdges[dir].resize(Index::set(-B), Index::set(numCells[dir] + B + 1));
    global_cellEdges[dir].resize(Index::set(-B), Index::set(global_numCells[dir] + B + 1));

    // List of dx values at cell centres:
    // dxCentres[q].resize(Index::set(-B, -B, -B),
    //                     Index::set(mx[0]+B, mx[1]+B, mx[2]+B));
    // dxCentres[q].resize(Index::set(-B, -B, -B),
    //                     Index::set(numCells[0]+B, numCells[1]+B,
    //                                numCells[2]+B));
    dxCentres[dir].resize(Index::set(-B), Index::set(numCells[dir] + B));
    global_dxCentres[dir].resize(Index::set(-B), Index::set(global_numCells[dir] + B));
    // List of inverse dx values at cell centres:
    // idxCentres[q].resize(Index::set(-B, -B, -B),
    //                      Index::set(mx[0]+B, mx[1]+B, mx[2]+B));
    // idxCentres[q].resize(Index::set(-B, -B, -B),
    //                      Index::set(numCells[0]+B, numCells[1]+B,
    //                                 numCells[2]+B));
    idxCentres[dir].resize(Index::set(-B), Index::set(numCells[dir] + B));

    // cellsDouble[q].resize(Index::set(-2*B-1, -2*B-1, -2*B-1),
    //                       Index::set(2*(mx[0]+B)+1,
    //                                  2*(mx[1]+B)+1,
    //                                  2*(mx[2]+B)+1));
    // cellsDouble[q].resize(Index::set(-2*B-1, -2*B-1, -2*B-1),
    //                       Index::set(2*(numCells[0]+B)+1,
    //                                  2*(numCells[1]+B)+1,
    //                                  2*(numCells[2]+B)+1));
    cellsDouble[dir].resize(Index::set(-2 * B - 1), Index::set(2 * (numCells[dir] + B) + 1));
    global_cellsDouble[dir].resize(Index::set(-2 * B - 1), Index::set(2 * (global_numCells[dir] + B) + 1));
    numDouble[0] = 2 * (numCells[0] + B) + 1;
    numDouble[1] = 2 * (numCells[1] + B) + 1;
    numDouble[2] = 2 * (numCells[2] + B) + 1;

    global_numDouble[0] = 2 * (global_numCells[0] + B) + 1;
    global_numDouble[1] = 2 * (global_numCells[1] + B) + 1;
    global_numDouble[2] = 2 * (global_numCells[2] + B) + 1;

    // List of dx values at all relevant positions
    // dxDouble[q].resize(Index::set(-2*B-1, -2*B-1, -2*B-1),
    //                    Index::set(2*(mx[0]+B)+1,
    //                               2*(mx[1]+B)+1,
    //                               2*(mx[2]+B)+1));
    // dxDouble[q].resize(Index::set(-2*B-1, -2*B-1, -2*B-1),
    //                    Index::set(2*(numCells[0]+B)+1,
    //                               2*(numCells[1]+B)+1,
    //                               2*(numCells[2]+B)+1));

    dxDouble[dir].resize(Index::set(-2 * B - 1), Index::set(2 * (numCells[dir] + B) + 1));

    // List of inverse dx values at all relevant positions
    // idxDouble[q].resize(Index::set(-2*B-1, -2*B-1, -2*B-1),
    //                     Index::set(2*(mx[0]+B)+1,
    //                                2*(mx[1]+B)+1,
    //                                2*(mx[2]+B)+1));
    // idxDouble[q].resize(Index::set(-2*B-1, -2*B-1, -2*B-1),
    //                     Index::set(2*(numCells[0]+B)+1,
    //                                2*(numCells[1]+B)+1,
    //                                2*(numCells[2]+B)+1));
    idxDouble[dir].resize(Index::set(-2 * B - 1), Index::set(2 * (numCells[dir] + B) + 1));
  }

  for (int d = 0; d < DIM; ++d) {
    Lx[d] = xe[d] - xb[d];
  }

  // Pick a non-linear grid function
  set_NonLinGridFunction();

  // Make the grid:
  generate_Grid();

#endif
}

#if (NON_LINEAR_GRID == CRONOS_ON)
void Grid::set_NonLinGridFunction() {

	// Grid parametrisations:
	// (0) - linear
	// (1) - exponential
	// (2) - logarithmic
	// (4) - user
	// (5) - percent-wise increase of cell size
	// (6) - slightly different version of (5)

	parametrisation[0] = static_cast<int>(value((char*)"grid_x"));
	parametrisation[1] = static_cast<int>(value((char*)"grid_y"));
	parametrisation[2] = static_cast<int>(value((char*)"grid_z"));
	



	if(rank == 0) {
		cout << " Building the grid " << endl;
	}

	for(int dir=0; dir<DIM; ++dir) {

		double ref_par = 0.;
		double cell_change = 0.;
		double ref_cellSize = 1.;
		if(parametrisation[dir]==3) {
			if(dir==0) {
				ref_par = value((char*)"gridRef_x");
			} else if (dir==1) {
				ref_par = value((char*)"gridRef_y");
			} else {
				ref_par = value((char*)"gridRef_z");
			}
		} else if(parametrisation[dir]==5 || parametrisation[dir]==6) {
			if(dir==0) {
				cell_change = value((char*)"cellChange_x");
			} else if (dir==1) {
				cell_change = value((char*)"cellChange_y");
			} else {
				cell_change = value((char*)"cellChange_z");
			}

			if(parametrisation[dir]==6) {
				if(dir==0) {
					ref_cellSize = value((char*)"refSize_x");
				} else if (dir==1) {
					ref_cellSize = value((char*)"refSize_y");
				} else {
					ref_cellSize = value((char*)"refSize_z");
				}
			}

		}

		switch (parametrisation[dir]) {
		case 0 :
			myNonlinGrid[dir] = new GridFunction_Lin();
			break;
		case 1 :
			myNonlinGrid[dir] = new GridFunction_Exp(global_xb[dir], 
			                                         global_xe[dir],
			                                         Lx[dir]); 
			break;
		case 2 :
			myNonlinGrid[dir] = new GridFunction_Log(global_xb[dir],global_xe[dir]);
			break;
		case 3:
			// Load necessary variables
			// Get reference point
			myNonlinGrid[dir] = new GridFunction_Ratio(ref_par,global_xb[dir],global_xe[dir]);
			break;
		case 5:
			// Load necessary variables
			// Get reference point
			myNonlinGrid[dir] = new GridFunction_Percent(global_xb[dir],global_xe[dir],
					cell_change, global_numCells[dir]);
			break;
		case 6:
			// Load necessary variables
			// Get reference point
			myNonlinGrid[dir] = new GridFunction_PercentVar(global_xb[dir], ref_cellSize,
					cell_change, global_numCells[dir]);
			break;
		case 4:
			set_UserGridFunction(dir, global_xb[dir], global_xe[dir], Lx[dir]);
			break;
		}

	}
}
#endif

int Grid::get_RankWidth(int dir) const {
	return numCellsEff[dir];
}

int Grid::get_numCells_global(int dir) const {
	assert(dir>=0 && dir<DIM);
	return global_numCells[dir];
}


int Grid::get_RankShift(int dir) const {
	return coords[dir]*numCellsEff[dir];
}

bool Grid::get_EdgeGridding() const {
	return EdgeGridding;
}


int Grid::get_idx_global(int dir, int ii) const {
	return get_RankShift(dir) + ii;
}



#if (NON_LINEAR_GRID == CRONOS_ON)

int Grid::get_CellIndex(int dir, REAL position) const {
	int index(-99);
	for(int ii=0; ii<=mx[dir]; ++ii) {
		if(position >= getEdgL(dir, ii) && position <= getEdgR(dir, ii)) {
			index = ii;
		}
	}
	return index;
}

REAL Grid::getCen_x(int ii) const {
	assert(ii <= numCells[0]+B-1);
	return cellCentres[0](ii);
}

double Grid::get_x(int ii, int shift_i) const {
	int iDouble = 2*ii + shift_i;
	assert(iDouble <= numDouble[0]);
	return cellsDouble[0](iDouble);
}

// Global version:
double Grid::get_x_global(int ii, int shift_i) const {
	int iDouble = 2*ii + shift_i;
	assert(iDouble <= global_numDouble[0]);
	return global_cellsDouble[0](iDouble);
}

REAL Grid::getCen_y(int jj) const {
	assert(jj <= numCells[1]+B-1);
	return cellCentres[1](jj);
}

double Grid::get_y(int jj, int shift_j) const {
	int jDouble = 2*jj + shift_j;
	assert(jDouble <= numDouble[1]);
	return cellsDouble[1](jDouble);
}

// Global version of the above:
double Grid::get_y_global(int jj, int shift_j) const {
	int jDouble = 2*jj + shift_j;
	assert(jDouble <= global_numDouble[1]);
	return global_cellsDouble[1](jDouble);
}

REAL Grid::getCen_z(int kk) const {
	assert(kk <= numCells[2]+B-1);
	return cellCentres[2](kk);
}

double Grid::get_z(int kk, int shift_k) const {
	int kDouble = 2*kk + shift_k;
	assert(kDouble <= numDouble[2]);
	return cellsDouble[2](kDouble);
}

// Global version of the above:
double Grid::get_z_global(int kk, int shift_k) const {
	int kDouble = 2*kk + shift_k;
	assert(kDouble <= global_numDouble[2]);
	return global_cellsDouble[2](kDouble);
}

REAL Grid::getCen(int dir, int ii) const {
	assert(ii <= numCells[dir]+B-1);
	return cellCentres[dir](ii);
}

double Grid::get_pos(int dir, int ii, int shift_i) const {
	int iDouble = 2*ii + shift_i;
	assert(iDouble <= numDouble[dir]);
	return cellsDouble[dir](iDouble);
}

// Global version of the above:
double Grid::get_pos_global(int dir, int iPos, int shift_i) const {
	int iDouble = 2*iPos + shift_i;
	assert(iDouble <= global_numDouble[dir]);
	return global_cellsDouble[dir](iDouble);
}
 


// Number zero is left side of cell 0
// REAL Grid::getEdge_x(int ii) {
// 	assert(ii <= numCells[0]+B-1);
// 	return cellEdges[0](ii);
// }

// REAL Grid::getEdge_y(int jj) {
// 	assert(jj <= numCells[1]+B-1);
// 	return cellEdges[1](jj);
// }

// REAL Grid::getEdge_z(int kk) {
// 	assert(kk <= numCells[2]+B-1);
// 	return cellEdges[2](kk);
// }

// REAL Grid::getEdge(int dir, int ii) {
// 	assert(ii <= numCells[dir]+B-1);
// 	return cellEdges[dir](ii);
// }


// Get position of left edge of cell i
REAL Grid::getEdgL_x(int ii) const {
	assert(ii <= numCells[0]+B-1);
	return cellEdges[0](ii);
}

// Get position of right edge of cell i
REAL Grid::getEdgR_x(int ii) const {
	assert(ii <= numCells[0]+B-1);
	return cellEdges[0](ii+1);
}

REAL Grid::getEdgL_y(int jj) const {
	assert(jj <= numCells[1]+B-1);
	return cellEdges[1](jj);
}

REAL Grid::getEdgR_y(int jj) const {
	assert(jj <= numCells[1]+B-1);
	return cellEdges[1](jj+1);
}

REAL Grid::getEdgL_z(int kk) const {
	assert(kk <= numCells[2]+B-1);
	return cellEdges[2](kk);
}

REAL Grid::getEdgR_z(int kk) const {
	assert(kk <= numCells[2]+B-1);
	return cellEdges[2](kk+1);
}

REAL Grid::getEdgL(int dir, int ii) const {
	assert(ii <= numCells[dir]+B-1);
	return cellEdges[dir](ii);
}

REAL Grid::getEdgR(int dir, int ii) const {
	assert(ii <= numCells[dir]+B-1);
	return cellEdges[dir](ii+1);
}


REAL Grid::getCen_dx(int dir, int ii) const {
	assert(ii <= numCells[dir]+B-1);
	// return (cellEdges[dir](ii+1) - cellEdges[dir](ii));
	return dxCentres[dir](ii);
}

// GLOBAL version of the above:
REAL Grid::getCen_dx_global(int dir, int ii) const {
	assert(ii <= global_numCells[dir]+B-1);
	// return (cellEdges[dir](ii+1) - cellEdges[dir](ii));
	return global_dxCentres[dir](ii);
}

REAL Grid::getCen_hx(int dir, int ii) const {
	assert(ii <= numCells[dir]+B-1);
	return 0.5*getCen_idx(dir, ii);
}

REAL Grid::getCen_idx(int dir, int ii) const {
	assert(ii <= numCells[dir]+B-1);
	return idxCentres[dir](ii);
}

REAL Grid::get_dx(int dir, int ii, int shift_i) {
	int iDouble = 2*ii + shift_i;
	assert(iDouble <= numDouble[dir]);
	return dxDouble[dir](iDouble);
}

REAL Grid::get_idx(int dir, int ii, int shift_i) {
	int iDouble = 2*ii + shift_i;
	assert(iDouble <= numDouble[dir]);
	return idxDouble[dir](iDouble);
}

REAL Grid::get_hx(int dir, int ii, int shift_i) {
	return 0.5*get_idx(dir,ii,shift_i);
}

REAL Grid::get_CellVolume(int ii, int jj, int kk) const {
	REAL volume((getEdgR_x(ii) - getEdgL_x(ii))*
	            (getEdgR_y(jj) - getEdgL_y(jj))*
	            (getEdgR_z(kk) - getEdgL_z(kk)));
#if (GEOM != CARTESIAN)
	volume *= (getCen_h0(ii,jj,kk)*getCen_h1(ii,jj,kk)*getCen_h2(ii,jj,kk));
#endif
	return volume;
}



#else // Linear grid:




int Grid::get_CellIndex(int dir, REAL position) const {
	int index;
	if(position < xb[dir] || position > xe[dir]) {
		index = -99;
	} else {
		index = static_cast<int>((position - xb[dir])/dx[dir]);
	}
	return index;
}


int Grid::Shift_Rank(int dir, int ii) const {
	// return ii + coords[dir]*numCells[dir];
	// return coords[dir]*mx[dir];
	return coords[dir]*numCellsEff[dir];
}



double Grid::get_x(REAL ii) const
{
	return x0[0]+ii*dx[0];
}

double Grid::get_x(int ii, int shift_i) const
{
	return x0[0]+(ii+0.5*shift_i)*dx[0];
}

// Global version of the above
double Grid::get_x_global(int ii, int shift_i) const
{
	return x0[0]+(ii+0.5*shift_i)*dx[0];
}

double Grid::getCen_x(int ii) const {
	return x0[0]+ii*dx[0];
}

// Position of left edge of cell i
// double Grid::getEdge_x(int ii) {
// #ifdef parallel
// 	ii += Shift_Rank(0, ii);
// #endif
// 	return x0[0]+(ii-0.5)*dx[0];
// }

// Position of left edge of cell i
double Grid::getEdgL_x(int ii) const {
	return x0[0]+(ii-0.5)*dx[0];
}

double Grid::getEdgR_x(int ii) const {
	return x0[0]+(ii+0.5)*dx[0];
}

double Grid::get_y(REAL ii) const {
	return x0[1]+ii*dx[1];
}

double Grid::get_y(int iy, int shift_y) const
{
	return x0[1]+(iy+0.5*shift_y)*dx[1];
}

// Global version of the above:
double Grid::get_y_global(int iy, int shift_y) const
{
	return x0[1]+(iy+0.5*shift_y)*dx[1];
}

double Grid::getCen_y(int jj) const {
	return x0[1]+jj*dx[1];
}

// double Grid::getEdge_y(int jj) {
// #ifdef parallel
// 	jj += Shift_Rank(1, jj);
// #endif
// 	return x0[1]+(jj-0.5)*dx[1];
// }

double Grid::getEdgL_y(int jj) const {
	return x0[1]+(jj-0.5)*dx[1];
}

double Grid::getEdgR_y(int jj) const {
	return x0[1]+(jj+0.5)*dx[1];
}

double Grid::get_z(REAL ii) const {
	return x0[2]+ii*dx[2];
}

double Grid::get_z(int iz, int shift_z) const
{
	return x0[2]+(iz+0.5*shift_z)*dx[2];
}

// Global version of the above
double Grid::get_z_global(int iz, int shift_z) const
{
	return x0[2]+(iz+0.5*shift_z)*dx[2];
}

double Grid::getCen_z(int kk) const {
	return x0[2]+kk*dx[2];
}

// double Grid::getEdge_z(int kk) {
// #ifdef parallel
// 	kk += Shift_Rank(2, kk);
// #endif
// 	return x0[2]+(kk-0.5)*dx[2];
// }

// Left edge of z in 3-direction
double Grid::getEdgL_z(int kk) const {
	return x0[2]+(kk-0.5)*dx[2];
}

// Right edge of z in 3-direction
double Grid::getEdgR_z(int kk) const {
	return x0[2]+(kk+0.5)*dx[2];
}

// double Grid::get_x(int dir, REAL ii) 
// {
// 	assert((dir >=0) && (dir < DIM)); 
// #ifdef parallel
// 	ii += Shift_Rank(dir, ii);
// #endif
// 	return x0[dir]+ii*dx[dir];
// }

REAL Grid::getCen(int dir, int ii) const {
	assert(ii <= numCells[dir]+B-1);
	return x0[dir]+ii*dx[dir];
}

// linear version of general get_pos operator
double Grid::get_pos(int dir, int iPos, int shift_i) const {
	assert(iPos >= -B && iPos <= numCells[dir]+B-1);
	return x0[dir]+(iPos+0.5*shift_i)*dx[dir];
}

// global version of the above
double Grid::get_pos_global(int dir, int iPos, int shift_i) const
{
	return x0[dir]+(iPos+0.5*shift_i)*dx[dir];
}

// REAL Grid::getEdge(int dir, int ii) {
// 	assert(ii <= numCells[dir]+B-1);
// #ifdef parallel
// 	ii += Shift_Rank(dir, ii);
// #endif
// 	return x0[dir]+(ii-0.5)*dx[dir];
// }

// Left edge of cell in direction dir
REAL Grid::getEdgL(int dir, int ii) const {
	assert(ii <= numCells[dir]+B-1);
	return x0[dir]+(ii-0.5)*dx[dir];
}

// Right edge of cell in direction dir
REAL Grid::getEdgR(int dir, int ii) const {
	assert(ii <= numCells[dir]+B-1);
	return x0[dir]+(ii+0.5)*dx[dir];
}

REAL Grid::getCen_dx(int dir, int ii) const {
	return dx[dir];
}

// GLOBAL version of the above:
REAL Grid::getCen_dx_global(int dir, int ii) const {
	return dx[dir];
}

REAL Grid::getCen_hx(int dir, int ii) const {
	return hx[dir];
}

REAL Grid::getCen_idx(int dir, int ii) const {
	return idx[dir];
}

REAL Grid::get_dx(int dir, int ii, int shift_i) {
	return dx[dir];
}

REAL Grid::get_idx(int dir, int ii, int shift_i) {
	return idx[dir];
}

REAL Grid::get_hx(int dir, int ii, int shift_i) {
	return hx[dir];
}

REAL Grid::get_CellVolume(int ii, int jj, int kk) const {
	REAL volume = dx[0]*dx[1]*dx[2];
#if (GEOM != CARTESIAN)
	volume *= getCen_h0(ii,jj,kk)*getCen_h1(ii,jj,kk)*getCen_h2(ii,jj,kk);
#endif
	return volume;
}


#endif


// REAL Grid::get_GridFunction(ProblemType &Problem, int dir, REAL ratio) {
	
// 	if(parametrisation[dir] == 0) {
// 		// Linear grid:
// 		return ratio;
// 	} else if(parametrisation[dir] == 1) {
// 		// Exponential grid:
// 		return get_ExpGrid(dir, ratio);
// 	} else if(parametrisation[dir] == 4) {
// 		if(dir == 0) {
// 			return Problem.grid_user_x(ratio);
// 		} else if (dir == 1) {
// 			return Problem.grid_user_y(ratio);
// 		} else {
// 			return Problem.grid_user_z(ratio);
// 		}
// 	} else {
// 		cerr << " No such gridding type yet " << endl;
// 		exit(-24);
// 		return -24.;
// 	}
// }

// REAL Grid::get_GridFunction(int dir, REAL ratio) {

#if (NON_LINEAR_GRID == CRONOS_ON)
REAL Grid::get_GridFunction(int dir, REAL ratio, REAL (*func) (REAL)) {
	
	if(parametrisation[dir] == 0) {
		// Linear grid:
		return ratio;
	} else if(parametrisation[dir] == 1) {
		// Exponential grid:
		return get_ExpGrid(dir, ratio);
	} else if(parametrisation[dir] == 2) {
		// New logarithmic grid
//		return
	} else if (parametrisation[dir] == 3) {
		// New ratio grid (logarithmic)
//		return get_RatioGrid(dir, ratio);
	} else if(parametrisation[dir] == 4) {
		

		// if(dir == 0) {
		// 	return Problem.grid_user_x(ratio);
		// } else if (dir == 1) {
		// 	return Problem.grid_user_y(ratio);
		// } else {
		// 	return Problem.grid_user_z(ratio);
		// }
		cerr << " To be included some time " << endl;
		exit(-24);
		return -24;
	} else {
		cerr << " No such gridding type yet " << endl;
		exit(-24);
		return -24;
	}
}

REAL Grid::get_ExpGrid(int dir, REAL ratio) {
	return global_xb[dir]/Lx[dir]*
		(pow(global_xe[dir]/global_xb[dir],ratio) - 1.);
}


// Compute grid

//void Grid::generate_Grid(REAL (*func) (int, REAL)) {
void Grid::generate_Grid() {

	for(int dir=0; dir<DIM; ++dir) {

		// Compute position of cell edges -- grid is defined via those
		for(int ii=-B; ii<=numCells[dir]+B+1; ++ii) {

			// New griding - xb at cell edge:
			REAL ratio = (1.*(ii+coords[dir]*numCells[dir]))/
				(1.*global_numCells[dir]);
			if(!EdgeGridding) {
				// Revert to old gridding scheme:
				ratio -= 0.5/(1.*global_numCells[dir]);
			} 
				
			// Compute position of cell edges
			// cellEdges[dir](ii) = global_xb[dir] + 
			// 	get_GridFunction(dir, ratio, *func)*Lx[dir];
			// cout << " Doing the grid " << dir << " " << ii << " " << numCells[dir]+B+1 << " " << cellEdges[dir].getHigh(0) << " " << mx[dir]+B+1 << endl;
			cellEdges[dir](ii) = global_xb[dir] + myNonlinGrid[dir]->get_gridFunc(ratio)*Lx[dir];
			if(parametrisation[dir] == 2 || parametrisation[dir] == 3 || parametrisation[dir] >= 5) {
				cellEdges[dir](ii) = myNonlinGrid[dir]->get_gridFunc(ratio);
			}
			//				get_GridFunction(dir, ratio, *func)*Lx[dir];

		}



		// Compute position of GLOBAL cell edges:
		for(int ii=-B; ii<=global_numCells[dir]+B+1; ++ii) {

			// New griding - xb at cell edge:
			REAL ratio = (1.*(ii))/(1.*global_numCells[dir]);
			
			if(!EdgeGridding) {
				// Revert to old gridding scheme:
				ratio -= 0.5/(1.*global_numCells[dir]);
			} 

			global_cellEdges[dir](ii) = global_xb[dir] +
				myNonlinGrid[dir]->get_gridFunc(ratio)*Lx[dir];

			if(parametrisation[dir] == 2 || parametrisation[dir] == 3 || parametrisation[dir] >= 5) {
				global_cellEdges[dir](ii) = myNonlinGrid[dir]->get_gridFunc(ratio);
			}

		}


		// Compute position of Cell centres as middle between cell
		// faces
		for(int ii=-B; ii<=numCells[dir]+B; ++ii) {
			cellCentres[dir](ii) = 0.5*(cellEdges[dir](ii  ) +
			                            cellEdges[dir](ii+1));
		}
		

		// Compute position of GLOBAL cell centres:
		for(int ii=-B; ii<=global_numCells[dir]+B; ++ii) {
			global_cellCentres[dir](ii) = 0.5*(global_cellEdges[dir](ii  ) +
			                                   global_cellEdges[dir](ii+1));
		}



		// Write array holding edges and centres:
		for(int ii=-B; ii<=numCells[dir]+B; ++ii) {
			cellsDouble[dir](2*ii  ) = cellCentres[dir](ii);
			cellsDouble[dir](2*ii-1) = cellEdges[dir](ii);
		}
		cellsDouble[dir](2*(numCells[dir]+B)+1) = cellEdges[dir](numCells[dir]+B + 1);

		// Write array holding GLOBAL edges and centres:
		for(int ii=-B; ii<=global_numCells[dir]+B; ++ii) {
			global_cellsDouble[dir](2*ii  ) = global_cellCentres[dir](ii);
			global_cellsDouble[dir](2*ii-1) = global_cellEdges[dir](ii);
		}

		// Compute List of dx-Values
		for(int ii=-B; ii<=numCells[dir]+B; ++ii) {
			dxCentres[dir](ii) = cellEdges[dir](ii+1)-cellEdges[dir](ii);
			idxCentres[dir](ii) = 1./dxCentres[dir](ii);

			dxDouble[dir](2*ii) = dxCentres[dir](ii);
			idxDouble[dir](2*ii) = 1./dxDouble[dir](ii);

			if(ii<numCells[dir]+B) {
				dxDouble[dir](2*ii+1) = (cellCentres[dir](ii+1) -
				                         cellCentres[dir](ii));
				idxDouble[dir](2*ii+1) = 1./dxDouble[dir](2*ii+1);
			}
		}

		// Compute List of GLOBAL dx values:
		for(int ii=-B; ii<=global_numCells[dir]+B; ++ii) {
			global_dxCentres[dir](ii) = global_cellEdges[dir](ii+1) -
				global_cellEdges[dir](ii);
		}

	}

	// Reset array bounds according to parametrisation (only necessary
	// locally)
	for(int dir=0; dir<DIM; ++dir) {
		xb[dir] = getEdgL(dir, 0);
		xe[dir] = getEdgR(dir, mx[dir]);
	}

}
#endif


/*
  Definition of different geometries
  Type           Direction          Meaning          Geom factor
  Cartesian:        0                 x                h0 = 1
                    1                 y                h1 = 1
                    2                 z                h2 = 1

  Plane Polar       0                 r                h0 = 1
                    1                phi               h1 = r
                    2                 z                h2 = 1

  Spherical         0                 r                h0 = 1
                    1                theta             h1 = r
                    2                phi               h2 = r sin theta
*/

#ifdef GEOM


REAL Grid::get_CellGeomTrafo(int ii, int jj, int kk) {
	//! Returns geometrical factors for cell centre
	return getCen_h0(ii,jj,kk)*getCen_h1(ii,jj,kk)*getCen_h2(ii,jj,kk);
}


REAL Grid::get_AreaGeomTrafo(int dir, int ii, int jj, int kk) {
	//! Returns geometrical factors for left-handed cell face
	if(dir==0) {
		return get_AreaGeomTrafo_x(ii, jj, kk);
	} else if (dir==1) {
		return get_AreaGeomTrafo_y(ii, jj, kk);
	} else {
		return get_AreaGeomTrafo_z(ii, jj, kk);
	}

	return getCen_h0(ii,jj,kk)*getCen_h1(ii,jj,kk)*getCen_h2(ii,jj,kk);
}


REAL Grid::get_AreaGeomTrafo_x(int ii, int jj, int kk) {
	//! Returns geometrical factors for left-handed cell face in x-direction
	return getEdgL_h1(ii,jj,kk)*getEdgL_h2(ii,jj,kk);
}

REAL Grid::get_AreaGeomTrafo_y(int ii, int jj, int kk) {
	//! Returns geometrical factors for left-handed cell face in y-direction
	return getEdgL_h0(ii,jj,kk)*getEdgL_h2(ii,jj,kk);
}

REAL Grid::get_AreaGeomTrafo_z(int ii, int jj, int kk) {
	//! Returns geometrical factors for left-handed cell face in z-direction
	return getEdgL_h0(ii,jj,kk)*getEdgL_h1(ii,jj,kk);
}


REAL Grid::get_CellArea_x(int ii, int jj, int kk) {
	//! Return area of left-handed cell face in x-direction
#if (NON_LINEAR_GRID == CRONOS_OFF)
	REAL Area = dx[1]*dx[2];
#else
	REAL Area = getCen_dx(1,jj)*getCen_dx(2,kk);
#endif
	Area *= getEdgL_h1(ii,jj,kk)*getEdgL_h2(ii,jj,kk);
	return Area;
}

REAL Grid::get_CellArea_y(int ii, int jj, int kk) {
	//! Return area of left-handed cell face in x-direction
#if (NON_LINEAR_GRID == CRONOS_OFF)
	REAL Area = dx[0]*dx[2];
#else
	REAL Area = getCen_dx(0,ii)*getCen_dx(2,kk);
#endif
	Area *= getEdgL_h0(ii,jj,kk)*getEdgL_h2(ii,jj,kk);
	return Area;
}

REAL Grid::get_CellArea_z(int ii, int jj, int kk) {
	//! Return area of left-handed cell face in x-direction
#if (NON_LINEAR_GRID == CRONOS_OFF)
	REAL Area = dx[0]*dx[1];
#else
	REAL Area = getCen_dx(0,ii)*getCen_dx(1,jj);
#endif
	Area *= getEdgL_h0(ii,jj,kk)*getEdgL_h1(ii,jj,kk);
	return Area;
}




#if (NON_LINEAR_GRID == CRONOS_OFF)
double Grid::h0(REAL ii, REAL jj, REAL kk) const {
#if (GEOM == 1) // Cartesian
	return 1.;
#endif
#if (GEOM == 2) // plane polar
	return 1.;
#endif
#if (GEOM == 3) // spherical
	return 1.;
#endif
}
#endif


double Grid::h0(int ii, int jj, int kk, int shift_i, int shift_j, int shift_k) const {
#if (GEOM == 1)   // Cartesian
	return 1.;
#elif (GEOM == 2) // plane polar
	return 1.;
#elif (GEOM == 3) // spherical
	return 1.;
#endif
}

double Grid::getCen_h0(int ii, int jj, int kk) const {
#if (GEOM == 1)   // Cartesian
	return 1.;
#elif (GEOM == 2) // plane polar
	return 1.;
#elif (GEOM == 3) // spherical
	return 1.;
#endif
}

double Grid::getEdgL_h0(int ii, int jj, int kk) const {
#if (GEOM == 1)   // Cartesian
	return 1.;
#elif (GEOM == 2) // plane polar
	return 1.;
#elif (GEOM == 3) // spherical
	return 1.;
#endif
}


#if (NON_LINEAR_GRID == CRONOS_OFF)
double Grid::h1(REAL ii, REAL jj, REAL kk) const {
#if (GEOM == 1) // Cartesian
	return 1.;
#endif
#if (GEOM == 2) // plane polar
	// return (x0[0]+ii*dx[0]);
	return (get_x(ii));
#endif
#if (GEOM == 3) // spherical
	return (get_x(ii));
	//  return (0.5*(xb[0]+ii*dx[0] + xe[0]-(mx[0]-ii)*dx[0]));
#endif
}
#endif

double Grid::h1(int ii, int jj, int kk, int shift_i, int shift_j, int shift_k) const {
#if (GEOM == 1)   // Cartesian
	return 1.;
#elif (GEOM == 2) // plane polar
	return (get_x(ii, shift_i));
#elif (GEOM == 3) // spherical
	return (get_x(ii, shift_i));
#endif
}

double Grid::getCen_h1(int ii, int jj, int kk) const {
#if (GEOM == 1)   // Cartesian
	return 1.;
#elif (GEOM == 2) // plane polar
	return (getCen_x(ii));
#elif (GEOM == 3) // spherical
	return (getCen_x(ii));
#endif
}


double Grid::getEdgL_h1(int ii, int jj, int kk) const {
#if (GEOM == 1)   // Cartesian
	return 1.;
#elif (GEOM == 2) // plane polar
	return (getEdgL_x(ii));
#elif (GEOM == 3) // spherical
	return (getEdgL_x(ii));
#endif
}


#if (NON_LINEAR_GRID == CRONOS_OFF)

double Grid::h2(REAL ii, REAL jj, REAL kk) const {
#if (GEOM == 1)
	return 1.;
#endif
#if (GEOM == 2)
	return 1.;
	//  return (xb[0]+ii*dx[0]);
#endif
#if (GEOM == 3)
	//return (xb[0]+ii*dx[0])*sin(xb[1]+jj*dx[1]);
	// return (x0[0]+ii*dx[0])*sin(0.5*(x0[1]+jj*dx[1] +
	//                                  xe[1]-(mx[1]-jj)*dx[1]));
	return (get_x(ii))*sin(get_y(jj));
#endif
}
#endif

double Grid::h2(int ii, int jj, int kk, int shift_i, int shift_j, int shift_k) const {
#if (GEOM == 1)   // Cartesian
	return 1.;
#elif (GEOM == 2) // plane polar
	return 1.;
#elif (GEOM == 3) // spherical
	return (get_x(ii, shift_i))*sin(get_y(jj, shift_j));
#endif
}

double Grid::getCen_h2(int ii, int jj, int kk) const {
#if (GEOM == 1)   // Cartesian
	return 1.;
#elif (GEOM == 2) // plane polar
	return 1.;
#elif (GEOM == 3) // spherical
	return (getCen_x(ii))*sin(getCen_y(jj));
#endif
}

double Grid::getEdgL_h2(int ii, int jj, int kk) const {
#if (GEOM == 1)   // Cartesian
	return 1.;
#elif (GEOM == 2) // plane polar
	return 1.;
#elif (GEOM == 3) // spherical
	return (getEdgL_x(ii))*sin(getEdgL_y(jj));
#endif
}

#endif


void Grid::ScaleGrid(NumArray<double>& xb, NumArray<double>& xe, const int &dim, const REAL &scale) {
	xb[dim] *= scale;
	xe[dim] *= scale;
	//   dx[DIM] *= scale;
	//   Lx[DIM] *= scale;
}


void Grid::CheckScale(NumArray<double>& xb, NumArray<double>& xe) {
	string val;
  
	// x-direction
	val = svalue((char*)"scale_x");
	if(val == "pi" || val == "Pi") {
		ScaleGrid(xb, xe, 0, M_PI);
	} else  if (val == "e") {
		ScaleGrid(xb, xe, 0, exp(1.));
	}

	// y-direction
	val = svalue((char*)"scale_y");
	if(val == "pi" || val == "Pi") {
		ScaleGrid(xb, xe, 1, M_PI);
	} else  if (val == "e") {
		ScaleGrid(xb, xe, 1, exp(1.));
	}

	// z-direction
	val = svalue((char*)"scale_z");
	if(val == "pi" || val == "Pi") {
		ScaleGrid(xb, xe, 2, M_PI);
	} else  if (val == "e") {
		ScaleGrid(xb, xe, 2, exp(1.));
	}
}

void Grid::MakeGridDir() {

	// Write grid output into directory
	sprintf(griddirname, "%s/%s_grid", getenv("poub"),getenv("pname"));
	for (char *tmp=griddirname+strlen(griddirname)-1; *tmp=='0'; tmp--) *tmp=0;
	
	// struct stat st;
	// if(stat(griddirname,&st) != 0) {

		// Make directory
		//mkdir(griddirname, 511);
	filesystem::create_directory(griddirname);
	filesystem::permissions(griddirname, filesystem::perms::owner_read | filesystem::perms::owner_write | filesystem::perms::owner_exec | filesystem::perms::group_exec | filesystem::perms::others_exec);
	auto ftime = filesystem::last_write_time(griddirname);
	filesystem::file_time_type currentTime(decltype(ftime)::clock::now());
	filesystem::last_write_time(griddirname, currentTime);

	// } else {
		    
	// 	printf(" /tmp is present\n");

	// }
}


void Grid::WriteGridInfo() {

	// Produce information output for each process

	// Write grid output into directory
	MakeGridDir();

	ofstream fout;
	char fname[1024];
	sprintf(fname,"%s/%s_grid.info",griddirname,getenv("pname"));
	fout.open(fname);

	fout << " xbeg[0] = " << xb[0]
	     << " xbeg[1] = " << xb[1]
	     << " xbeg[2] = " << xb[2] << endl;
	fout << " xend[0] = " << xe[0]
	     << " xend[1] = " << xe[1]
	     << " xend[2] = " << xe[2] << endl;
	fout << " dx[0] = " << dx[0]
	     << " dx[1] = " << dx[1]
	     << " dx[2] = " << dx[2] << endl;

	fout.close();

}


bool Grid::WriteGridStructure() {
	#ifndef UNIT_TEST
	/*
	  Write the actual grid structure. Writing nodes, cells and
	  transformations
	*/

	// The grid is currently written only by main CPU - this is NOT
	// suited for very large computations due to memory problems and
	// writing efficiency. For the future we plan to include a
	// parallel output facility for the data - this, however, has to
	// be handled with care (i.e., might not work) on NFS mounted file
	// systems.
	if(rank == 0) {

		// // Make directory
		// MakeGridDir();
		// Open hdf5 file
		string gridfilename = string(griddirname) +"/"+ string(getenv("pname")) +"_grid.h5";
		Hdf5Stream gridh5file(gridfilename, 18, rank);

		// Write descriptive attribute to file:
		string geomName[] = { "Cartesian"  ,
		                      "Cylindrical",
		                      "Spherical"  };

		// Indicate usage of non-linear grid:
#if (NON_LINEAR_GRID == CRONOS_ON)
		geomName[GEOM-1] += "_nonlin";
#endif
		gridh5file.AddGlobalAttr("Geom_Type", geomName[GEOM-1]);

	
		// Builidng matrix-arrays to hold location data:
		string NameNodes[3] = { "/Xn", "/Yn", "/Zn" };
		string NameCenters[3] = { "/Xc", "/Yc", "/Zc" };
		string NameTrafo[3] = { "/Tx", "/Ty", "/Tz" };
		string NameLines[3] = { "/Lx", "/Ly", "/Lz" };
		string NameLineCentres[3] = { "/xCen", "/yCen", "/zCen" };
		string NameLineEdges[3] = { "/xL", "/yL", "/zL" };
		NumMatrix<float,3> CartNodes(Index::set(0,0,0),
		                             Index::set(global_mx[0]+1,
		                                        global_mx[1]+1,
		                                        global_mx[2]+1));
		NumMatrix<float,3> CartCenters(Index::set(0,0,0),
		                               Index::set(global_mx[0],
		                                          global_mx[1],
		                                          global_mx[2]));
		NumMatrix<float,1> curvLines;  // size depends on dimension

		NumMatrix<float,1> lineCentres, lineEdgesL;

		// NumMatrix<float,3> TrafoVec[3];
		// for(int dim=0; dim<3; ++dim) {
		// 	TrafoVec[dim].resize(Index::set(0,0,0), Index::set(mx[0],mx[1],mx[2]));
		// }
		// No nice 4D matrix yet...
		int mxTrafoVec[] = {global_mx[2]+1, global_mx[1]+1, global_mx[0]+1, 3 };
		int sizeTrafoVec((global_mx[0]+1)*(global_mx[1]+1)*(global_mx[2]+1)*3);
 		float * TrafoVec = new float [sizeTrafoVec];
		
		// Compute locations for different directions.
		for (short did = 0; did < 3; ++did) {     // loop {X,Y,Z} direction
			// cn, cc, TT: arrays for current component (no ghost cells!) for...
			// Cartesian nodes (cn), Cartesian centers (cc),
			// transformation "vectors" (TT)

			// mimick 'cellsDouble' arrays (but note B -> BOUT_FLT!)
			const int ci_ran[] = {-2*                BOUT_FLT -1,
			                      +2*(global_mx[did]+BOUT_FLT)+1};
			// const int mx_ci = ci_ran[1]-ci_ran[0]+1;
			curvLines.resize(Index::set(ci_ran[0]),
			                 Index::set(ci_ran[1]));

			lineCentres.resize(Index::set(-BOUT_FLT),
					Index::set(global_mx[did]+BOUT_FLT));
			lineEdgesL.resize(Index::set(-BOUT_FLT),
					Index::set(global_mx[did]+BOUT_FLT+1));

			// Compute 1D coords at dx/2 intervals
			for (int ci = ci_ran[0]; ci <= ci_ran[1]; ci++) {
				curvLines(ci) = get_pos_global(did, ci/2, ci%2);
				// this should equal cellsDouble[did](ci) for nonlinear case
			}

			// Store global cell centres
			for (int iPos=BOUT_FLT; iPos<=global_mx[did]+BOUT_FLT; ++iPos) {
				lineCentres(iPos) = get_pos_global(did, iPos, 0);
			}
			// Store global cell edges (left-handed ones)
			for (int iPos=BOUT_FLT; iPos<=global_mx[did]+BOUT_FLT+1; ++iPos) {
				lineEdgesL(iPos) = get_pos_global(did, iPos, -1);
			}

			for (short gid = 0; gid <= 1; ++gid) { // two grids: 0 for cells,
				//              1 for nodes (to avoid repeating GEOM equations)
				
				// Loop over the full grid:
				int ndx = 0; // Running index
				for (int k = 0; k <= global_mx[2]+gid; k++) {        // nodes are offset from
					//	float c3 = get_z (k - 0.5*gid);  //  centers by 1/2
					// float c3 = get_z(k, -gid);  //  centers by 1/2
					float c3 = get_z_global(k, -gid);  //  centers by 1/2
					for (int j = 0; j <= global_mx[1]+gid; j++) {
						// float c2 = get_y (j - 0.5*gid);
						// float c2 = get_y (j, -gid);
						float c2 = get_y_global (j, -gid);
						for (int i = 0; i <= global_mx[0]+gid; i++) {
							// float c1 = get_x (i - 0.5*gid);
							// float c1 = get_x (i, -gid);
							float c1 = get_x_global (i, -gid);

#if (GEOM == 1)
							// Cartesian coordinates
							float pos_local[] = { c1, c2, c3 };    // C=[x,y,z]
							float Tt[] = { 1, 0, 0,
							               0, 1, 0,
							               0, 0, 1  };
#elif (GEOM == 2)
							// Cyldindrical grid
							float pos_local[] = { c1 * cos(c2),    // C=[r,phi,z]
							                      c1 * sin(c2),
							                      c3            };
							float Tt[] = { cos(c2), -sin(c2), 0,
							               sin(c2),  cos(c2), 0,
							               0      ,  0      , 1  };
#elif (GEOM == 3)
							// Spherical grid
							float pos_local[] = { c1 * sin(c2) * cos(c3),  // C=[r,theta,phi]
							                      c1 * sin(c2) * sin(c3),
							                      c1 * cos(c2)            };
							float Tt[] = { sin(c2)*cos(c3),  cos(c2)*cos(c3), -sin(c3),
							               sin(c2)*sin(c3),  cos(c2)*sin(c3),  cos(c3),
							               cos(c2)        , -sin(c2)        ,  0        };
#else
							cerr << "error: GEOM must be from {1,2,3}, is: " << GEOM << endl;
							exit(-62);
#endif
							if (gid == 0) {
								CartCenters(i,j,k) =  pos_local[did]; // build centers from c_temp
								for (int coid = 0; coid < 3; ++coid) {// TT is at centers
									assert((3*ndx)+coid < sizeTrafoVec);
									TrafoVec[ (3*ndx)+coid ]  = Tt [ (3*did)+coid ];
								}
							} else {
								CartNodes(i,j,k) = pos_local[did]; // build nodes from c_temp
							}
							++ndx;
						}  // end i loop
					}     // end j loop
				}        // end k loop
			}         // end gid loop
			// END { build arrays }
		
			// BEGIN { write arrays to HDF5 file }
		
			// --> write 1d coord marks:
			cout << " Writing Lines: " << NameLines[did] << endl;
			gridh5file.Write1DMatrix(NameLines[did], curvLines);
			gridh5file.Write1DMatrix(NameLineCentres[did], lineCentres);
			gridh5file.Write1DMatrix(NameLineEdges[did], lineEdgesL);
// int mx_ci = ci_ran[1]-ci_ran[0]+1;
//   gridh5file.WriteNDArray(NameLines[did], curvLines, &mx_ci, 1);

			// --> write nodes:
			cout << " Writing Nodes: " << NameNodes[did] << endl;
			gridh5file.Write3DMatrix(NameNodes[did], CartNodes);

			// --> write cells:
			cout << " Writing Cells: " << NameCenters[did] << endl;
			gridh5file.Write3DMatrix(NameCenters[did], CartCenters);

			// --> write trafo coefficients
			cout << " Writing Multi " << NameTrafo[did] << endl;
			gridh5file.WriteNDArray(NameTrafo[did], TrafoVec,
			                        mxTrafoVec, 4);
				//			gridh5file.WriteArray(NameTrafo[did], TrafoVec, sizeTrafoVec);
			cout << " Done " << endl;

		}
 		delete[] TrafoVec;
		gridh5file.close();
// 		cout << "Grid written to file " << gridfilename << endl;
	}
	return true;

#endif
}


bool Grid::WriteGridData() {
	#ifndef UNIT_TEST
	/*************************************************************
       Writing the full grid into output file.
	 ************************************************************/

	NumMatrix<float,3> gridcells[DIM];
	int lbound[DIM], ubound[DIM];
	// Setting size of grid-matrix:
	int rim = BOUT_FLT;
	for(int i=0;i<3;i++){
		lbound[i]=-rim;
		ubound[i]=mx[i]+rim;
	}
	for(int i=0; i<DIM; ++i) {
		gridcells[i].resize(lbound,ubound);
	}
	
	// Compute position of centre for all grid cells
	for(int k = -rim; k <= mx[2]+rim; k++) {
		for (int j = -rim; j <= mx[1]+rim; j++) {
			for (int i = -rim; i <= mx[0]+rim; i++) {
				gridcells[0](i,j,k) = getCen_x(i);
				gridcells[1](i,j,k) = getCen_y(j);
				gridcells[2](i,j,k) = getCen_z(k);
			}
		}
	}
	
	string filename = getenv("poub");
	filename += "/";
	filename += getenv("pname");
	filename += "_grid";

	filename += ".h5";
	Hdf5Stream h5out(filename, 1, rank);
	h5out.Write3DVecMatrix("Grid", gridcells[0], gridcells[1], gridcells[2]);

	// Reduce size of grid-matrix
	for(int i=0;i<3;i++){
		lbound[i]= 0;
		ubound[i]= 0;
	}
	for(int i=0; i<DIM; ++i) {
		gridcells[i].resize(lbound,ubound);
	}
	return true;

#endif
}

#if (NON_LINEAR_GRID == CRONOS_ON)
void Grid::add_nonLinGridToHdf(Hdf5Stream &h5out, int num_ghost) {
	//!
	/* h5out hdf5 output file to add grid properties to
	 * num_ghost number of ghost cells
	 * */

	// We start by storing the grid
	int num[3];
	for(int iDir=0; iDir<3; ++iDir) {
		num[iDir] = 2*num_ghost + global_numCells[iDir];
	}

	NumMatrix<float,1> xCen, yCen, zCen;
	NumMatrix<float,1> xL, yL, zL;

	xCen.resize(Index::set(-num_ghost), Index::set(global_mx[0]+num_ghost));
	yCen.resize(Index::set(-num_ghost), Index::set(global_mx[1]+num_ghost));
	zCen.resize(Index::set(-num_ghost), Index::set(global_mx[2]+num_ghost));

	xL.resize(Index::set(-num_ghost), Index::set(global_mx[0]+num_ghost+1));
	yL.resize(Index::set(-num_ghost), Index::set(global_mx[1]+num_ghost+1));
	zL.resize(Index::set(-num_ghost), Index::set(global_mx[2]+num_ghost+1));


	for(int ix=-num_ghost; ix<=global_mx[0]+num_ghost; ++ix) {
		xCen(ix) = global_cellCentres[0](ix);
	}
	for(int ix=-num_ghost; ix<=global_mx[0]+num_ghost+1; ++ix) {
		xL(ix) = global_cellEdges[0](ix);
	}

	for(int iy=-num_ghost; iy<=global_mx[1]+num_ghost; ++iy) {
		yCen(iy) = global_cellCentres[1](iy);
	}
	for(int iy=-num_ghost; iy<=global_mx[1]+num_ghost+1; ++iy) {
		yL(iy) = global_cellEdges[1](iy);
	}

	for(int iz=-num_ghost; iz<=global_mx[2]+num_ghost; ++iz) {
		zCen(iz) = global_cellCentres[2](iz);
	}
	for(int iz=-num_ghost; iz<=global_mx[2]+num_ghost+1; ++iz) {
		zL(iz) = global_cellEdges[2](iz);
	}

	// Now write all arrays to h5-File

	// First add a grid group
	hid_t group = h5out.AddGroup("Data/NonlinGrid");

	// Now write grid stuff
	h5out.Write1DMatrix("xCentres", xCen, group);
	h5out.Write1DMatrix("xEdgesL", xL, group);

	h5out.Write1DMatrix("yCentres", yCen, group);
	h5out.Write1DMatrix("yEdgesL", yL, group);

	h5out.Write1DMatrix("zCentres", zCen, group);
	h5out.Write1DMatrix("zEdgesL", zL, group);

	// Now close the group
	h5out.CloseGroup(group);




}
#endif



REAL Grid::get_Volume() {
	return ((global_xe[0] - global_xb[0])*
	        (global_xe[1] - global_xb[1])*
	        (global_xe[2] - global_xb[2]));
}

int Grid::axis_treatment(int dir) const {
#if (GEOM == CARTESIAN) 
	return 0.;
#endif
#if (GEOM == CYLINDRICAL)
	if(dir == 1) {
		return include_coordAxisPhi;
	} else {
		return 0;
	}
#elif (GEOM == SPHERICAL)
	if(dir == 2) {
		return include_coordAxisPhi;
	} else {
		return 0;
	}
#endif
}


int Grid::get_singularity_treatment(int dir, int top) const {
	int bound = top + 2*dir;
	return get_singularity_treatment(bound);
}

int Grid::get_singularity_treatment(int bound) const {

	return include_Singularity[bound];

}