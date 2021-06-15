#ifndef CRONOS_MHD_BLAST_3D_H
#define CRONOS_MHD_BLAST_3D_H = 1

#include "problem.H"

// using namespace std;


#if(GEOM < 3)

class MHDBlast3D: public ProblemType {
public:
	MHDBlast3D(Data &);
	void init_fields(Data &, int [], int []);
// #ifdef parallel
// 	void get_velo(Data &gdata);
// #endif
private:
	REAL InnerPressure, OuterPressure;
	REAL Blast_Radius, Blast_Volume;
	REAL Mag0;
	RandomNumberGenerator* random;
};


MHDBlast3D::MHDBlast3D(Data &gdata) : ProblemType(gdata)
{
	name = " 3D MHD blast wave ";
#if (GEOM == 1)
	name += "(Cartesian) ";
#elif (GEOM == 2)
	name += "(Cylindrical) ";
#endif
#if (GEOM > 2)
	cerr << " Error: Supposed to be a cylindrical or Cartesian problem ";
	cerr << endl;
	exit(2);
#endif
	this->OuterPressure = value((char*)"OuterPressure");
	this->InnerPressure = value((char*)"InnerPressure");
	this->Blast_Radius = 3.5*gdata.dx[0];
	this->Blast_Radius = value((char*)"Blast_Radius");
	this->Mag0 = value((char*)"Magnetic_Field")/(sqrt(4.*pi));
	this->Mag0 = value((char*)"Magnetic_Field");

}

// #ifdef parallel
// void MHDBlast3D::get_velo(Data &gdata) {
// 	//! Get velo on the line y==0. and z==0.

// 	bool on_XAxis(false);
// 	if(gdata.get_y(0.) <= 0. && 0 < gdata.get_y(gdata.mx[1]) &&
// 	   gdata.get_z(0.) <= 0. && 0 < gdata.get_z(gdata.mx[2]) ){
// 		on_XAxis = true;
// 	}

// 	if(on_XAxis) {
// 		int nproc_x = gdata.nproc[0];

// 		// Local data
// 		NumMatrix<double,1> velo_own;
// 		velo_own.resize(Index::set(0), Index::set(gdata.mx[0]));

// 		// Find index for cell holding x-Axis
// 		int iy = gdata.get_CellIndex(1, 0.);
// 		int iz = gdata.get_CellIndex(2, 0.);
// 		// Assign local data:
// 		for(int ipos=0; ipos<=gdata.mx[0]; ++ipos) {
// 			velo_own(ipos) = gdata.om[q_sx](ipos, iy, iz);
// 		}

// 		// one array for each rank in x-direction
// 		NumMatrix<double,1> velo_other[nproc_x];
// 		for(int irank=0; irank<nproc_x; ++irank) {
// 			velo_other[irank].resize(Index::set(0), Index::set(gdata.mx[0]));
// 			velo_other[irank].clear();
// 		}
		
// 		// Size of array for single rank
// 		int size = (gdata.om[q_sx].getHigh(0) -
// 		            gdata.om[q_sx].getLow(0) + 1);
		
// 		// Global data array:
// 		NumMatrix<double,1> velo_global;
// 		velo_global.resize(Index::set(0), Index::set(gdata.global_mx[0]));
		
// 		int num = 2*(gdata.nproc[0]); // Communication amount
// 		MPI_Request requests[num];
// 		MPI_Status statusrl[num];
		
// 		for (int i=0; i<num; i++){
// 			requests[i] = MPI_REQUEST_NULL;
// 		}

// 		// Loop over all processes at the same y and z position
// 		int number=0;
// 		for(int irank=0; irank<nproc_x; ++irank) {
			
// 			// Get rank of other process
// 			int rankother = gdata.rankpos(irank,gdata.coords[1],
// 			                              gdata.coords[2]);


// 			if(gdata.rank != rankother) {

// 				int tag = rankother; // message tag
				
// 				// Receive data:
// 				MPI_Irecv((double *)velo_other[irank], size, MPI_DOUBLE,
// 				          rankother, tag, gdata.comm3d, &requests[number]);
// 				number++;

// 				tag = gdata.rank;
// 				// Send data:
// 				MPI_Isend((double *)velo_own, size, MPI_DOUBLE, rankother, tag,
// 				          gdata.comm3d, &requests[number]);
// 				number++;

// 			} else {
				
// 				// Use local data:
// 				for(int ipos=0; ipos<=gdata.mx[0]; ++ipos) {
// 					velo_other[irank](ipos) = velo_own(ipos);
// 				}
				
// 			}
			
// 		}
		
// 		// Wait for communication to complete
// 		MPI_Waitall(num, requests, statusrl);

// 		// Save global data
// 		for(int irank=0; irank<nproc_x; ++irank) {
// 			int shift = irank*gdata.get_RankWidth(0);
			
// 			for(int ipos=0; ipos<=gdata.mx[0]; ++ipos) {
// 				velo_global(ipos+shift) = velo_other[irank](ipos);
// 			}
// 		}
		
// 	}
// 	MPI_Barrier(gdata.comm3d);
// }
// #endif



void MHDBlast3D::init_fields(Data &gdata, int ibeg[3], int iend[3]) 
{

	for(int k=ibeg[2]; k<=iend[2]; ++k) {
		REAL zval = gdata.getCen_z(k);
		for(int j=ibeg[1]; j<=iend[1]; ++j) {
#if (GEOM == 1)
			REAL yval = gdata.getCen_y(j);
#elif (GEOM == 2)
			REAL phi = gdata.getCen_y(j);
#endif
			for(int i=ibeg[0]; i<=iend[0]; ++i) {
#if (GEOM == 1)
				REAL xval = gdata.getCen_x(i);
				REAL r_cyl = sqrt(sqr(xval) + sqr(yval));
				REAL phi = atan2(yval, xval);
#elif (GEOM == 2)
				REAL r_cyl = gdata.getCen_x(i);
#endif
				REAL r_sph = sqrt(sqr(r_cyl) + sqr(zval));


				gdata.om[q_rho](i,j,k) = rho0;
				gdata.om[q_sx](i,j,k) = 0.;
				gdata.om[q_sy](i,j,k) = 0.;
				gdata.om[q_sz](i,j,k) = 0.;
				gdata.om[q_Eges](i,j,k) = OuterPressure;

				if(r_sph <= Blast_Radius) {
					gdata.om[q_Eges](i,j,k) = InnerPressure;
				}
				gdata.om[q_Eges](i,j,k) *= 1./(gamma-1.);
				
				// Vector potential:
				gdata.om[q_Bx](i,j,k) = 0.;
				gdata.om[q_By](i,j,k) = 0.;
			}
		}
	}

	// z-component of Vector potential -- shift taken into account!
	for(int k=ibeg[2]; k<=iend[2]; ++k) {
		for(int j=ibeg[1]; j<=iend[1]; ++j) {
			for(int i=ibeg[0]; i<=iend[0]; ++i) {
#if (GEOM == 1)
				REAL yval = gdata.getEdgL_y(j);

#elif (GEOM == 2)
				REAL phi = gdata.getEdgL_y(j);
				REAL r_cyl = gdata.getEdgL_x(i);
				REAL yval = r_cyl*sin(phi);
#endif
				gdata.om[q_Bz](i,j,k) = Mag0*yval;

			}
		}
	}


	// cout << " Mag0 " << Mag0 << endl;
	// exit(2);

	gdata.om[q_Bx].rename("A_x");
	gdata.om[q_By].rename("A_y");
	gdata.om[q_Bz].rename("A_z");

}


#endif

#endif
