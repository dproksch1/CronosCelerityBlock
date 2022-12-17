#ifndef CRONOS_COMPUTE_STEP_H
#define CRONOS_COMPUTE_STEP_H

#include "reconst_block.H"
#include "transformations.H"
#include "RiemannSolverHD.H"
#include "PhysFluxesHD.H"
#include "utils.H"

namespace gpu {

    template<typename AccessorType, typename NumFluxAccessorType, typename NumPtotalAccessorType, typename IdxConstType, typename FluidConstType>
    void computeStep(const AccessorType& om_acc, size_t ix, size_t iy, size_t iz, int q, double* cfl_lin, NumFluxAccessorType& numFlux,
                    NumPtotalAccessorType& num_ptotal, const bool thermal, const double gamma, const double problem_cs2, const double denominator,
                    const double half_beta, const int fluidType, IdxConstType idx, FluidConstType fluidConst) {

        //no vector in device code
        double uPri[gpu::FaceMax][N_OMINT] = {};
        double uCon[gpu::FaceMax][N_OMINT] = {};
        double physFlux[gpu::FaceMax][N_OMINT] = {};
        double physValPtherm[gpu::FaceMax] = {};
        double physValPtotal[gpu::FaceMax] = {};

        double uPriOld[gpu::FaceMax][N_OMINT] = {};

        //if (ix == ixEnd && iy == gdata.mx[1]+n_ghost[1]-1 && iz == gdata.mx[2]+n_ghost[2]-1) cout << "reach limit: " << ix << "." <<  iy << "." << iz << "\n";
        //const int fluidType = Riemann[DirX]->get_Fluid_Type();

        //pointwise reconstruction (incomplete)
        gpu::compute(om_acc, uPri, ix, iy, iz, q, DirAll);
        
        if (ix >= 2 && iy >= 2 && iz >= 2) {

            gpu::compute(om_acc, uPriOld, ix -1, iy, iz, q, DirX);
            gpu::compute(om_acc, uPriOld, ix, iy -1, iz, q, DirY);
            gpu::compute(om_acc, uPriOld, ix, iy, iz -1, q, DirZ);

            for (int i = 0; i < N_OMINT; i++) {
                uPri[gpu::FaceEast][i] = uPri[gpu::FaceEast][i];
                uPri[gpu::FaceNorth][i] = uPri[gpu::FaceNorth][i];
                uPri[gpu::FaceTop][i] = uPri[gpu::FaceTop][i];
            }
        
        }
        
        for (int dir = 0; dir < DirMax; ++dir) {
            int ixOff = (dir == DirX) ? ix : ix - 1;
            int iyOff = (dir == DirY) ? iy : iy - 1;
            int izOff = (dir == DirZ) ? iz : iz - 1;

            int face = dir * 2;

            gpu::get_Cons(om_acc, uPri[face], uCon[face], physValPtherm[face], physValPtotal[face], ix, iy, iz, face,
                            thermal, gamma, problem_cs2, denominator, half_beta, fluidType, fluidConst);

            gpu::get_PhysFlux(om_acc, physFlux[face], uPri[face], uCon[face], physValPtherm[face], face, fluidConst);

            gpu::get_Cons(om_acc, uPri[face+1], uCon[face+1], physValPtherm[face+1], physValPtotal[face+1], ix, iy, iz, face+1,
                            thermal, gamma, problem_cs2, denominator, half_beta, fluidType, fluidConst);

            gpu::get_PhysFlux(om_acc, physFlux[face+1], uPri[face+1], uCon[face+1], physValPtherm[face+1], face+1, fluidConst);
        }

        double numVals_Ch[DirMax][gpu::NumV_Max] = {1., 1.};

        if (ix >= 2 && iy >= 2 && iz >= 2) {
            for (int dir = 0; dir < DirMax; ++dir) {
                int face = dir * 2;
                *cfl_lin = cl::sycl::fmax(gpu::get_vChar(om_acc, uPri[face], uCon[face], physValPtherm[face], uPri[face+1], uCon[face+1],
                                physValPtherm[face+1], numVals_Ch[dir], dir, gamma, idx, fluidConst), *cfl_lin);
                num_ptotal[dir] = gpu::get_NumFlux(om_acc, uPri[face+1], uCon[face+1], physFlux[face+1], physValPtherm[face+1],
                                        physValPtotal[face+1], uPri[face], uCon[face], physFlux[face], physValPtherm[face], 
                                        physValPtotal[face], numFlux[dir], numVals_Ch[dir], dir, gamma, fluidConst);
            }
        }
    }

    template<typename AccessorType, typename NumFluxAccessorType, typename NumPtotalAccessorType, typename IdxConstType, typename FluidConstType, typename OutputAccType>
    void computeStep(const AccessorType& om_acc, size_t ix, size_t iy, size_t iz, int q, double* cfl_lin, NumFluxAccessorType& numFlux,
                    NumPtotalAccessorType& num_ptotal, const bool thermal, const double gamma, const double problem_cs2, const double denominator,
                    const double half_beta, const int fluidType, IdxConstType idx, FluidConstType fluidConst, OutputAccType& OutAcc) {

        //no vector in device code
        double uPri[gpu::FaceMax][N_OMINT] = {1.0};
        double uCon[gpu::FaceMax][N_OMINT] = {};
        double physFlux[gpu::FaceMax][N_OMINT] = {};
        double physValPtherm[gpu::FaceMax] = {};
        double physValPtotal[gpu::FaceMax] = {};

        double uPriOld[gpu::FaceMax][N_OMINT] = {};

        //if (ix == ixEnd && iy == gdata.mx[1]+n_ghost[1]-1 && iz == gdata.mx[2]+n_ghost[2]-1) cout << "reach limit: " << ix << "." <<  iy << "." << iz << "\n";
        //const int fluidType = Riemann[DirX]->get_Fluid_Type();

        //pointwise reconstruction (incomplete)
        gpu::compute(om_acc, uPri, ix, iy, iz, q, DirAll);

        for (int i = 0; i < gpu::FaceMax; i++) {
            for (int j = 0; j < N_OMINT; j++) {
                OutAcc[i][j] = uPri[i][j];
            }
        }
        
        if (ix >= 3 && iy >= 3 && iz >= 3) {

            gpu::compute(om_acc, uPriOld, ix -1, iy, iz, q, DirX);
            gpu::compute(om_acc, uPriOld, ix, iy -1, iz, q, DirY);
            gpu::compute(om_acc, uPriOld, ix, iy, iz -1, q, DirZ);

            for (int i = 0; i < N_OMINT; i++) {
                uPri[gpu::FaceEast][i] = uPri[gpu::FaceEast][i];
                uPri[gpu::FaceNorth][i] = uPri[gpu::FaceNorth][i];
                uPri[gpu::FaceTop][i] = uPri[gpu::FaceTop][i];
            }
        
        }
        
        for (int dir = 0; dir < DirMax; ++dir) {
            int ixOff = (dir == DirX) ? ix : ix - 1;
            int iyOff = (dir == DirY) ? iy : iy - 1;
            int izOff = (dir == DirZ) ? iz : iz - 1;

            int face = dir * 2;

            gpu::get_Cons(om_acc, uPri[face], uCon[face], physValPtherm[face], physValPtotal[face], ix, iy, iz, face,
                            thermal, gamma, problem_cs2, denominator, half_beta, fluidType, fluidConst);

            gpu::get_PhysFlux(om_acc, physFlux[face], uPri[face], uCon[face], physValPtherm[face], face, fluidConst);

            gpu::get_Cons(om_acc, uPri[face+1], uCon[face+1], physValPtherm[face+1], physValPtotal[face+1], ix, iy, iz, face+1,
                            thermal, gamma, problem_cs2, denominator, half_beta, fluidType, fluidConst);

            gpu::get_PhysFlux(om_acc, physFlux[face+1], uPri[face+1], uCon[face+1], physValPtherm[face+1], face+1, fluidConst);
        }

        double numVals_Ch[DirMax][gpu::NumV_Max] = {1., 1.};

        if (ix >= 2 && iy >= 2 && iz >= 2) {
            for (int dir = 0; dir < DirMax; ++dir) {
                int face = dir * 2;
                *cfl_lin = cl::sycl::fmax(gpu::get_vChar(om_acc, uPri[face], uCon[face], physValPtherm[face], uPri[face+1], uCon[face+1],
                                physValPtherm[face+1], numVals_Ch[dir], dir, gamma, idx, fluidConst), *cfl_lin);
                num_ptotal[dir] = gpu::get_NumFlux(om_acc, uPri[face+1], uCon[face+1], physFlux[face+1], physValPtherm[face+1],
                                        physValPtotal[face+1], uPri[face], uCon[face], physFlux[face], physValPtherm[face], 
                                        physValPtotal[face], numFlux[dir], numVals_Ch[dir], dir, gamma, fluidConst);
                
            }
        }
        //double testArray[gpu::FaceMax][N_OMINT] = {};
        //compute(om_acc, testArray, ix, iy, iz, q, DirAll);
        //*cfl_lin = OutAcc[1][3];
        /*for (int i = 0; i < gpu::FaceMax; i++) {
            for (int j = 0; j < N_OMINT; j++) {
                OutAcc[i][j] = testArray[i][j];
            }
        }*/
        *cfl_lin = OutAcc[1][3];
    }


}

#endif