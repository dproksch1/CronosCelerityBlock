#ifndef CRONOS_COMPUTE_STEP_H
#define CRONOS_COMPUTE_STEP_H

#include "reconst_block.H"
#include "transformations.H"
#include "RiemannSolverHD.H"
#include "PhysFluxesHD.H"
#include "utils.H"

namespace gpu {

    template<typename uPriAccessorType, typename NumFluxAccessorType, typename NumPtotalAccessorType, typename IdxConstType, typename FluidConstType>
    void computeStep(const uPriAccessorType& uPri_acc, size_t ix, size_t iy, size_t iz, size_t ixyz, double* cfl_lin, NumFluxAccessorType& numFlux,
                    NumPtotalAccessorType& num_ptotal, const bool thermal, const double gamma, const double problem_cs2, const double denominator,
                    const double half_beta, const int fluidType, IdxConstType idx, FluidConstType fluidConst) {

        //no vector in device code
        double uPri[gpu::FaceMax][N_OMINT] = {};
        double uCon[gpu::FaceMax][N_OMINT] = {};
        double physFlux[gpu::FaceMax][N_OMINT] = {};
        double physValPtherm[gpu::FaceMax] = {};
        double physValPtotal[gpu::FaceMax] = {};

        for (int i = 0; i < gpu::FaceMax; i++) {
            for (int j = 0; j < N_OMINT; j++) {
                uPri[i][j] = uPri_acc[ixyz][i][j];
            }
        }
        
        for (int dir = 0; dir < DirMax; ++dir) {
            int ixOff = (dir == DirX) ? ix : ix - 1;
            int iyOff = (dir == DirY) ? iy : iy - 1;
            int izOff = (dir == DirZ) ? iz : iz - 1;

            int face = dir * 2;

            gpu::get_Cons(uPri[face], uCon[face], physValPtherm[face], physValPtotal[face], ix, iy, iz, face,
                            thermal, gamma, problem_cs2, denominator, half_beta, fluidType, fluidConst);

            gpu::get_PhysFlux(physFlux[face], uPri[face], uCon[face], physValPtherm[face], face, fluidConst);

            gpu::get_Cons(uPri[face+1], uCon[face+1], physValPtherm[face+1], physValPtotal[face+1], ix, iy, iz, face+1,
                            thermal, gamma, problem_cs2, denominator, half_beta, fluidType, fluidConst);

            gpu::get_PhysFlux(physFlux[face+1], uPri[face+1], uCon[face+1], physValPtherm[face+1], face+1, fluidConst);
        }

        double numVals_Ch[DirMax][gpu::NumV_Max] = {};

        if (ix >= 2 && iy >= 2 && iz >= 2) {
            for (int dir = 0; dir < DirMax; ++dir) {
                int face = dir * 2;
                *cfl_lin = cl::sycl::fmax(gpu::get_vChar(uPri[face], uCon[face], physValPtherm[face], uPri[face+1], uCon[face+1],
                                physValPtherm[face+1], numVals_Ch[dir], dir, gamma, idx, fluidConst), *cfl_lin);
                num_ptotal[dir] = gpu::get_NumFlux(uPri[face+1], uCon[face+1], physFlux[face+1], physValPtherm[face+1],
                                        physValPtotal[face+1], uPri[face], uCon[face], physFlux[face], physValPtherm[face], 
                                        physValPtotal[face], numFlux[dir], numVals_Ch[dir], dir, gamma, fluidConst);
            }
        }
    }

    template<typename uPriAccessorType, typename uPriAccessorType2, typename NumFluxAccessorType, typename NumPtotalAccessorType, typename IdxConstType, typename FluidConstType>
    void computeStep2(const uPriAccessorType& uPri_acc, uPriAccessorType2& uPri, size_t ix, size_t iy, size_t iz, size_t ixyz,double* cfl_lin, NumFluxAccessorType& numFlux,
                    NumPtotalAccessorType& num_ptotal, const bool thermal, const double gamma, const double problem_cs2, const double denominator,
                    const double half_beta, const int fluidType, IdxConstType idx, FluidConstType fluidConst) {

        //no vector in device code
        double uCon[gpu::FaceMax][N_OMINT] = {};
        double physFlux[gpu::FaceMax][N_OMINT] = {};
        double physValPtherm[gpu::FaceMax] = {};
        double physValPtotal[gpu::FaceMax] = {};

        double uPri2[gpu::FaceMax][N_OMINT] = {};
        for (int i = 0; i < gpu::FaceMax; i++) {
            for (int j = 0; j < N_OMINT; j++) {
                uPri2[i][j] = uPri_acc[ixyz][i][j];
            }
        }

        // for (int i = 0; i < gpu::FaceMax; i++) {
        //     for (int j = 0; j < N_OMINT; j++) {
        //         uPri[i][j] = uPri2[i][j];
        //     }
        // }
        
        for (int dir = 0; dir < DirMax; ++dir) {
            int ixOff = (dir == DirX) ? ix : ix - 1;
            int iyOff = (dir == DirY) ? iy : iy - 1;
            int izOff = (dir == DirZ) ? iz : iz - 1;

            int face = dir * 2;

            gpu::get_Cons(uPri[face], uCon[face], physValPtherm[face], physValPtotal[face], ix, iy, iz, face,
                            thermal, gamma, problem_cs2, denominator, half_beta, fluidType, fluidConst);

            gpu::get_PhysFlux(physFlux[face], uPri[face], uCon[face], physValPtherm[face], face, fluidConst);

            gpu::get_Cons(uPri[face+1], uCon[face+1], physValPtherm[face+1], physValPtotal[face+1], ix, iy, iz, face+1,
                            thermal, gamma, problem_cs2, denominator, half_beta, fluidType, fluidConst);

            gpu::get_PhysFlux(physFlux[face+1], uPri[face+1], uCon[face+1], physValPtherm[face+1], face+1, fluidConst);
        }

        double numVals_Ch[DirMax][gpu::NumV_Max] = {};

        if (ix >= 2 && iy >= 2 && iz >= 2) {
            for (int dir = 0; dir < DirMax; ++dir) {
                int face = dir * 2;
                *cfl_lin = cl::sycl::fmax(gpu::get_vChar(uPri[face], uCon[face], physValPtherm[face], uPri[face+1], uCon[face+1],
                                physValPtherm[face+1], numVals_Ch[dir], dir, gamma, idx, fluidConst), *cfl_lin);
                num_ptotal[dir] = gpu::get_NumFlux(uPri[face+1], uCon[face+1], physFlux[face+1], physValPtherm[face+1],
                                        physValPtotal[face+1], uPri[face], uCon[face], physFlux[face], physValPtherm[face], 
                                        physValPtotal[face], numFlux[dir], numVals_Ch[dir], dir, gamma, fluidConst);
            }
        }
    }


    template<typename uPriAccessorType, typename uPriAccessorType3, typename NumFluxAccessorType, typename NumPtotalAccessorType, typename IdxConstType, typename FluidConstType>
    void computeStep3(const uPriAccessorType& uPri_acc, uPriAccessorType3& uPri2, size_t ix, size_t iy, size_t iz, size_t ixyz,double* cfl_lin, NumFluxAccessorType& numFlux,
                    NumPtotalAccessorType& num_ptotal, const bool thermal, const double gamma, const double problem_cs2, const double denominator,
                    const double half_beta, const int fluidType, IdxConstType idx, FluidConstType fluidConst) {

        double uPri[gpu::FaceMax][N_OMINT] = {};
        for (int i = 0; i < gpu::FaceMax; i++) {
            for (int j = 0; j < N_OMINT; j++) {
                uPri[i][j] = uPri_acc[ixyz][i][j];
            }
        }

        computeStep2(uPri_acc, uPri, ix, iy, iz + 1, ixyz + 1, cfl_lin, numFlux, num_ptotal, thermal, gamma, problem_cs2,
													denominator, half_beta, fluidType, idx, fluidConst);
    }

    template<typename uPriAccessorType, typename omAccessorType, typename ghostType, typename rangeEndType>
    void reconstFunc(uPriAccessorType& uPri_acc, const omAccessorType& om_acc, size_t ix, size_t iy,
                     size_t iz, size_t ixyz, int q, const ghostType& n_ghost, const rangeEndType& rangeEnd) {

        if (ix >= n_ghost[0] && ix <= rangeEnd[0] && iy >= n_ghost[1] && iy <= rangeEnd[1] && iz >= n_ghost[2] && iz <= rangeEnd[2]) {
            if (ix >= n_ghost[0] && iy >= n_ghost[1] && iz >= n_ghost[2]) {
                auto [uPriWest, uPriSouth, uPriBottom] = gpu::computeWSB(om_acc, ix, iy, iz);
                auto [uPriEast, uPriNorth, uPriTop] = gpu::computeENT(om_acc, ix, iy, iz);
                uPri_acc[ixyz][gpu::FaceWest][q] = uPriWest;
                uPri_acc[ixyz][gpu::FaceEast][q] = uPriEast;
                uPri_acc[ixyz][gpu::FaceNorth][q] = uPriNorth;
                uPri_acc[ixyz][gpu::FaceSouth][q] = uPriSouth;
                uPri_acc[ixyz][gpu::FaceBottom][q] = uPriBottom;
                uPri_acc[ixyz][gpu::FaceTop][q] = uPriTop;
            }
        }
    }

}

#endif