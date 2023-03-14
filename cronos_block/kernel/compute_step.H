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

}

#endif