/**
 * @file
 * @brief The GPU-compatible compute step for the main kernel
 *
 * This file provides the compute step for the main kernel defined and executed in HyperbolicSolver::singlestep(Data&, GridFunc&,ProblemType&, int, Queue&).
 * The function is executed 4 times per kernel iteration and utilizes additional functions defined in solver_func.H and reconstruction_block.H.
 */

#ifndef CRONOS_COMPUTE_STEP_H
#define CRONOS_COMPUTE_STEP_H

#include "reconst_block.H"
#include "transformations_block.H"
#include "solver_func.H"
#include "utils.H"

namespace gpu {

    //! @brief Reconstructs the vertices and computes the numerical fluxes for the given cell
    template<typename omAccessorType, typename NumFluxAccessorType, typename NumPtotalAccessorType, typename CFlagAccessorType, typename IdxConstType, typename FluidConstType>
    void computeStep(const omAccessorType& om_rho_acc, const omAccessorType& om_sx_acc, const omAccessorType& om_sy_acc, const omAccessorType& om_sz_acc, const omAccessorType& om_Eges_acc,
                    size_t ix, size_t iy, size_t iz, double* cfl_lin, NumFluxAccessorType& numFlux,
                    NumPtotalAccessorType& num_ptotal, const CFlagAccessorType& carbuncle_flag, const bool thermal, const double gamma, const double problem_cs2, const double denominator,
                    const double half_beta, const int fluidType, const double Theta, bool use_carbuncle, IdxConstType idx, FluidConstType fluidConst) {

        double uPri[gpu::FaceMax][N_OMINT] = {};
        double uCon[gpu::FaceMax][N_OMINT] = {};
        double physFlux[gpu::FaceMax][N_OMINT] = {};
        double physValPtherm[gpu::FaceMax] = {};
        double physValPtotal[gpu::FaceMax] = {};
        int cflag[gpu::FaceMax] = {};

        {
            auto [uPriWest0, uPriSouth0, uPriBottom0] = gpu::computeWSB(om_rho_acc, Theta, ix, iy, iz);
            auto [uPriEast0, uPriNorth0, uPriTop0] = gpu::computeENT(om_rho_acc, Theta, ix, iy, iz);
            uPri[gpu::FaceWest][0] = uPriWest0;
            uPri[gpu::FaceEast][0] = uPriEast0;
            uPri[gpu::FaceSouth][0] = uPriSouth0;
            uPri[gpu::FaceNorth][0] = uPriNorth0;
            uPri[gpu::FaceBottom][0] = uPriBottom0;
            uPri[gpu::FaceTop][0] = uPriTop0;

            auto [uPriWest1, uPriSouth1, uPriBottom1] = gpu::computeWSB(om_sx_acc, Theta, ix, iy, iz);
            auto [uPriEast1, uPriNorth1, uPriTop1] = gpu::computeENT(om_sx_acc, Theta, ix, iy, iz);
            uPri[gpu::FaceWest][1] = uPriWest1;
            uPri[gpu::FaceEast][1]= uPriEast1;
            uPri[gpu::FaceNorth][1] = uPriNorth1;
            uPri[gpu::FaceSouth][1] = uPriSouth1;
            uPri[gpu::FaceBottom][1] = uPriBottom1;
            uPri[gpu::FaceTop][1] = uPriTop1;

            auto [uPriWest2, uPriSouth2, uPriBottom2] = gpu::computeWSB(om_sy_acc, Theta, ix, iy, iz);
            auto [uPriEast2, uPriNorth2, uPriTop2] = gpu::computeENT(om_sy_acc, Theta, ix, iy, iz);
            uPri[gpu::FaceWest][2] = uPriWest2;
            uPri[gpu::FaceEast][2] = uPriEast2;
            uPri[gpu::FaceNorth][2] = uPriNorth2;
            uPri[gpu::FaceSouth][2] = uPriSouth2;
            uPri[gpu::FaceBottom][2] = uPriBottom2;
            uPri[gpu::FaceTop][2] = uPriTop2;

            auto [uPriWest3, uPriSouth3, uPriBottom3] = gpu::computeWSB(om_sz_acc, Theta, ix, iy, iz);
            auto [uPriEast3, uPriNorth3, uPriTop3] = gpu::computeENT(om_sz_acc, Theta, ix, iy, iz);
            uPri[gpu::FaceWest][3] = uPriWest3;
            uPri[gpu::FaceEast][3] = uPriEast3;
            uPri[gpu::FaceNorth][3] = uPriNorth3;
            uPri[gpu::FaceSouth][3] = uPriSouth3;
            uPri[gpu::FaceBottom][3] = uPriBottom3;
            uPri[gpu::FaceTop][3] = uPriTop3;

            auto [uPriWest4, uPriSouth4, uPriBottom4] = gpu::computeWSB(om_Eges_acc, Theta, ix, iy, iz);
            auto [uPriEast4, uPriNorth4, uPriTop4] = gpu::computeENT(om_Eges_acc, Theta, ix, iy, iz);
            uPri[gpu::FaceWest][4] = uPriWest4;
            uPri[gpu::FaceEast][4] = uPriEast4;
            uPri[gpu::FaceNorth][4] = uPriNorth4;
            uPri[gpu::FaceSouth][4] = uPriSouth4;
            uPri[gpu::FaceBottom][4] = uPriBottom4;
            uPri[gpu::FaceTop][4] = uPriTop4;
        }
        
        for (int dir = 0; dir < DirMax; ++dir) {

            int face = dir * 2;

            gpu::get_Cons(uPri[face], uCon[face], physValPtherm[face], physValPtotal[face], face,
                            thermal, gamma, problem_cs2, denominator, half_beta, fluidType, fluidConst);

            gpu::get_PhysFlux(physFlux[face], uPri[face], uCon[face], physValPtherm[face], face, fluidConst);

            gpu::get_Cons(uPri[face+1], uCon[face+1], physValPtherm[face+1], physValPtotal[face+1], face+1,
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
                                        physValPtotal[face+1], cflag[face+1], uPri[face], uCon[face], physFlux[face], physValPtherm[face], 
                                        physValPtotal[face], cflag[face], numFlux[dir], numVals_Ch[dir], use_carbuncle, dir, gamma, fluidConst);
            }
        }
    }
}

#endif