#!/bin/bash

# Name of your job.
#SBATCH --job-name=Cronos_Block_Celerity_Test
#SBATCH --output=build/output.txt

# Send status information to this email address. 
#SBATCH --mail-user=daniel.proksch@student.uibk.ac.at

# Send an e-mail when the job has finished or failed. 
#SBATCH --mail-type=END,FAIL

#SBATCH --partition=nv75

#SBATCH --nodes=2
#SBATCH --ntasks-per-node=8

#SBATCH --time=0-00:30



cd build/
##nvidia-smi
##./proj tests/SedovExplosionHllcRef_tiny SedovCartesian_hllc
mpirun -n 4 proj tests/SedovExplosionHllcRef_small SedovCartesian_hllc
##ctest .
###nsys profile --capture-range cudaProfilerApi --capture-range-end stop --cuda-memory-usage=true --trace cuda,mpi --gpu-metrics-device=all proj tests/SedovExplosionHllcRef_medium SedovCartesian_hllc


###mpirun -n 4 nsys profile --cuda-memory-usage=true --trace cuda --gpu-metrics-device=all proj tests/SedovExplosionHllcRef_small SedovCartesian_hllc
