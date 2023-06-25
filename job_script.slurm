#!/bin/bash

# Name of your job.
#SBATCH --job-name=Cronos_Block_Celerity_Test

# Send status information to this email address. 
#SBATCH --mail-user=daniel.proksch@student.uibk.ac.at

# Send an e-mail when the job has finished or failed. 
#SBATCH --mail-type=END,FAIL

#SBATCH --time=0-00:30

#SBATCH --gpus=4

cmake -B build -DHIPSYCL_TARGETS="cuda:sm_86" -DDIRECTLINK_CNL=true -DCMAKE_BUILD_TYPE=Release
cmake --build ./build --config Release
cd build/
./proj tests/SedovExplosionHllcRef_small SedovCartesian_hllc


