export CUDA_HOME=/usr/local/cuda
#export CXX="clang++-9"
#export CC="clang-9"
export CXX="clang++-10"
export CC="clang-10"

cmake -G Ninja .. -DCMAKE_PREFIX_PATH="/home/philipp/software-local/gsl-2.6;/home/philipp/software-local/hdf5-1.12.0;/lib-installs-local/hipsycl/lib" -DHIPSYCL_PLATFORM=cuda -DHIPSYCL_GPU_ARCH=sm_75 -DCMAKE_BUILD_TYPE=Debug
