# CronosCode Block Proxyapp

#### Standalone usage of the proxyapp:
1. Install the dependencies GSL, MPI, HDF5, Git and either hipSYCL or ComputeCpp
2. Either install CronosNumLib directly or decide to use the prebuild version in the dependency folder
3. Build the CronosCode Block Project by executing (from the directory containing this file):
            `cmake -B [PATH_TO_BUILD_DIR] -DCMAKE_BUILD_TYPE=Release`
    1. If one of the dependencies isn't found on its own, add its path using `-DCMAKE_PREFIX_PATH`
    2. If CronosNumLib is not installed, add `-DDIRECTLINK_CNL=true` to enable direct linkage
    3. If your MPI-HDF5 integration results in a compilation error, try adding `-DHDF5_PREFER_PARALLEL=true`
5. After the build succeeds, compile the project by executing:
            `cmake -build [PATH_TO_BUILD_DIR] --config Release`
6. For testing the project you can execute `ctest .` from inside the build directory
7. In regards to simulation setup, refer back to the Cronos manual


#### Usage as a submodule of Cronos:
1. Add to Cronos using Git Submodule
2. Add directories to include list in CMakeList.txt:
    `cronos_block/configuration/standalone/cronos`
    `cronos_block/data`
    `cronos_block/stepfunctions`
    `cronos_block/kernel`
    `cronos_block/reconstruction`
    `cronos_block/utils`
    `cronos_block/queue`
3. Add all `.C`-files in directories to SOURCES in CMakeList.txt:
    `cronos_block/data`
    `cronos_block/kernel`
    `cronos_block/stepfunctions`