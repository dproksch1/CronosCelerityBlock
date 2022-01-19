# CronosCode Block Proxyapp

#### Standalone usage of the proxyapp:
1. Install the dependencies CronosNumLib, GSL, MPI, HDF5, Git and either hipSYCL or ComputeCpp
2. Build the CronosCode Block Project by executing (from the directory containing this file):
            `cmake -B [PATH_TO_BUILD_DIR] -DCMAKE_BUILD_TYPE=Release`
3. If one of the dependencies isn't found on its own, add its path using `-DCMAKE_PREFIX_PATH`
4. After the build succeeds, compile the project by executing:
            `cmake -build [PATH_TO_BUILD_DIR] --config Release`
5. For testing the project you can execute `ctest .` from inside the build directory


#### Usage as a submodule of Cronos:
1. Add to Cronos using Git Submodule
2. Add directories to include list in CMakeList.txt:
    `cronos_block/configuration/standalone/cronos`
    `cronos_block/kernel`
    `cronos_block/reconstruction`
    `cronos_block/riemann_solver`
3. Add all `.C`-files in directories to SOURCES in CMakeList.txt:
    `cronos_block/kernel`
    `cronos_block/reconstruction`