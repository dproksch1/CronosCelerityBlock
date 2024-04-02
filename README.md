# CronosCelerityBlock (CronosCelerity Proxyapp)

This is the main implementation proxyapp of [CronosCelerity](https://github.com/philippgs/CronosCelerity), a partial port of the astrophysics simulation library **Cronos** to the **Celerity** framework. This port allows for an execution of simulations of single fluid hydrodynamics problems on GPUs and distributed-memory accelerator clustors.

CronosCelerityBlock is usable both as a part of the CronosCelerity codebase and as a standalone application. There are no major differences between these execution types as the only goal of the bigger CronosCelerity codebase is to provide a better start-off point for future port extensions, to cover more of Cronos featureset. The proxyapp on the other hand is considerally stripped down, as it is fully focussed on the scope of our current port.

## Cronos

[Cronos](https://www.uibk.ac.at/astro/research_groups/ralf-kissmann/home/cronos-code.html.en) is a project by Kissmann et al. at the University of Innsbruck. The goal of the project is to provide an implementation for the solution of partial-differential equations representing hyperbolic conversion law, that form the basis of hydrodynamical and magnetohydrodynamical computations. Using these calculations, Cronos is able to perform simluations for a variety of HD and MHD problem types.

## Celerity

[Celerity](https://celerity.github.io/) is a programming model and runtime environment developed as a research project by the *Distributed and Parallel Systems Group* at the University of Innsbruck. The goal of the project is to provide an intuitive, high-level framework for the creation of distributed memory parallel programs. To achieve this, Celerity extends the royalty-free heterogeneous device programming standard [SYCL](https://www.khronos.org/sycl/), defined by the *Khronos Group*, by adding an MPI-based multi-process distributed-memory execution layer and an intelligent workload and data dependency management system, that allows for the automatic asynchrounous and parallel scheduling of tasks and data-transfers. Celerity is not a standalone SYCL implementation, but instead utilizes another installed implementation as its basis, extending and modifying its functionality.

## Utilization of CronosCelerity

Refer to the README in the [CronosCelerity repository](https://github.com/philippgs/CronosCelerity).

## Utilization of CronosCelerityBlock as a Standalone Application

### Dependencies

CronosCelerityBlock was implemented using Celerity 0.3.2 and a modified version of hipSYCL 0.9.3 (since then renamed to AdaptiveCpp) that enables the utilization of SYCL 2020 reductions. We recommend using these versions, as the newer versions of Celerity differ in syntax and do not support this method of achieving reductions anymore. We did not test the application for DPC++.

#### Full dependency list:
- [Git](https://git-scm.com/)
- [GNU Scientific Library (GSL)](https://www.gnu.org/software/gsl/) &nbsp; *we recommend version 2.5 or newer*
- Message Passing Interface (MPI): [MPICH](https://www.mpich.org/), [Open MPI](https://www.open-mpi.org/)
- [Parallel HDF5](https://www.hdfgroup.org/solutions/hdf5/) *we recommend version 1.10.4 or newer*
- [hipSYCL 0.9.3 (now called AdaptiveCpp)](https://github.com/AdaptiveCpp/AdaptiveCpp/releases/tag/v0.9.3) and the [SYCL 2020 reduction patch](https://github.com/AdaptiveCpp/AdaptiveCpp/pull/578) by Fabian Knorr &nbsp; *there is a [repository fork](https://github.com/fknorr/hipSYCL) with these changes, but we did not test it*
    - [LLVM with Clang](https://llvm.org/) &nbsp; *for hipSYCL; we recommend version 15.0.4 or newer*
    - CUDA, OpenCL or rocM &nbsp; *for hipSYCL; we tested our version using CUDA and exclusively on NVIDIA GPUs*
- [Celerity 0.3.2](https://github.com/celerity/celerity-runtime/releases/tag/v0.3.2) &nbsp; *we recommend cloning the [repository](https://github.com/celerity/celerity-runtime) and fetching commit [59e7d61](https://github.com/celerity/celerity-runtime/tree/59e7d61d3e5ca96fe6cdc38df42ea8b16f76762c) or an adjacent commit*

If you want to use a newer version of Celerity, you have to update the syntax and parts of the implementation structure according to the changes made. Additionally, you have to swap the underlying SYCL implementation to [DPC++](https://www.intel.com/content/www/us/en/developer/tools/oneapi/dpc-compiler.html#gs.76awfe) or another implementation, that supports SYCL 2020 reductions. Refer to the [Celerity](https://celerity.github.io/docs/getting-started) and [DPC++](https://www.intel.com/content/www/us/en/developer/tools/oneapi/dpc-compiler-documentation.html) documentations for any attempted changes.

Due the reliance of hipSYCL (AdaptiveCpp) on LLVM/Clang, we do not support windows for the execution of our implementation out of the box. For local executions outside of clusters we recommend Linux or [WSL2](https://learn.microsoft.com/en-us/windows/wsl/install).

### Preperation and Initial Compilation

Fetch the git submodule using:
```
git submodule update --init --recursive
```

Setup all depencendies for your system and then setup the build directory of CronosCelerityBlock using:
```
cmake -B build -DHIPSYCL_TARGETS="<targets>" -DDIRECTLINK_CNL=true -DHDF5_PREFER_PARALLEL=true -DCMAKE_BUILD_TYPE=Release
```
To find your respective compilation targets (used for the setups of SYCL and Celerity as well as the compilation of Celerity applications), refer to the [AdaptiveCPP documentation](https://github.com/AdaptiveCpp/AdaptiveCpp/blob/develop/doc/using-hipsycl.md) and the [CUDA compatibility](https://developer.nvidia.com/cuda-gpus) or [AMDGPU processor](https://llvm.org/docs/AMDGPUUsage.html#processors) overview. If you want to use another SYCL version, refer to the documentation of that implementation.

`-DDIRECTLINK_CNL=true` activates the directly linkage of the CronosNumLib version provided in the dependencies folder. As an alterative you can install CronosNumLib directly on your system.

Whether `-DHDF5_PREFER_PARALLEL=true` is needed can depend on your installed HDF5 version.

If one of the dependencies isn't found on its own, add its path using `-DCMAKE_PREFIX_PATH`.

Then perform the compilation of CronosCelerityBlock using
```
cmake --build ./build --config Release
```
You can execute `ctest .` from inside of the build directory to test the application.

### Simulation Setup

To setup your own simulation in CronosCelerityBlock, add a new folder to the `cronos/configuration` directory containing:

- a `modules.C` file that relates simulation classes to their types (see Section 3.2 in the [Cronos Manual](external/CronosDocs/Cronos-Manual.pdf)).
- a `mod_<SimulationName>.H` file that defines a class for your simulation, extending the class `ProblemType` and implementing the problem's constructor and user functions. For reference see the mod-files for the Sod shock tube or Sedov explosion simulation setups.
- a `constants.H` that defines certain parameters and global variables in the program. The file still follows the general structure layed out in Section 3.2 of the [Cronos Manual](external/CronosDocs/Cronos-Manual.pdf), but with the addition of `CRONOS_DISTR_OUTPUT` that sets the Celerity output type as eather distributed (`CRONOS_ON`) or performed on master-node (`CRONOS_OFF`). For reference see the constants-files for the Sod shock tube or Sedov explosion simulation setups.

The folder does **not** need to contain a `Makefile` or `machinefile.local` as both of these aspects are managed by cmake and Celerity respectively.

The next important step is to add the configuration folder for your simulation to the linker via `include_directories` in the file `CMakeLists.txt`. `modules.C` and all other non-header C++ files additionally have to be added to the `set(SOURCES)` command in the same file. Remove or comment out all configurations files and folders from these two sections in `CMakeLists.txt` to prevent linking issues and redefinitions.

The final step of the simulation setup is the creation of a cat file `<pname>.cat`, that handles the simulations input parameters. As reference see Section 3.10 of the [Cronos Manual](external/CronosDocs/Cronos-Manual.pdf) and the premade cat files for the Sod shock tube or Sedov explosion simulation setups, found in the directory `tests`. The source directory of our own cat file does not matter, as it will be provided to the application via a command line argument.

### Usage

Recompile CronosCelerityBlock using the commands from the ***Preperation and Initial Compilation*** Section of this README. The execute the program via
```
build/proj <datadir> <pname>
```
where *pname* is your chosen title of your cat file and *datadir* is both the relative path to the directory containing the cat file (so that *datadir/pname.cat* is the relative path to your cat file). *datadir* is also reused by the application as an output folder.

## Context

This project was created as part of a master's thesis by Daniel Proksch at the University of Innbruck. The thesis and project were supervised by Phillip Gschwandtner.

