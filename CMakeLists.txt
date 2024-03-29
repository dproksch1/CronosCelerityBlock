cmake_minimum_required (VERSION 3.16.3)
project (CronosCodeBlock C CXX)
enable_language(Fortran)

set(CMAKE_CXX_STANDARD 17)
#set(CMAKE_VERBOSE_MAKEFILE ON)

##################
# Compiler Setup #
##################

if(MSVC)

	# disable "register no longer supported storage class" warning
	add_compile_options("/wd5033")
	# disable "sprintf is unsafe, use sprintf_s" warning
	add_compile_options("/wd4996")
	# disable "possible loss of data in numeric conversion" warning
	add_compile_options("/wd4244")
	# enable parallel build in Visual Studio
	add_compile_options("/MP")
	# required for some symbols to be available, such as SZ_encoder_enabled
	add_compile_options("-DH5_BUILT_AS_DYNAMIC_LIB")
	# enable debug information in object files
	add_compile_options("$<$<CONFIG:DEBUG>:/Z7;/DEBUG>")
	# set iterator checking
	add_compile_options("$<$<CONFIG:DEBUG>:-D_ITERATOR_DEBUG_LEVEL=2>")
	add_compile_options("$<$<CONFIG:RELEASE>:-D_ITERATOR_DEBUG_LEVEL=0>")
	# enable debug information in object files
	#set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} /MD")

else()

  add_compile_options("$<$<CONFIG:RELEASE>:-D_ITERATOR_DEBUG_LEVEL=0>")
endif()

############################
# Third-party Dependencies #
############################

list(APPEND CMAKE_MODULE_PATH ${CMAKE_SOURCE_DIR}/cmake)
include(ExternalProject)
find_package(GSL REQUIRED)
find_package(MPI REQUIRED)
find_package(HDF5 REQUIRED COMPONENTS C) # CXX HL

if(NOT DEFINED DIRECTLINK_CNL)
  SET(DIRECTLINK_CNL false)
endif()

if(NOT DIRECTLINK_CNL)
  find_package(CronosNumLib REQUIRED)
endif()

# SYCL-specific dependency setup

option(CELERITY_SYCL_IMPL "hipSYCL|ComputeCpp")

# Find SYCL implementation.
if(CELERITY_SYCL_IMPL STREQUAL "hipSYCL")
  find_package(hipSYCL CONFIG REQUIRED)
elseif(CELERITY_SYCL_IMPL STREQUAL "ComputeCpp")
  find_package(ComputeCpp REQUIRED)
else()
  # We first check if hipSYCL can be found, otherwise we try ComputeCpp.
  find_package(hipSYCL CONFIG)
  if(hipSYCL_FOUND)
    set(CELERITY_SYCL_IMPL "hipSYCL")
  else()
    find_package(ComputeCpp REQUIRED)
    set(CELERITY_SYCL_IMPL "ComputeCpp")
  endif()
endif()

if(CELERITY_SYCL_IMPL STREQUAL "hipSYCL")
  message(STATUS "Found hipSYCL: ${hipSYCL_DIR}")
else()
  message(STATUS "Found ComputeCpp: ${ComputeCpp_DIR}")
endif()

if(CELERITY_SYCL_IMPL STREQUAL "ComputeCpp" AND NOT COMPUTECPP_USER_FLAGS MATCHES "-D_CRT_SECURE_NO_WARNINGS")
  set(COMPUTECPP_USER_FLAGS "${COMPUTECPP_USER_FLAGS};-D_CRT_SECURE_NO_WARNINGS=1"
    CACHE STRING "" FORCE)
endif()

find_package(Celerity CONFIG REQUIRED)

if(CMAKE_GENERATOR STREQUAL "Ninja")
  # Force colored warnings in Ninja's output, if the compiler has -fdiagnostics-color support.
  # Rationale in https://github.com/ninja-build/ninja/issues/814
  set(CMAKE_SYCL_FLAGS "${CMAKE_SYCL_FLAGS} -fdiagnostics-color=always")
endif()


##############################
# Application-specific Setup #
##############################

include_directories(${HDF5_INCLUDE_DIRS})
include_directories(${GSL_INCLUDE_DIRS})


if(DIRECTLINK_CNL)
  add_subdirectory(external/CronosNumLib)
else()
  include_directories(${CRONOSNUMLIB_INCLUDE_DIR}/CronosNumLib)
endif()

# include_directories(cronos_block/test)
include_directories(cronos_block/data)
include_directories(cronos_block/general)
include_directories(cronos_block/interface_functions)
include_directories(cronos_block/kernel)
include_directories(cronos_block/queue)
include_directories(cronos_block/reconstruction)
include_directories(cronos_block/riemann_solver)
include_directories(cronos_block/stepfunctions)
# include_directories(cronos_block/utils)
include_directories(cronos_block/configuration/standalone/cronos_block)
# include_directories(cronos_block/configuration/ShockTubeSod)
include_directories(cronos_block/configuration/SedovExplosion)
include_directories(cronos_block/configuration/SedovExplosion/randgen)

if(DIRECTLINK_CNL)
  message(STATUS "Direct Inclusion of CronosNumLib")
  include_directories(external/CronosNumLib/Matrix)
  include_directories(external/CronosNumLib/util)
endif()

set(SOURCES
    cronos_block/general/buildinfo.C
    cronos_block/general/CException.C
    cronos_block/general/cfl.C
    cronos_block/general/fluid.C
    cronos_block/general/grid.C
    cronos_block/general/gridParam.C
    cronos_block/general/movie.C
    cronos_block/general/normalisation.C
    cronos_block/general/proj.C
    cronos_block/general/specific.C
    cronos_block/general/units.C
    cronos_block/data/data.C
    cronos_block/interface_functions/problem.C
    cronos_block/interface_functions/solver.C
    cronos_block/kernel/fields_local.C
    cronos_block/kernel/singlestep_block.C
    cronos_block/kernel/transformations_block.C
    cronos_block/kernel/solver_func.C
    cronos_block/riemann_solver/RiemannSolver.C
    cronos_block/stepfunctions/gridfunc.C
    cronos_block/stepfunctions/gridgen.C
    cronos_block/stepfunctions/Hdf5File_cbase.C
    cronos_block/stepfunctions/timestepping.C
    #test/proj.C
    #test/test_data.C
    #test/singlestep_stub.C
    # cronos_block/configuration/ShockTubeSod/modules.C
    cronos_block/configuration/SedovExplosion/modules.C
    cronos_block/configuration/SedovExplosion/randgen/randgen.C

    #cronos_block/utils/celerity_test.C
)

FILE(GLOB_RECURSE HEADERS "*.h")

add_executable(proj ${SOURCES} cronos_block/general/proj.C ${HEADERS})
#add_executable(proj ${SOURCES} test/proj.C ${HEADERS})

set_property(TARGET proj PROPERTY CXX_STANDARD 17)
#target_link_options (proj PRIVATE "$<$<BOOL:$<TARGET_PROPERTY:POSITION_INDEPENDENT_CODE>>:-fPIC>")

# properly set up SYCL compiler in case of hipSYCL

set(DEVICE_SOURCES "")
if(CELERITY_SYCL_IMPL STREQUAL "hipSYCL")
  # For hipSYCL we have to pass all source files into add_sycl_to_target
  # Don't just do it in general (yields better ComputeCpp build performance)
  set(DEVICE_SOURCES ${SOURCES} cronos_block/general/proj.C)
endif()
add_sycl_to_target(TARGET proj SOURCES ${DEVICE_SOURCES})

add_celerity_to_target(TARGET proj SOURCES cronos_block/general/proj.C ${HEADERS})

if(NOT MSVC)
	target_link_libraries(proj PRIVATE stdc++fs)
endif()

if(DIRECTLINK_CNL)
  target_link_libraries(proj PRIVATE matrix util)
else()
  target_link_libraries(proj PUBLIC ${CRONOSNUMLIB_MATRIX} ${CRONOSNUMLIB_UTIL})
endif()
target_link_libraries(proj PRIVATE ${MPI_C_LIBRARIES})
target_link_libraries(proj PRIVATE GSL::gsl GSL::gslcblas)
target_link_libraries(proj PRIVATE ${HDF5_LIBRARIES})

#########
# TESTS #
#########

enable_testing()

# configure_file(${CMAKE_SOURCE_DIR}/tests/SodToro5HllcRef_tiny/SodToro5HllcRef.cat ${CMAKE_BINARY_DIR}/tests/SodToro5HllcRef_tiny/SodToro5HllcRef.cat COPYONLY)
# configure_file(${CMAKE_SOURCE_DIR}/tests/SodToro5HllcRef_small/SodToro5HllcRef.cat ${CMAKE_BINARY_DIR}/tests/SodToro5HllcRef_small/SodToro5HllcRef.cat COPYONLY)
# configure_file(${CMAKE_SOURCE_DIR}/tests/SodToro5HllcRef_medium/SodToro5HllcRef.cat ${CMAKE_BINARY_DIR}/tests/SodToro5HllcRef_medium/SodToro5HllcRef.cat COPYONLY)

# add_test(NAME tiny_run COMMAND $<TARGET_FILE:proj> ${CMAKE_BINARY_DIR}/tests/SodToro5HllcRef_tiny SodToro5HllcRef)
# add_test(NAME tiny_compare COMMAND ${HDF5_DIFF_EXECUTABLE} -v -d 1.5e-11 --exclude-path /Data/version ${CMAKE_BINARY_DIR}/tests/SodToro5HllcRef_tiny/SodToro5HllcRef_float/SodToro5HllcRef_flt_step3.h5 ${CMAKE_SOURCE_DIR}/tests/SodToro5HllcRef_tiny/SodToro5HllcRef_float/SodToro5HllcRef_flt_step3.h5)
	
# add_test(NAME small_run COMMAND $<TARGET_FILE:proj> ${CMAKE_BINARY_DIR}/tests/SodToro5HllcRef_small SodToro5HllcRef)
# add_test(NAME small_compare COMMAND ${HDF5_DIFF_EXECUTABLE} -v -d 1.5e-11 --exclude-path /Data/version ${CMAKE_BINARY_DIR}/tests/SodToro5HllcRef_small/SodToro5HllcRef_float/SodToro5HllcRef_flt_step16.h5 ${CMAKE_SOURCE_DIR}/tests/SodToro5HllcRef_small/SodToro5HllcRef_float/SodToro5HllcRef_flt_step16.h5)

# add_test(NAME medium_run COMMAND $<TARGET_FILE:proj> ${CMAKE_BINARY_DIR}/tests/SodToro5HllcRef_medium SodToro5HllcRef)
# add_test(NAME medium_compare COMMAND ${HDF5_DIFF_EXECUTABLE} -v -d 1.5e-11 --exclude-path /Data/version ${CMAKE_BINARY_DIR}/tests/SodToro5HllcRef_medium/SodToro5HllcRef_float/SodToro5HllcRef_flt_step990.h5 ${CMAKE_SOURCE_DIR}/tests/SodToro5HllcRef_medium/SodToro5HllcRef_float/SodToro5HllcRef_flt_step990.h5)

# set_tests_properties(tiny_run PROPERTIES DEPENDS proj)
# set_tests_properties(tiny_run PROPERTIES FIXTURES_SETUP EXEC_TINY)
# set_tests_properties(tiny_compare PROPERTIES FIXTURES_REQUIRED EXEC_TINY)
	
# set_tests_properties(small_run PROPERTIES DEPENDS proj)
# set_tests_properties(small_run PROPERTIES FIXTURES_SETUP EXEC_SMALL)
# set_tests_properties(small_compare PROPERTIES FIXTURES_REQUIRED EXEC_SMALL)
	
# set_tests_properties(medium_run PROPERTIES DEPENDS proj)
# set_tests_properties(medium_run PROPERTIES FIXTURES_SETUP EXEC_MEDIUM)
# set_tests_properties(medium_compare PROPERTIES FIXTURES_REQUIRED EXEC_MEDIUM)

configure_file(${CMAKE_SOURCE_DIR}/tests/SedovExplosionHllcRef_tiny/SedovCartesian_hllc.cat ${CMAKE_BINARY_DIR}/tests/SedovExplosionHllcRef_tiny/SedovCartesian_hllc.cat COPYONLY)
configure_file(${CMAKE_SOURCE_DIR}/tests/SedovExplosionHllcRef_small/SedovCartesian_hllc.cat ${CMAKE_BINARY_DIR}/tests/SedovExplosionHllcRef_small/SedovCartesian_hllc.cat COPYONLY)
# configure_file(${CMAKE_SOURCE_DIR}/tests/SedovExplosionHllcRef_medium/SedovCartesian_hllc.cat ${CMAKE_BINARY_DIR}/tests/SedovExplosionHllcRef_medium/SedovCartesian_hllc.cat COPYONLY)

# add_test(NAME tiny_run COMMAND mpirun -n 2 $<TARGET_FILE:proj> ${CMAKE_BINARY_DIR}/tests/SedovExplosionHllcRef_tiny SedovCartesian_hllc)
# add_test(NAME tiny_compare COMMAND ${HDF5_DIFF_EXECUTABLE} -v -d 1.5e-11 --exclude-path /Data/version ${CMAKE_BINARY_DIR}/tests/SedovExplosionHllcRef_tiny/SedovCartesian_hllc_float/SedovCartesian_hllc_flt_step1.h5 ${CMAKE_SOURCE_DIR}/tests/SedovExplosionHllcRef_tiny/SedovCartesian_hllc_float/SedovCartesian_hllc_flt_step1.h5)
add_test(NAME small1_run COMMAND mpirun -n 1 $<TARGET_FILE:proj> ${CMAKE_BINARY_DIR}/tests/SedovExplosionHllcRef_small SedovCartesian_hllc)
add_test(NAME small1_compare COMMAND ${HDF5_DIFF_EXECUTABLE} --verbose=0 -d 1.5e-11 --exclude-path /Data/version ${CMAKE_BINARY_DIR}/tests/SedovExplosionHllcRef_small/SedovCartesian_hllc_float/SedovCartesian_hllc_flt_step4.h5 ${CMAKE_SOURCE_DIR}/tests/SedovExplosionHllcRef_small/SedovCartesian_hllc_float/SedovCartesian_hllc_flt_step4.h5)

add_test(NAME small2_run COMMAND mpirun -n 2 $<TARGET_FILE:proj> ${CMAKE_BINARY_DIR}/tests/SedovExplosionHllcRef_small SedovCartesian_hllc)
add_test(NAME small2_compare COMMAND ${HDF5_DIFF_EXECUTABLE} --verbose=0 -d 1.5e-11 --exclude-path /Data/version ${CMAKE_BINARY_DIR}/tests/SedovExplosionHllcRef_small/SedovCartesian_hllc_float/SedovCartesian_hllc_flt_step4.h5 ${CMAKE_SOURCE_DIR}/tests/SedovExplosionHllcRef_small/SedovCartesian_hllc_float/SedovCartesian_hllc_flt_step4.h5)

add_test(NAME small4_run COMMAND mpirun -n 4 $<TARGET_FILE:proj> ${CMAKE_BINARY_DIR}/tests/SedovExplosionHllcRef_small SedovCartesian_hllc)
add_test(NAME small4_compare COMMAND ${HDF5_DIFF_EXECUTABLE} --verbose=0 -d 1.5e-11 --exclude-path /Data/version ${CMAKE_BINARY_DIR}/tests/SedovExplosionHllcRef_small/SedovCartesian_hllc_float/SedovCartesian_hllc_flt_step4.h5 ${CMAKE_SOURCE_DIR}/tests/SedovExplosionHllcRef_small/SedovCartesian_hllc_float/SedovCartesian_hllc_flt_step4.h5)

add_test(NAME small8_run COMMAND mpirun -n 8 $<TARGET_FILE:proj> ${CMAKE_BINARY_DIR}/tests/SedovExplosionHllcRef_small SedovCartesian_hllc)
add_test(NAME small8_compare COMMAND ${HDF5_DIFF_EXECUTABLE} --verbose=0 -d 1.5e-11 --exclude-path /Data/version ${CMAKE_BINARY_DIR}/tests/SedovExplosionHllcRef_small/SedovCartesian_hllc_float/SedovCartesian_hllc_flt_step4.h5 ${CMAKE_SOURCE_DIR}/tests/SedovExplosionHllcRef_small/SedovCartesian_hllc_float/SedovCartesian_hllc_flt_step4.h5)
# add_test(NAME small_compare COMMAND ${HDF5_DIFF_EXECUTABLE} -v -d 1.5e-11 --exclude-path /Data/version ${CMAKE_SOURCE_DIR}/tests/SedovExplosionHllcRef_small/SedovCartesian_hllc_float/SedovCartesian_hllc_flt_step4.h5 ${CMAKE_BINARY_DIR}/tests/SedovExplosionHllcRef_small/SedovCartesian_hllc_float/SedovCartesian_hllc_flt_step4.h5)
# add_test(NAME small_compare COMMAND ${HDF5_DIFF_EXECUTABLE} -v -d 1.5e-11 --exclude-path /Data/version ${CMAKE_BINARY_DIR}/tests/SedovExplosionHllcRef_small/SedovCartesian_hllc_double/SedovCartesian_hllc_dbl_step1.h5 ${CMAKE_SOURCE_DIR}/tests/SedovExplosionHllcRef_small/SedovCartesian_hllc_double/SedovCartesian_hllc_dbl_step1.h5)
# add_test(NAME medium_run COMMAND mpirun -n 2 $<TARGET_FILE:proj> ${CMAKE_BINARY_DIR}/tests/SedovExplosionHllcRef_medium SedovCartesian_hllc)
# add_test(NAME medium_compare COMMAND ${HDF5_DIFF_EXECUTABLE} -v -d 1.5e-11 --exclude-path /Data/version ${CMAKE_BINARY_DIR}/tests/SedovExplosionHllcRef_medium/SedovCartesian_hllc_float/SedovCartesian_hllc_flt_step29.h5 ${CMAKE_SOURCE_DIR}/tests/SedovExplosionHllcRef_medium/SedovCartesian_hllc_float/SedovCartesian_hllc_flt_step29.h5)
# add_test(NAME medium_compare COMMAND ${HDF5_DIFF_EXECUTABLE} -v -d 1.5e-11 --exclude-path /Data/version ${CMAKE_BINARY_DIR}/tests/SedovExplosionHllcRef_medium/SedovCartesian_hllc_double/SedovCartesian_hllc_dbl_step1.h5 ${CMAKE_SOURCE_DIR}/tests/SedovExplosionHllcRef_medium/SedovCartesian_hllc_double/SedovCartesian_hllc_dbl_step1.h5)

# set_tests_properties(tiny_run PROPERTIES DEPENDS proj)
# set_tests_properties(tiny_run PROPERTIES FIXTURES_SETUP EXEC_TINY)
# set_tests_properties(tiny_compare PROPERTIES FIXTURES_REQUIRED EXEC_TINY)
set_tests_properties(small1_run PROPERTIES DEPENDS proj)
set_tests_properties(small1_run PROPERTIES FIXTURES_SETUP EXEC1_SMALL)
set_tests_properties(small1_compare PROPERTIES FIXTURES_REQUIRED EXEC1_SMALL)

set_tests_properties(small2_run PROPERTIES DEPENDS proj)
set_tests_properties(small2_run PROPERTIES FIXTURES_SETUP EXEC2_SMALL)
set_tests_properties(small2_compare PROPERTIES FIXTURES_REQUIRED EXEC2_SMALL)

set_tests_properties(small4_run PROPERTIES DEPENDS proj)
set_tests_properties(small4_run PROPERTIES FIXTURES_SETUP EXEC4_SMALL)
set_tests_properties(small4_compare PROPERTIES FIXTURES_REQUIRED EXEC4_SMALL)

set_tests_properties(small8_run PROPERTIES DEPENDS proj)
set_tests_properties(small8_run PROPERTIES FIXTURES_SETUP EXEC8_SMALL)
set_tests_properties(small8_compare PROPERTIES FIXTURES_REQUIRED EXEC8_SMALL)
# set_tests_properties(medium_run PROPERTIES DEPENDS proj)
# set_tests_properties(medium_run PROPERTIES FIXTURES_SETUP EXEC_MEDIUM)
# set_tests_properties(medium_compare PROPERTIES FIXTURES_REQUIRED EXEC_MEDIUM)