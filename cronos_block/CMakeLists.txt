cmake_minimum_required (VERSION 3.16.3)
project (CronosCodeBlock C CXX)
enable_language(Fortran)

set(CMAKE_CXX_STANDARD 17)

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


endif()

############################
# Third-party Dependencies #
############################

list(APPEND CMAKE_MODULE_PATH ${CMAKE_SOURCE_DIR}/cmake)
find_package(GSL REQUIRED)
find_package(MPI REQUIRED)
find_package(HDF5 REQUIRED COMPONENTS C) # CXX HL
find_package(CronosNumLib REQUIRED)

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
include_directories(${CRONOSNUMLIB_INCLUDE_DIR}/CronosNumLib)

include_directories(test)
include_directories(general)
include_directories(interface_functions)
include_directories(kernel)
include_directories(riemann_solver)
include_directories(reconstruction)
include_directories(utils)
include_directories(configuration/ShockTubeSod)
include_directories(configuration/standalone/cronos_block)

set(SOURCES
    general/buildinfo.C
    general/CException.C
    general/cfl.C
    general/data.C
    general/eos.C
    general/fluid.C
    general/grid.C
    general/gridfunc.C
    general/gridgen.C
    general/gridParam.C
    general/Hdf5File_cbase.C
    general/movie.C
    general/normalisation.C
    general/proj.C
    general/sources.C
    general/specific.C
    general/timestepping.C
    general/units.C
    interface_functions/problem.C
    interface_functions/solver.C
    kernel/fields_local.C
    kernel/singlestep_block.C
    kernel/transformations.C
    riemann_solver/changes.C
    riemann_solver/hllc_hd.C
    riemann_solver/PhysFluxesHD.C
    riemann_solver/RiemannSolver.C
    riemann_solver/vChar_hd.C
    reconstruction/reconst.C
    #test/proj.C
    #test/test_data.C
    #test/singlestep_stub.C
    configuration/ShockTubeSod/modules.C
)

FILE(GLOB_RECURSE HEADERS "*.h")

add_executable(proj ${SOURCES} general/proj.C ${HEADERS})
#add_executable(proj ${SOURCES} test/proj.C ${HEADERS})

set_property(TARGET proj PROPERTY CXX_STANDARD 17)

# properly set up SYCL compiler in case of hipSYCL

set(DEVICE_SOURCES "")
if(CELERITY_SYCL_IMPL STREQUAL "hipSYCL")
  # For hipSYCL we have to pass all source files into add_sycl_to_target
  # Don't just do it in general (yields better ComputeCpp build performance)
  set(DEVICE_SOURCES ${SOURCES} cronos/generic/proj.C cronos/generic/playground.C)
endif()
add_sycl_to_target(TARGET proj SOURCES ${DEVICE_SOURCES})

#add_celerity_to_target(TARGET proj SOURCES cronos/generic/proj.C ${HEADERS})

if(NOT MSVC)
	target_link_libraries(proj PRIVATE stdc++fs)
endif()

target_link_libraries(proj PUBLIC ${CRONOSNUMLIB_MATRIX} ${CRONOSNUMLIB_UTIL})
target_link_libraries(proj PRIVATE ${MPI_C_LIBRARIES})
target_link_libraries(proj PRIVATE GSL::gsl GSL::gslcblas)
target_link_libraries(proj PRIVATE ${HDF5_LIBRARIES})

#########
# TESTS #
#########

enable_testing()

configure_file(${CMAKE_SOURCE_DIR}/tests/SodToro5HllcRef_tiny/SodToro5HllcRef.cat ${CMAKE_BINARY_DIR}/tests/SodToro5HllcRef_tiny/SodToro5HllcRef.cat COPYONLY)
configure_file(${CMAKE_SOURCE_DIR}/tests/SodToro5HllcRef_small/SodToro5HllcRef.cat ${CMAKE_BINARY_DIR}/tests/SodToro5HllcRef_small/SodToro5HllcRef.cat COPYONLY)
configure_file(${CMAKE_SOURCE_DIR}/tests/SodToro5HllcRef_medium/SodToro5HllcRef.cat ${CMAKE_BINARY_DIR}/tests/SodToro5HllcRef_medium/SodToro5HllcRef.cat COPYONLY)

add_test(NAME tiny_run COMMAND $<TARGET_FILE:proj> ${CMAKE_BINARY_DIR}/tests/SodToro5HllcRef_tiny SodToro5HllcRef)
add_test(NAME tiny_compare COMMAND ${HDF5_DIFF_EXECUTABLE} --exclude-path /Data/version ${CMAKE_BINARY_DIR}/tests/SodToro5HllcRef_tiny/SodToro5HllcRef_float/SodToro5HllcRef_flt_step3.h5 ${CMAKE_SOURCE_DIR}/tests/SodToro5HllcRef_tiny/SodToro5HllcRef_float/SodToro5HllcRef_flt_step3.h5)
	
add_test(NAME small_run COMMAND $<TARGET_FILE:proj> ${CMAKE_BINARY_DIR}/tests/SodToro5HllcRef_small SodToro5HllcRef)
add_test(NAME small_compare COMMAND ${HDF5_DIFF_EXECUTABLE} --exclude-path /Data/version ${CMAKE_BINARY_DIR}/tests/SodToro5HllcRef_small/SodToro5HllcRef_float/SodToro5HllcRef_flt_step16.h5 ${CMAKE_SOURCE_DIR}/tests/SodToro5HllcRef_small/SodToro5HllcRef_float/SodToro5HllcRef_flt_step16.h5)

add_test(NAME medium_run COMMAND $<TARGET_FILE:proj> ${CMAKE_BINARY_DIR}/tests/SodToro5HllcRef_medium SodToro5HllcRef)
add_test(NAME medium_compare COMMAND ${HDF5_DIFF_EXECUTABLE} --exclude-path /Data/version ${CMAKE_BINARY_DIR}/tests/SodToro5HllcRef_medium/SodToro5HllcRef_float/SodToro5HllcRef_flt_step990.h5 ${CMAKE_SOURCE_DIR}/tests/SodToro5HllcRef_medium/SodToro5HllcRef_float/SodToro5HllcRef_flt_step990.h5)

set_tests_properties(tiny_run PROPERTIES DEPENDS proj)
set_tests_properties(tiny_run PROPERTIES FIXTURES_SETUP EXEC_TINY)
set_tests_properties(tiny_compare PROPERTIES FIXTURES_REQUIRED EXEC_TINY)
	
set_tests_properties(small_run PROPERTIES DEPENDS proj)
set_tests_properties(small_run PROPERTIES FIXTURES_SETUP EXEC_SMALL)
set_tests_properties(small_compare PROPERTIES FIXTURES_REQUIRED EXEC_SMALL)
	
set_tests_properties(medium_run PROPERTIES DEPENDS proj)
set_tests_properties(medium_run PROPERTIES FIXTURES_SETUP EXEC_MEDIUM)
set_tests_properties(medium_compare PROPERTIES FIXTURES_REQUIRED EXEC_MEDIUM)