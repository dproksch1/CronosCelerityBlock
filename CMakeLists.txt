cmake_minimum_required (VERSION 3.16.3)
project (CronosCode C CXX)
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


#message(${CMAKE_SOURCE_DIR}/cmake)
list(APPEND CMAKE_MODULE_PATH ${CMAKE_SOURCE_DIR}/cmake)
include(ExternalProject)
find_package(GSL REQUIRED)
find_package(MPI REQUIRED)
find_package(HDF5 REQUIRED COMPONENTS C) # CXX HL

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

execute_process(COMMAND "git" "--no-pager" "describe" "--dirty" "--always" "--tags" WORKING_DIRECTORY ${PROJECT_SOURCE_DIR} OUTPUT_VARIABLE CRONOS_GIT_VERSION OUTPUT_STRIP_TRAILING_WHITESPACE)
execute_process(COMMAND "git" "rev-parse" "--verify" "HEAD" WORKING_DIRECTORY ${PROJECT_SOURCE_DIR} OUTPUT_VARIABLE CRONOS_GIT_COMMIT OUTPUT_STRIP_TRAILING_WHITESPACE)
execute_process(COMMAND "git" "remote" "get-url" "--all" "origin" WORKING_DIRECTORY ${PROJECT_SOURCE_DIR} OUTPUT_VARIABLE APP_GIT_REPO OUTPUT_STRIP_TRAILING_WHITESPACE)
execute_process(COMMAND "git" "--no-pager" "describe" "--dirty" "--always" "--tags" WORKING_DIRECTORY ${PROJECT_SOURCE_DIR} OUTPUT_VARIABLE APP_GIT_VERSION OUTPUT_STRIP_TRAILING_WHITESPACE)
execute_process(COMMAND "git" "rev-parse" "--verify" "HEAD" WORKING_DIRECTORY ${PROJECT_SOURCE_DIR} OUTPUT_VARIABLE APP_GIT_COMMIT OUTPUT_STRIP_TRAILING_WHITESPACE)
string(TIMESTAMP BUILD_DATE)

add_compile_definitions(CRONOS_ROOT_DIR=\"${PROJECT_SOURCE_DIR}\")
add_compile_definitions(CRONOS_GIT_VERSION=\"${CRONOS_GIT_VERSION}\")
add_compile_definitions(CRONOS_GIT_COMMIT=\"${CRONOS_GIT_COMMIT}\")
add_compile_definitions(APP_GIT_REPO=\"${APP_GIT_REPO}\")
add_compile_definitions(APP_GIT_VERSION=\"${APP_GIT_VERSION}\")
add_compile_definitions(APP_GIT_COMMIT=\"${APP_GIT_COMMIT}\")
add_compile_definitions(BUILD_DATE=\"${BUILD_DATE}\")

# add_compile_definitions(parallel)

add_compile_definitions(CL_TARGET_OPENCL_VERSION=220)

# import numeric library

add_subdirectory(external/CronosNumLib)

# add include directories

include_directories(${HDF5_INCLUDE_DIRS})
include_directories(${GSL_INCLUDE_DIRS})

# TODO: properly lookup package via find cmake or at least a root_path variable
#include_directories(../CronosNumLib/Matrix)
#include_directories(../CronosNumLib/util)
#include_directories(external/CronosNumLib/Matrix)
#include_directories(external/CronosNumLib/util)
find_package(CronosNumLib REQUIRED)
include_directories(${CRONOSNUMLIB_INCLUDE_DIR}/CronosNumLib)

include_directories(cronos/DataReader        )
include_directories(cronos/Geom              )
include_directories(cronos/InterfaceFunctions)
include_directories(cronos/Limiter           )
include_directories(cronos/RiemannSolvers    )
include_directories(cronos/RiemannSolvers/HD    )
include_directories(cronos/RiemannSolvers/MHD    )
include_directories(cronos/StandardTests     )
include_directories(cronos/Testing           )
# You can always only include 1 -> find solutions for that
include_directories(cronos/Testing/ShockTubeSod)
#include_directories(cronos/Testing/SedovExplosion)
#include_directories(cronos/Testing/BrioWu)
include_directories(cronos/UnitTests         )
include_directories(cronos/User              )
include_directories(cronos/Util              )
include_directories(cronos/generic           )
include_directories(cronos/release_docs      )

include_directories(SYSTEM ${MPI_INCLUDE_PATH})

# add source directories
# potentially look to remove some unused sources
set(SOURCES
	cronos/Limiter/limiter.C
	cronos/generic/Hdf5File_cbase.C
	cronos/generic/buildinfo.C
	cronos/generic/movie.C
	cronos/generic/CException.C
	cronos/generic/gridParam.C
	cronos/generic/eos.C
	cronos/generic/normalisation.C
	cronos/generic/units.C
	cronos/generic/gridgen.C
	cronos/generic/data.C
	cronos/generic/multifluid.C
	cronos/generic/specific.C
	cronos/generic/grid.C
	cronos/generic/gridfunc.C
	cronos/generic/timestepping.C
	cronos/generic/randgen.C
	cronos/generic/user_equations.C
	cronos/generic/fluid.C
	cronos/generic/queue.C
	cronos/InterfaceFunctions/problem.C
	#cronos/InterfaceFunctions/singlestep.C
	cronos/Geom/sources.C
	#cronos/RiemannSolvers/reconst_const.C
	cronos/RiemannSolvers/transformations.C
	cronos/RiemannSolvers/save_data.C
	cronos/RiemannSolvers/reconst_2nd.C
	cronos/RiemannSolvers/fieldLists.C
	cronos/RiemannSolvers/reconst.C
	cronos/RiemannSolvers/fields_local.C
	cronos/RiemannSolvers/HD/hllc_hd.C
	cronos/RiemannSolvers/HD/PhysFluxesHD.C
	cronos/RiemannSolvers/HD/hll_hd.C
	cronos/RiemannSolvers/HD/vChar_hd.C
	cronos/RiemannSolvers/RiemannSolver.C
	cronos/RiemannSolvers/MHD/DissipationMHD.C
	cronos/RiemannSolvers/MHD/ConstrainedTransport.C
	cronos/RiemannSolvers/MHD/CTStone.C
	cronos/RiemannSolvers/MHD/CTLondrilloDelZanna.C
	cronos/RiemannSolvers/MHD/vChar.C
	cronos/RiemannSolvers/MHD/hll.C
	cronos/RiemannSolvers/MHD/hlld.C
	cronos/RiemannSolvers/MHD/PhysFluxesMHD.C #commented out
	cronos/RiemannSolvers/MHD/bc_emf.C
	cronos/RiemannSolvers/cfl.C
	cronos/RiemannSolvers/changes.C
	#cronos/RiemannSolvers/reconst_WENO.C
	#cronos/RiemannSolvers/singlestep_general.C
	cronos/RiemannSolvers/singlestep_block.C
	# application-specific source files
	cronos/Testing/ShockTubeSod/modules.C
	#cronos/Testing/SedovExplosion/modules.C
	#cronos/Testing/BrioWu/modules.C
	# fall-back files if no application-specific versions found
	cronos/InterfaceFunctions/solver.C
	cronos/InterfaceFunctions/user_pde.C
	#cronos/generic/Hdf5File.C # deprecated, see comment at the beginning of the source file
	#cronos/RiemannSolvers/singlestep_CTSecond.C # produces compiler error due to non-existing SourceTerms constructor signature
	#cronos/RiemannSolvers/singlestep_block.C # produces compiler error due to non-existing SourceTerms constructor signature

)

# glob for all headers (used for better MSVC integration)
# TODO: only include required headers, not from all apps?
FILE(GLOB_RECURSE HEADERS "*.h")

add_executable(proj ${SOURCES} cronos/generic/proj.C ${HEADERS})
add_executable(playground ${SOURCES} cronos/generic/playground.C ${HEADERS})

# set C++ language standard

set_property(TARGET proj PROPERTY CXX_STANDARD 17)
set_property(TARGET playground PROPERTY CXX_STANDARD 17)

# properly set up SYCL compiler in case of hipSYCL

set(DEVICE_SOURCES "")
if(CELERITY_SYCL_IMPL STREQUAL "hipSYCL")
  # For hipSYCL we have to pass all source files into add_sycl_to_target
  # Don't just do it in general (yields better ComputeCpp build performance)
  set(DEVICE_SOURCES ${SOURCES} cronos/generic/proj.C cronos/generic/playground.C)
endif()
add_sycl_to_target(TARGET proj SOURCES ${DEVICE_SOURCES})
add_sycl_to_target(TARGET playground SOURCES ${DEVICE_SOURCES})

add_celerity_to_target(TARGET proj SOURCES cronos/generic/proj.C ${HEADERS})
add_celerity_to_target(TARGET playground SOURCES cronos/generic/playground.C ${HEADERS})

# add std::filesystem linking, required by some compilers

if(NOT MSVC)
	target_link_libraries(proj PRIVATE stdc++fs)
	target_link_libraries(playground PRIVATE stdc++fs)
endif()

# link to numeric library

#target_link_libraries(proj PRIVATE matrix util)
#target_link_libraries(playground PRIVATE matrix util)
target_link_libraries(proj PUBLIC ${CRONOSNUMLIB_MATRIX} ${CRONOSNUMLIB_UTIL})
target_link_libraries(playground PUBLIC ${CRONOSNUMLIB_MATRIX} ${CRONOSNUMLIB_UTIL})

# link to third-party libraries

target_link_libraries(proj PRIVATE ${MPI_C_LIBRARIES})
target_link_libraries(proj PRIVATE GSL::gsl GSL::gslcblas)
target_link_libraries(proj PRIVATE ${HDF5_LIBRARIES})
target_link_libraries(playground PRIVATE ${MPI_C_LIBRARIES})
target_link_libraries(playground PRIVATE GSL::gsl GSL::gslcblas)
target_link_libraries(playground PRIVATE ${HDF5_LIBRARIES})
#target_link_libraries(proj PRIVATE 
#C:/dev/HDF5-1.12.0-win64/lib/libhdf5_cpp.lib 
#C:/dev/HDF5-1.12.0-win64/lib/libhdf5_hl_cpp.lib 
#C:/dev/HDF5-1.12.0-win64/lib/libhdf5_tools.lib 
#C:/dev/HDF5-1.12.0-win64/lib/libhdf5_hl.lib 
#C:/dev/HDF5-1.12.0-win64/lib/libhdf5.lib 
#C:/dev/HDF5-1.12.0-win64/lib/hdf5_cpp.lib 
#C:/dev/HDF5-1.12.0-win64/lib/hdf5_hl_cpp.lib 
#C:/dev/HDF5-1.12.0-win64/lib/hdf5_tools.lib 
#C:/dev/HDF5-1.12.0-win64/lib/hdf5_hl.lib 
#C:/dev/HDF5-1.12.0-win64/lib/hdf5.lib 
#)
#target_link_libraries(proj PRIVATE ${LINK_LIBS})


#########
# TESTS #
#########

enable_testing()

configure_file(${CMAKE_SOURCE_DIR}/tests/SodToro5HllcRef_tiny/SodToro5HllcRef.cat ${CMAKE_BINARY_DIR}/tests/SodToro5HllcRef_tiny/SodToro5HllcRef.cat COPYONLY)
configure_file(${CMAKE_SOURCE_DIR}/tests/SodToro5HllcRef_small/SodToro5HllcRef.cat ${CMAKE_BINARY_DIR}/tests/SodToro5HllcRef_small/SodToro5HllcRef.cat COPYONLY)
configure_file(${CMAKE_SOURCE_DIR}/tests/SodToro5HllcRef_medium/SodToro5HllcRef.cat ${CMAKE_BINARY_DIR}/tests/SodToro5HllcRef_medium/SodToro5HllcRef.cat COPYONLY)
#configure_file(${CMAKE_SOURCE_DIR}/tests/SedovCartHllcRef_medium/SedovCartHllcRef.cat ${CMAKE_BINARY_DIR}/tests/SedovCartHllcRef_medium/SedovCartHllcRef.cat COPYONLY)
#configure_file(${CMAKE_SOURCE_DIR}/tests/BrioWuHlldRef_medium/BrioWuHlldRef.cat ${CMAKE_BINARY_DIR}/tests/BrioWuHlldRef_medium/BrioWuHlldRef.cat COPYONLY)

#if(MSVC)
#	add_test(small ${CMAKE_BINARY_DIR}/proj.exe ${CMAKE_BINARY_DIR}/tests/SodToro5HllcRef SodToro5HllcRef && h5diff --exclude-path /Data/version ${CMAKE_BINARY_DIR}/tests/SodToro5HllcRef_float/SodToro5HllcRef_flt_step990_coord0-0-0.h5 ${CMAKE_SOURCE__DIR}/tests/SodToro5HllcRef/SodToro5HllcRef_float/SodToro5HllcRef_flt_step990.h5)
#else()

	#add_test(NAME tiny_clean REMOVE ${CMAKE_BINARY_DIR}/tests/SodToro5HllcRef_tiny/SodToro5HllcRef_float/SodToro5HllcRef_flt_step3.h5)
	add_test(NAME tiny_run COMMAND $<TARGET_FILE:proj> ${CMAKE_BINARY_DIR}/tests/SodToro5HllcRef_tiny SodToro5HllcRef)
	add_test(NAME tiny_compare COMMAND ${HDF5_DIFF_EXECUTABLE} --exclude-path /Data/version ${CMAKE_BINARY_DIR}/tests/SodToro5HllcRef_tiny/SodToro5HllcRef_float/SodToro5HllcRef_flt_step3.h5 ${CMAKE_SOURCE_DIR}/tests/SodToro5HllcRef_tiny/SodToro5HllcRef_float/SodToro5HllcRef_flt_step3.h5)
	
	add_test(NAME small_run COMMAND $<TARGET_FILE:proj> ${CMAKE_BINARY_DIR}/tests/SodToro5HllcRef_small SodToro5HllcRef)
	add_test(NAME small_compare COMMAND ${HDF5_DIFF_EXECUTABLE} --exclude-path /Data/version ${CMAKE_BINARY_DIR}/tests/SodToro5HllcRef_small/SodToro5HllcRef_float/SodToro5HllcRef_flt_step16.h5 ${CMAKE_SOURCE_DIR}/tests/SodToro5HllcRef_small/SodToro5HllcRef_float/SodToro5HllcRef_flt_step16.h5)

	add_test(NAME medium_run COMMAND $<TARGET_FILE:proj> ${CMAKE_BINARY_DIR}/tests/SodToro5HllcRef_medium SodToro5HllcRef)
	add_test(NAME medium_compare COMMAND ${HDF5_DIFF_EXECUTABLE} --exclude-path /Data/version ${CMAKE_BINARY_DIR}/tests/SodToro5HllcRef_medium/SodToro5HllcRef_float/SodToro5HllcRef_flt_step990.h5 ${CMAKE_SOURCE_DIR}/tests/SodToro5HllcRef_medium/SodToro5HllcRef_float/SodToro5HllcRef_flt_step990.h5)
	
	#add_test(NAME sedov_cart_medium_run COMMAND $<TARGET_FILE:proj> ${CMAKE_BINARY_DIR}/tests/SedovCartHllcRef_medium SedovCartHllcRef)
	#add_test(NAME sedov_cart_medium_compare COMMAND ${HDF5_DIFF_EXECUTABLE} --exclude-path /Data/version ${CMAKE_BINARY_DIR}/tests/SedovCartHllcRef_medium/SedovCartHllcRef_float/SedovCartHllcRef_flt_step2.h5 ${CMAKE_SOURCE_DIR}/tests/SedovCartHllcRef_medium/SedovCartHllcRef_float/SedovCartHllcRef_flt_step2.h5)

	#add_test(NAME brio_wu_medium_run COMMAND $<TARGET_FILE:proj> ${CMAKE_BINARY_DIR}/tests/BrioWuHlldRef_medium BrioWuHlldRef)

	set_tests_properties(tiny_run PROPERTIES DEPENDS proj)
	set_tests_properties(tiny_run PROPERTIES FIXTURES_SETUP EXEC_TINY)
	set_tests_properties(tiny_compare PROPERTIES FIXTURES_REQUIRED EXEC_TINY)
	
	set_tests_properties(small_run PROPERTIES DEPENDS proj)
	set_tests_properties(small_run PROPERTIES FIXTURES_SETUP EXEC_SMALL)
	set_tests_properties(small_compare PROPERTIES FIXTURES_REQUIRED EXEC_SMALL)
	
	set_tests_properties(medium_run PROPERTIES DEPENDS proj)
	set_tests_properties(medium_run PROPERTIES FIXTURES_SETUP EXEC_MEDIUM)
	set_tests_properties(medium_compare PROPERTIES FIXTURES_REQUIRED EXEC_MEDIUM)

	#set_tests_properties(sedov_cart_medium_run PROPERTIES DEPENDS proj)
	#set_tests_properties(sedov_cart_medium_run PROPERTIES FIXTURES_SETUP EXEC_SEDOV_CART_MEDIUM)
	#set_tests_properties(sedov_cart_medium_compare PROPERTIES FIXTURES_REQUIRED EXEC_SEDOV_CART_MEDIUM)

	#set_tests_properties(brio_wu_medium_run PROPERTIES DEPENDS proj)
	
	#test currently failing
	#add_test(NAME small COMMAND ${CMAKE_COMMAND} -DCMD1=$<TARGET_FILE:proj> -P ${CMAKE_SOURCE_DIR}/tests/SodToro5HllcRef/run_test.cmake)
	#set_tests_properties(small PROPERTIES FIXTURES_REQUIRED proj)
	
#endif()

##set_tests_properties(small_run PROPERTIES FAIL_REGULAR_EXPRESSION "Error;error;ERROR")
# The following command required CMake 3.17 or higher
##list(APPEND CMAKE_CTEST_ARGUMENTS "--output-on-failure")

#set_target_properties(proj PROPERTIES COMPILE_DEFINITIONS BUILDER_STATIC_DEFINE)

# TODO:
# * fix GSL debug linking for MSVC
