cmake_minimum_required (VERSION 3.8)
project (CronosNumLib CXX)

set(CMAKE_CXX_STANDARD 14)

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
	# add_compile_options("-DH5_BUILT_AS_DYNAMIC_LIB")
	# enable debug information in object files
	add_compile_options("$<$<CONFIG:DEBUG>:/Z7;/DEBUG>")
	# set iterator checking
	add_compile_options("$<$<CONFIG:DEBUG>:-D_ITERATOR_DEBUG_LEVEL=2>")
	add_compile_options("$<$<CONFIG:RELEASE>:-D_ITERATOR_DEBUG_LEVEL=0>")
	# enable debug information in object files
	#set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} /MD")

endif()


find_package(MPI REQUIRED)

include_directories(Matrix)
include_directories(util)
include_directories(SYSTEM ${MPI_INCLUDE_PATH})

#Matrix/matrix_1d_inc.C 
#Matrix/matrix_2d_inc.C 
#Matrix/matrix_3d_inc.C 
#Matrix/num_boundary_1d_inc.C 
#Matrix/num_boundary_2d_inc.C 
#Matrix/num_boundary_3d_inc.C 
#Matrix/num_matrix_1d_inc.C 
#Matrix/num_matrix_2d_inc.C 
#Matrix/num_matrix_3d_inc.C

set(MATRIX_SOURCES
	Matrix/array.C 
	Matrix/boundary.C 
	Matrix/fftmatrix.C 
	Matrix/matrix.C 
)

FILE(GLOB_RECURSE MATRIX_HEADERS "Matrix/*.h")

add_library(matrix ${MATRIX_SOURCES} ${MATRIX_HEADERS})
target_link_libraries(matrix PRIVATE ${MPI_C_LIBRARIES})

set_property(TARGET matrix PROPERTY CXX_STANDARD 14)

set(UTIL_SOURCES
	util/CronosOstream.C
	util/ParameterFileReader.C
	util/util.C
)

FILE(GLOB_RECURSE UTIL_HEADERS "util/*.h")

add_library(util ${UTIL_SOURCES} ${UTIL_HEADERS})
target_link_libraries(util PRIVATE ${MPI_C_LIBRARIES})

set_property(TARGET util PROPERTY CXX_STANDARD 14)