###############################################################################
# Find CronosNumLib
#
# This sets the following variables:
# CRONOSNUMLIB_FOUND - True if CronosNumLib was found.
# CRONOSNUMLIB_INCLUDE_DIR - Directories containing the CronosNumLib include files.
# CRONOSNUMLIB_MATRIX - Directory containing the matrix library for CronosNumLib.
# CRONOSNUMLIB_UTIL - Directory containing the util library for CronosNumLib.

cmake_minimum_required (VERSION 3.16.3)
include(FindPackageHandleStandardArgs)

SET(CMAKE_FIND_LIBRARY_PREFIXES "lib")
SET(CMAKE_FIND_LIBRARY_SUFFIXES ".so" ".a")

macro(CRONOSNUMLIB_REPORT_NOT_FOUND REASON_MSG)
    unset(CRONOSNUMLIB_FOUND)
    unset(CRONOSNUMLIB_INCLUDE_DIR)
    unset(CRONOSNUMLIB_MATRIX)
    unset(CRONOSNUMLIB_UTIL)

    if(CRONOSNUMLIB_FIND_QUIETLY)
        message(STATUS "Failed to find CronosNumLib - " ${REASON_MSG} ${ARGN})
    elseif(CRONOSNUMLIB_FIND_REQUIRED)
        message(FATAL_ERROR "Failed to find CronosNumLib - " ${REASON_MSG} ${ARGN})
    else()
        message("-- Failed to find CronosNumLib - " ${REASON_MSG} ${ARGN})
    endif()
endmacro(CRONOSNUMLIB_REPORT_NOT_FOUND)

list(APPEND CRONOSNUMLIB_CHECK_INCLUDE_DIRS
  /usr/local/include
  /usr/local/homebrew/include
  /opt/local/var/macports/software
  /opt/local/include
  /usr/include
  ${CMAKE_CURRENT_SOURCE_DIR}/../external/CronosNumLib/include
  ${CMAKE_CURRENT_SOURCE_DIR}/dependencies/include
  )

list(APPEND CRONOSNUMLIB_CHECK_LIBRARY_DIRS
  /usr/local/lib
  /usr/local/lib/CronosNumLib/Linux-amd64
  /usr/local/homebrew/lib
  /usr/local/homebrew/lib/CronosNumLib/Linux-amd64
  /opt/local/lib
  /opt/local/lib/CronosNumLib/Linux-amd64
  /usr/lib
  /usr/lib/CronosNumLib/Linux-amd64
  /usr/lib/x86_64-linux-gnu
  /usr/lib/x86_64-linux-gnu/CronosNumLib/Linux-amd64
  ${CMAKE_CURRENT_SOURCE_DIR}/../external/CronosNumLib/lib
  ${CMAKE_SOURCE_DIR}/dependencies/lib
  ${CMAKE_CURRENT_SOURCE_DIR}/dependencies/lib/CronosNumLib/Linux-amd64
  )

find_path(CRONOSNUMLIB_INCLUDE_DIR
  "CronosNumLib"
  PATHS ${CRONOSNUMLIB_CHECK_INCLUDE_DIRS}
  NO_DEFAULT_PATH
)

find_library(CRONOSNUMLIB_MATRIX
  matrix_mt
  HINTS ${CRONOSNUMLIB_CHECK_LIBRARY_DIRS}
  NO_DEFAULT_PATH
)

find_library(CRONOSNUMLIB_UTIL
  util_mt
  HINTS ${CRONOSNUMLIB_CHECK_LIBRARY_DIRS}
  NO_DEFAULT_PATH
)

set(CRONOSNUMLIB_FOUND TRUE)

if(NOT CRONOSNUMLIB_INCLUDE_DIR OR NOT EXISTS ${CRONOSNUMLIB_INCLUDE_DIR})
  CRONOSNUMLIB_REPORT_NOT_FOUND("Could not find CronosNumLib include directory.")
endif()
if(NOT CRONOSNUMLIB_MATRIX OR NOT EXISTS ${CRONOSNUMLIB_MATRIX})
  CRONOSNUMLIB_REPORT_NOT_FOUND("Could not find CronosNumLib Matrix library.")
endif()
if(NOT CRONOSNUMLIB_UTIL OR NOT EXISTS ${CRONOSNUMLIB_UTIL})
  CRONOSNUMLIB_REPORT_NOT_FOUND("Could not find CronosNumLib Util library.")
endif()

Find_package_handle_standard_args( CronosNumLib DEFAULT_MSG
  CRONOSNUMLIB_INCLUDE_DIR CRONOSNUMLIB_MATRIX CRONOSNUMLIB_UTIL
)

if(CRONOSNUMLIB_FOUND)
  message("-- Found CronosNumLib: " ${CRONOSNUMLIB_INCLUDE_DIR} "/CronosNumLib")
  mark_as_advanced(FORCE CRONOSNUMLIB_INCLUDE_DIR CRONOSNUMLIB_MATRIX CRONOSNUMLIB_UTIL)
endif()