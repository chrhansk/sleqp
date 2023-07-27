# Once done, this will define
#
#  SPQR_INCLUDE_DIRS   - where to find SuiteSparseQR_C.h, etc.
#  SPQR_LIBRARIES      - List of libraries when using SPQR.
#  SPQR_FOUND          - True if SPQR found.

find_path(SPQR_INCLUDE_DIRS
  NAMES SuiteSparseQR_C.h
  PATHS $ENV{SPQR_DIR} ${INCLUDE_INSTALL_DIR}
  PATH_SUFFIXES suitesparse ufsparse)

if(SPQR_INCLUDE_DIRS)
  define_extract_int("${SPQR_INCLUDE_DIRS}/SuiteSparseQR_definitions.h"
    SPQR_MAIN_VERSION)
  define_extract_int("${SPQR_INCLUDE_DIRS}/SuiteSparseQR_definitions.h"    "SPQR_SUB_VERSION"
    SPQR_SUB_VERSION)
  define_extract_int("${SPQR_INCLUDE_DIRS}/SuiteSparseQR_definitions.h"
    SPQR_SUBSUB_VERSION)

  set(SPQR_VERSION "${SPQR_MAIN_VERSION}.${SPQR_SUB_VERSION}.${SPQR_SUBSUB_VERSION}")
endif()

find_library(SPQR_LIBRARY spqr PATHS $ENV{SPQR_DIR} ${LIB_INSTALL_DIR})
find_library(CHOLMOD_LIBRARY cholmod PATHS $ENV{SPQR_DIR} ${LIB_INSTALL_DIR})

set(SPQR_LIBRARIES "${SPQR_LIBRARY};${CHOLMOD_LIBRARY}")

if(SPQR_LIBRARY)
  if(NOT(OpenMP_FOUND))
    find_package(OpenMP)
  endif()

  if(OpenMP_FOUND)
    set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${OpenMP_C_FLAGS}")
  endif()

endif(SPQR_LIBRARY)

include(FindPackageHandleStandardArgs)
find_package_handle_standard_args(SPQR
  REQUIRED_VARS SPQR_INCLUDE_DIRS SPQR_LIBRARIES
  VERSION_VAR SPQR_VERSION)

mark_as_advanced(SPQR_INCLUDE_DIRS
  SPQR_LIBRARIES
  SPQR_VERSION)
