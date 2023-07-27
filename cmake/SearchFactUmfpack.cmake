# Umfpack lib usually requires linking to a BLAS library.
# It is up to the user of this module to find a BLAS and link to it.
#
# Once done, this will define
#
#  UMFPACK_INCLUDE_DIRS   - where to find umfpack.h, etc.
#  UMFPACK_LIBRARIES      - List of libraries when using Umfpack.
#  UMFPACK_FOUND          - True if Umfpack found.

find_path(UMFPACK_INCLUDE_DIRS
  NAMES umfpack.h
  PATHS $ENV{UMFPACKDIR} ${INCLUDE_INSTALL_DIR}
  PATH_SUFFIXES suitesparse ufsparse)

if(UMFPACK_INCLUDE_DIRS)
  define_extract_int("${UMFPACK_INCLUDE_DIRS}/umfpack.h"
    UMFPACK_MAIN_VERSION)
  define_extract_int("${UMFPACK_INCLUDE_DIRS}/umfpack.h"
    UMFPACK_SUB_VERSION)
  define_extract_int("${UMFPACK_INCLUDE_DIRS}/umfpack.h"
    UMFPACK_SUBSUB_VERSION)

  set(UMFPACK_VERSION "${UMFPACK_MAIN_VERSION}.${UMFPACK_SUB_VERSION}.${UMFPACK_SUBSUB_VERSION}")
endif()

find_library(UMFPACK_LIBRARY umfpack PATHS $ENV{UMFPACKDIR} ${LIB_INSTALL_DIR})

set(UMFPACK_LIBRARIES "${UMFPACK_LIBRARY}")

if(UMFPACK_LIBRARY)

  if(NOT UMFPACK_LIBDIR)
    get_filename_component(UMFPACK_LIBDIR ${UMFPACK_LIBRARY} PATH)
  endif(NOT UMFPACK_LIBDIR)

  if(NOT(OpenMP_FOUND))
    find_package(OpenMP)
  endif()

  if(OpenMP_FOUND)
    set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${OpenMP_C_FLAGS}")
  endif()

endif(UMFPACK_LIBRARY)

include(FindPackageHandleStandardArgs)
find_package_handle_standard_args(Umfpack
  REQUIRED_VARS UMFPACK_INCLUDE_DIRS UMFPACK_LIBRARIES
  VERSION_VAR UMFPACK_VERSION)

mark_as_advanced(UMFPACK_INCLUDE_DIRS
  UMFPACK_LIBRARIES
  UMFPACK_VERSION)
