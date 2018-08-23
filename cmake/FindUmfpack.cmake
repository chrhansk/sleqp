# Umfpack lib usually requires linking to a BLAS library.
# It is up to the user of this module to find a BLAS and link to it.
#
# Once done, this will define
#
#  UMFPACK_INCLUDE_DIRS   - where to find umfpack.h, etc.
#  UMFPACK_LIBRARIES      - List of libraries when using Umfpack.
#  UMFPACK_FOUND          - True if Umfpack found.

function(extract_define file name result)
  file(STRINGS "${file}"
    file_result
    REGEX "^#define ${name} +[0-9]+")
  string(REGEX REPLACE "^#define ${name} +([0-9]+).*" "\\1" replace_result ${file_result})
  set(${result} ${replace_result} PARENT_SCOPE)
endfunction()

find_path(UMFPACK_INCLUDE_DIRS
  NAMES umfpack.h
  PATHS $ENV{UMFPACKDIR} ${INCLUDE_INSTALL_DIR}
  PATH_SUFFIXES suitesparse ufsparse)

if(UMFPACK_INCLUDE_DIRS)
  extract_define("${UMFPACK_INCLUDE_DIRS}/umfpack.h"
    "UMFPACK_MAIN_VERSION"
    UMFPACK_MAIN_VERSION)
  extract_define("${UMFPACK_INCLUDE_DIRS}/umfpack.h"
    "UMFPACK_SUB_VERSION"
    UMFPACK_SUB_VERSION)
  extract_define("${UMFPACK_INCLUDE_DIRS}/umfpack.h"
    "UMFPACK_SUBSUB_VERSION"
    UMFPACK_SUBSUB_VERSION)

  set(UMFPACK_VERSION "${UMFPACK_MAIN_VERSION}.${UMFPACK_SUB_VERSION}.${UMFPACK_SUBSUB_VERSION}")
endif()

find_library(UMFPACK_LIBRARIES umfpack PATHS $ENV{UMFPACKDIR} ${LIB_INSTALL_DIR})

if(UMFPACK_LIBRARIES)

  if(NOT UMFPACK_LIBDIR)
    get_filename_component(UMFPACK_LIBDIR ${UMFPACK_LIBRARIES} PATH)
  endif(NOT UMFPACK_LIBDIR)

  find_library(COLAMD_LIBRARY colamd PATHS ${UMFPACK_LIBDIR} $ENV{UMFPACKDIR} ${LIB_INSTALL_DIR})
  if(COLAMD_LIBRARY)
    set(UMFPACK_LIBRARIES ${UMFPACK_LIBRARIES} ${COLAMD_LIBRARY})
  endif()

  find_library(AMD_LIBRARY amd PATHS ${UMFPACK_LIBDIR} $ENV{UMFPACKDIR} ${LIB_INSTALL_DIR})
  if(AMD_LIBRARY)
    set(UMFPACK_LIBRARIES ${UMFPACK_LIBRARIES} ${AMD_LIBRARY})
  endif()

  find_library(SUITESPARSE_LIBRARY SuiteSparse PATHS ${UMFPACK_LIBDIR} $ENV{UMFPACKDIR} ${LIB_INSTALL_DIR})
  if(SUITESPARSE_LIBRARY)
    set(UMFPACK_LIBRARIES ${UMFPACK_LIBRARIES} ${SUITESPARSE_LIBRARY})
  endif()

  find_library(CHOLMOD_LIBRARY cholmod PATHS $ENV{UMFPACK_LIBDIR} $ENV{UMFPACKDIR} ${LIB_INSTALL_DIR})
  if(CHOLMOD_LIBRARY)
    set(UMFPACK_LIBRARIES ${UMFPACK_LIBRARIES} ${CHOLMOD_LIBRARY})
  endif()

endif(UMFPACK_LIBRARIES)

include(FindPackageHandleStandardArgs)
find_package_handle_standard_args(Umfpack
  REQUIRED_VARS UMFPACK_INCLUDE_DIRS UMFPACK_LIBRARIES
  VERSION_VAR UMFPACK_VERSION)

mark_as_advanced(UMFPACK_INCLUDE_DIRS
  UMFPACK_LIBRARIES
  UMFPACK_VERSION
  AMD_LIBRARY
  COLAMD_LIBRARY
  CHOLMOD_LIBRARY
  SUITESPARSE_LIBRARY)
