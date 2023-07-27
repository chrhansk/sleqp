# Once done, this will define
#
#  CHOLMOD_INCLUDE_DIRS   - where to find cholmod.h, etc.
#  CHOLMOD_LIBRARIES      - List of libraries when using CHOLMOD.
#  CHOLMOD_FOUND          - True if CHOLMOD found.

find_path(CHOLMOD_INCLUDE_DIRS
  NAMES cholmod.h
  PATHS $ENV{CHOLMOD_DIR} ${INCLUDE_INSTALL_DIR}
  PATH_SUFFIXES suitesparse)

find_library(CHOLMOD_LIBRARY cholmod PATHS $ENV{CHOLMOD_DIR} ${LIB_INSTALL_DIR})

set(CHOLMOD_LIBRARIES "${CHOLMOD_LIBRARY}")

include(FindPackageHandleStandardArgs)

find_package_handle_standard_args(CHOLMOD
  REQUIRED_VARS CHOLMOD_INCLUDE_DIRS CHOLMOD_LIBRARIES)

mark_as_advanced(CHOLMOD_INCLUDE_DIRS
  CHOLMOD_LIBRARIES)
