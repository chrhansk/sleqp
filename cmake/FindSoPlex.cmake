# - Try to find SOPLEX
# See http://soplex.zib.de/ for more information on SoPlex
#
# Once done, this will define
#
#  SOPLEX_INCLUDE_DIRS   - where to find soplex.h, etc.
#  SOPLEX_LIBRARIES      - List of libraries when using soplex.
#  SOPLEX_FOUND          - True if soplex found.
#
# A maintainer may set SOPLEX_ROOT to a soplex installation root to tell
# this module where to look.

find_path(SOPLEX_INCLUDE_DIR
  NAMES soplex.h
  PATHS ${SOPLEX_ROOT} ${SCIP_ROOT} "${SCIP_ROOT}/soplex"
  PATH_SUFFIXES src include soplex)

# Try to find a PIC-version first
find_library(SOPLEX_LIBRARY
  NAMES soplex-pic soplex
  PATHS ${SOPLEX_ROOT} ${SCIP_ROOT}
  PATH_SUFFIXES lib)

if(SOPLEX_INCLUDE_DIR)
  find_file(SPX_HEADER
    NAMES "spxdefines.h"
    PATHS ${SOPLEX_INCLUDE_DIR}
    PATH_SUFFIXES "soplex")

  if(SPX_HEADER)
    file(STRINGS "${SPX_HEADER}" SOPLEX_DEF_H REGEX "^#define SOPLEX_VERSION +[0-9]+")
    string(REGEX REPLACE "^#define SOPLEX_VERSION +([0-9]+).*" "\\1" SVER ${SOPLEX_DEF_H})

    string(REGEX REPLACE "([0-9]).*" "\\1" SOPLEX_VERSION_MAJOR ${SVER})
    string(REGEX REPLACE "[0-9]([0-9]).*" "\\1" SOPLEX_VERSION_MINOR ${SVER})
    string(REGEX REPLACE "[0-9][0-9]([0-9]).*" "\\1" SOPLEX_VERSION_PATCH ${SVER})
    set(SOPLEX_VERSION "${SOPLEX_VERSION_MAJOR}.${SOPLEX_VERSION_MINOR}.${SOPLEX_VERSION_PATCH}")
  endif()
endif()

find_package(GMP QUIET)
find_package(ZLIB QUIET)

find_package_handle_standard_args(SOPLEX
  FOUND_VAR SOPLEX_FOUND
  REQUIRED_VARS SOPLEX_INCLUDE_DIR SOPLEX_LIBRARY GMP_LIBRARIES ZLIB_LIBRARIES
  VERSION_VAR SOPLEX_VERSION)

set(SOPLEX_INCLUDE_DIRS ${SOPLEX_INCLUDE_DIR} ${GMP_INCLUDE_DIR} ${ZLIB_INCLUDE_DIRS})
set(SOPLEX_LIBRARIES ${SOPLEX_LIBRARY} ${GMP_LIBRARIES} ${ZLIB_LIBRARIES})

mark_as_advanced(SOPLEX_INCLUDE_DIR SOPLEX_INCLUDE_DIRS SOPLEX_LIBRARIES SOPLEX_LIBRARY GMP_LIBRARIES ZLIB_LIBRARIES)
