# - Try to find Gurobi
# See http://soplex.zib.de/ for more information on SoPlex
#
# Once done, this will define
#
#  GUROBI_INCLUDE_DIRS   - where to find gurobi_c.h, etc.
#  GUROBI_LIBRARIES      - List of libraries when using gurobi.
#  GUROBI_FOUND          - True if Gurobi found.
#
# A maintainer may set GUROBI_ROOT to a soplex installation root to tell
# this module where to look.

find_path(GUROBI_INCLUDE_DIR
  NAMES gurobi_c.h
  PATHS ${GUROBI_ROOT}
  PATH_SUFFIXES include)

# TODO: FInd a better way to get the correct version
find_library(GUROBI_LIBRARY
  NAMES gurobi gurobi80
  PATHS ${GUROBI_ROOT}
  PATH_SUFFIXES lib)

if(GUROBI_INCLUDE_DIR)
  find_file(GUROBI_HEADER
    NAMES "gurobi_c.h"
    PATHS ${GUROBI_INCLUDE_DIR})

  if(GUROBI_HEADER)
    file(STRINGS "${GUROBI_HEADER}" GUROBI_MAJOR REGEX "^#define GRB_VERSION_MAJOR +[0-9]+")
    string(REGEX REPLACE "^#define GRB_VERSION_MAJOR +([0-9]+)" "\\1" MAJOR ${GUROBI_MAJOR})

    file(STRINGS "${GUROBI_HEADER}" GUROBI_MINOR REGEX "^#define GRB_VERSION_MINOR +[0-9]+")
    string(REGEX REPLACE "^#define GRB_VERSION_MINOR +([0-9]+)" "\\1" MINOR ${GUROBI_MINOR})

    file(STRINGS "${GUROBI_HEADER}" GUROBI_TECHNICAL REGEX "^#define GRB_VERSION_TECHNICAL +[0-9]+")
    string(REGEX REPLACE "^#define GRB_VERSION_TECHNICAL +([0-9]+)" "\\1" TECHNICAL ${GUROBI_TECHNICAL})

    set(GUROBI_VERSION "${MAJOR}.${MINOR}.${TECHNICAL}")
  endif()
endif()

set(GUROBI_INCLUDE_DIRS ${GUROBI_INCLUDE_DIR})
set(GUROBI_LIBRARIES ${GUROBI_LIBRARY})

find_package_handle_standard_args(Gurobi
  FOUND_VAR GUROBI_FOUND
  REQUIRED_VARS GUROBI_INCLUDE_DIR GUROBI_LIBRARY
  VERSION_VAR GUROBI_VERSION)
