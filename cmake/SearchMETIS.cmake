# Once done, this will define
#
#  METIS_LIBRARIES      - List of libraries when using METIS.
#  METIS_FOUND          - True if MA27 found.

find_path(METIS_INCLUDE_DIR
  NAMES metis.h
  PATHS $ENV{METISDIR})

if(METIS_INCLUDE_DIR)
  define_extract_int("${METIS_INCLUDE_DIR}/metis.h"
    METIS_VER_MAJOR)
  define_extract_int("${METIS_INCLUDE_DIR}/metis.h"
    METIS_VER_MINOR)
  define_extract_int("${METIS_INCLUDE_DIR}/metis.h"
    METIS_VER_SUBMINOR)

  set(METIS_VERSION "${METIS_VER_MAJOR}.${METIS_VER_MINOR}.${METIS_VER_SUBMINOR}")
endif()

find_library(METIS_LIBRARY metis PATHS $ENV{METISDIR})

find_package_handle_standard_args(METIS
  REQUIRED_VARS METIS_LIBRARY
  VERSION_VAR METIS_VERSION)

mark_as_advanced(METIS_LIBRARY)

set(METIS_LIBRARIES ${METIS_LIBRARY})
