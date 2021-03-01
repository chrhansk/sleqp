# Once done, this will define
#
#  METIS_LIBRARIES      - List of libraries when using METIS.
#  METIS_FOUND          - True if MA27 found.

function(extract_define file name result)
  file(STRINGS "${file}"
    file_result
    REGEX "^#define ${name} +[0-9]+")
  string(REGEX REPLACE "^#define ${name} +([0-9]+).*" "\\1" replace_result ${file_result})
  set(${result} ${replace_result} PARENT_SCOPE)
endfunction()

find_path(METIS_INCLUDE_DIR
  NAMES metis.h
  PATHS $ENV{METISDIR})

if(METIS_INCLUDE_DIR)
  extract_define("${UMFPACK_INCLUDE_DIRS}/metis.h"
    "METIS_VER_MAJOR"
    METIS_VER_MAJOR)
  extract_define("${UMFPACK_INCLUDE_DIRS}/metis.h"
    "METIS_VER_MINOR"
    METIS_VER_MINOR)
  extract_define("${UMFPACK_INCLUDE_DIRS}/metis.h"
    "METIS_VER_SUBMINOR"
    METIS_VER_SUBMINOR)

  set(METIS_VERSION "${METIS_VER_MAJOR}.${METIS_VER_MINOR}.${METIS_VER_SUBMINOR}")
endif()

find_library(METIS_LIBRARY metis PATHS $ENV{METISDIR})

find_package_handle_standard_args(METIS
  REQUIRED_VARS METIS_LIBRARY
  VERSION_VAR METIS_VERSION)

mark_as_advanced(METIS_LIBRARY)

set(METIS_LIBRARIES ${METIS_LIBRARY})
