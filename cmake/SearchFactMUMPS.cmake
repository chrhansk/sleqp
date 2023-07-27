# Once done, this will define
#
#  MUMPS_INCLUDE_DIRS   - where to find mumps_c_types.h, etc.
#  MUMPS_LIBRARIES      - List of libraries when using Umfpack.
#  MUMPS_FOUND          - True if Umfpack found.

find_path(MUMPS_INCLUDE_DIR
  NAMES dmumps_c.h
  PATHS $ENV{MUMPSDIR}
  PATH_SUFFIXES mumps mumps-seq-shared)

if(MUMPS_INCLUDE_DIR)
  define_extract_version(
    "${MUMPS_INCLUDE_DIR}/dmumps_c.h"
    MUMPS_VERSION)
endif()

find_library(MUMPS_COMMON NAMES mumps_common coinmumps)
find_library(MUMPS_LIBRARY NAMES dmumps coinmumps)

include(FindPackageHandleStandardArgs)

set(MUMPS_REQUIRED_VARS
  MUMPS_INCLUDE_DIR
  MUMPS_LIBRARY
  MUMPS_COMMON)

set(MUMPS_LIBRARIES
  ${MUMPS_LIBRARY}
  ${MUMPS_COMMON})

set(MUMPS_INCLUDE_DIRS
  ${MUMPS_INCLUDE_DIR})

# If we need to initialize MPI manually, we need
# access to MPI itself (i.e., link against the library)
if(SLEQP_MUMPS_INIT_MPI)
  find_package(MPI)

  list(APPEND MUMPS_REQUIRED_VARS MPI_C_LIBRARIES MPI_C_INCLUDE_DIRS)
  list(APPEND MUMPS_LIBRARIES ${MPI_C_LIBRARIES})
  list(APPEND MUMPS_INCLUDE_DIRS ${MPI_C_INCLUDE_DIRS})

  set(SLEQP_WITH_MPI On)
endif()

find_package_handle_standard_args(MUMPS
  REQUIRED_VARS ${MUMPS_REQUIRED_VARS}
  VERSION_VAR MUMPS_VERSION)

mark_as_advanced(
  MUMPS_INCLUDE_DIRS
  MUMPS_LIBRARIES)
