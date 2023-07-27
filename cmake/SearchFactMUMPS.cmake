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

if(SLEQP_MUMPS_NEEDS_MPI)
    find_package(MPI)
endif()

find_library(MUMPS_COMMON NAMES mumps_common coinmumps)
find_library(MUMPS_LIBRARY NAMES dmumps coinmumps)

include(FindPackageHandleStandardArgs)

if(SLEQP_MUMPS_NEEDS_MPI)
    find_package_handle_standard_args(MUMPS
      REQUIRED_VARS MUMPS_INCLUDE_DIR MUMPS_LIBRARY MUMPS_COMMON MPI_C_LIBRARIES MPI_C_INCLUDE_DIRS
      VERSION_VAR MUMPS_VERSION)
else()
    find_package_handle_standard_args(MUMPS
      REQUIRED_VARS MUMPS_INCLUDE_DIR MUMPS_LIBRARY MUMPS_COMMON
      VERSION_VAR MUMPS_VERSION)
endif()

set(MUMPS_LIBRARIES
  ${MUMPS_LIBRARY}
  ${MUMPS_COMMON}
  ${MPI_C_LIBRARIES})

set(MUMPS_INCLUDE_DIRS
  ${MUMPS_INCLUDE_DIR}
  ${MPI_C_INCLUDE_DIRS})

mark_as_advanced(
  MUMPS_INCLUDE_DIRS
  MUMPS_LIBRARIES)


if(SLEQP_MUMPS_NEEDS_MPI)
  set(SLEQP_WITH_MPI ON)
endif()
