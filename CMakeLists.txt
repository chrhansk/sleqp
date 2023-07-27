cmake_minimum_required(VERSION 3.14)

project(sleqp
  VERSION 1.0.1
  LANGUAGES C CXX)

set(PROJECT_DESCRIPTION "Active set-based NLP solver")
set(PROJECT_MAINTAINER_NAME "Christoph Hansknecht")
set(PROJECT_MAINTAINER_EMAIL "c.hansknecht@tu-braunschweig.de")

set(PROJECT_LONG_DESCRIPTION
"This package provides SLEQP, an \
open-source solver for large-scale nonlinear continuous \
optimization. SLEQP is available as a callable library with interfaces \
to C, Python, AMPL, and MATLAB / Octave. It uses an active set method based on a \
combination of successive linear programming and equality constrained \
quadratic programming.")

include(FeatureSummary)
include(GNUInstallDirs)
include(CMakePackageConfigHelpers)
include(FindPackageHandleStandardArgs)
include(GenerateExportHeader)
include(CTest)

list(APPEND CMAKE_MODULE_PATH "${PROJECT_SOURCE_DIR}/cmake/")

set(VERSION_FILE "${PROJECT_BINARY_DIR}/sleqp_version")
write_file(${VERSION_FILE} "${PROJECT_VERSION}")

include(CompilerOptions)
include(ExtractDefine)
include(Sanitizers)

# Get the current working branch
execute_process(
  COMMAND git rev-parse --abbrev-ref HEAD
  WORKING_DIRECTORY ${PROJECT_SOURCE_DIR}
  OUTPUT_VARIABLE "PROJECT_GIT_BRANCH"
  OUTPUT_STRIP_TRAILING_WHITESPACE
)

# Get the latest abbreviated commit hash of the working branch
execute_process(
  COMMAND git log -1 --format=%h
  WORKING_DIRECTORY ${PROJECT_SOURCE_DIR}
  OUTPUT_VARIABLE "PROJECT_GIT_COMMIT_HASH"
  OUTPUT_STRIP_TRAILING_WHITESPACE
)

set(PROJECT_NAMESPACE "${PROJECT_NAME}::")
set(LICENSE_FILE_NAME "LICENSE")
set(LICENSE_FILE "${PROJECT_SOURCE_DIR}/${LICENSE_FILE_NAME}")

string(TOUPPER "${PROJECT_NAME}" PROJECT_PREFIX)

# Options / user variables

option(SLEQP_ENABLE_UNIT_TESTS "Whether or not to enable unit tests" ON)

option(SLEQP_ENABLE_C_UNIT_TESTS "Whether or not to enable C unit tests" ON)

option(SLEQP_ENABLE_NUM_ASSERTS "Whether or not to numerical assertions" OFF)

option(SLEQP_ENABLE_CUTEST "Whether or not to include the CUTest suite" OFF)

option(SLEQP_ENABLE_PYTHON "Whether or not to enable python bindings" ON)

option(SLEQP_DEBUG "Whether or not to enable debug messages" OFF)

option(SLEQP_GENERATE_COVERAGE "Whether or not to generate coverage" OFF)

option(SLEQP_FORMAT_CODES "Whether or not to enable ANSI format codes" ON)

option(SLEQP_ENABLE_OCTAVE_MEX "Compile with Octave MEX interface." OFF)
option(SLEQP_ENABLE_MATLAB_MEX "Compile with MATLAB MEX interface." OFF)

option(SLEQP_ENABLE_AMPL "Compile with AMPL interface." OFF)

option(SLEQP_MUMPS_NEEDS_MPI "Does MUMPS require MPI?" ON)

# Required for MUMPS
set(SLEQP_WITH_MPI OFF)

set(SLEQP_LPS ""
  CACHE STRING "The LP solver used as a backend")

# Find all packages

find_package(PkgConfig)

include(SearchLPS)

include(SearchFact)

find_package(trlib REQUIRED)

set(SLEQP_TRLIB_VERSION "${trlib_VERSION}")

find_package(Octave)

if(SLEQP_ENABLE_OCTAVE_MEX AND NOT(Octave_FOUND))
  message(FATAL_ERROR "Could not find Octave")
endif()

find_package(Matlab COMPONENTS MAIN_PROGRAM)

if(SLEQP_ENABLE_MATLAB_MEX AND NOT(Matlab_FOUND))
  message(FATAL_ERROR "Could not find MATLAB")
endif()

if(${SLEQP_ENABLE_MATLAB_MEX} OR ${SLEQP_ENABLE_OCTAVE_MEX})
  set(SLEQP_ENABLE_MEX On)
endif()

pkg_check_modules(COINASL coinasl)

if(SLEQP_ENABLE_AMPL AND NOT(COINASL_FOUND))
  message(FATAL_ERROR "Could not find ASL")
endif()

find_package(CHECK)
find_package(Threads REQUIRED)

if(NOT(OpenMP_FOUND))
  find_package(OpenMP)
endif()

add_library(trlib INTERFACE IMPORTED)

set(SLEQP_DEPENDENCIES
  Threads::Threads
  trlib::trlib
  ${SLEQP_LPS_LIBRARIES}
  ${SLEQP_FACT_LIBRARIES}
  m)

set(SLEQP_LIBRARY_DIRS
  ${SLEQP_FACT_LIBRARY_DIRS}
  ${SLEQP_LPS_LIBRARY_DIRS})

if(UNIX AND NOT APPLE)
  list(APPEND SLEQP_DEPENDENCIES rt)
endif()

add_custom_target(doc)

# Config files

set(PKGCONFIG_FILE "${PROJECT_BINARY_DIR}/sleqp.pc")

configure_file("${PROJECT_SOURCE_DIR}/cmake/sleqp.pc.in"
  "${PKGCONFIG_FILE}"
  @ONLY)

include_directories("${PROJECT_BINARY_DIR}")

# Testing

enable_testing()

set(TEST_DIR "${CMAKE_CURRENT_SOURCE_DIR}/src/test" CACHE INTERNAL "")

if(SLEQP_ENABLE_UNIT_TESTS)
  add_custom_target(build_tests)
  if(CHECK_FOUND)
    include_directories(${CHECK_INCLUDE_DIR})
  else()
    message(STATUS "Could not find Check, C unit tests will be disabled")
    set(SLEQP_ENABLE_C_UNIT_TESTS OFF)
  endif()
else()
  set(SLEQP_ENABLE_C_UNIT_TESTS OFF)
endif()

if(SLEQP_ENABLE_CUTEST)
  set(SLEQP_ENABLE_UNIT_TESTS ON)
endif()

if(SLEQP_GENERATE_COVERAGE)
  include(Coverage)

  if(NOT SLEQP_ENABLE_UNIT_TESTS)
    message(FATAL_ERROR "Cannot generate coverage without unit tests")
  endif()

  enable_coverage()

  add_test_coverage_target(NAME "test_coverage"
    OUTPUT "${PROJECT_BINARY_DIR}/coverage.xml"
    DEPENDS build_tests
    EXCLUDE ${TEST_DIR})
endif()

# Subdirectories

add_subdirectory(src)

add_subdirectory(bindings)

# Documentation

find_package(Doxygen)

if(DOXYGEN_FOUND)
  add_subdirectory(doc)
endif()

# Coverage

# Miscellaneous installs

set(CONFIG_DESTINATION "${CMAKE_INSTALL_LIBDIR}/cmake/sleqp")

set(PKGCONFIG_DESTINATION "${CMAKE_INSTALL_LIBDIR}/pkgconfig")

set(CONFIG_FILE "${PROJECT_SOURCE_DIR}/cmake/${PROJECT_NAME}-config.cmake")
set(CONFIG_VERSION_FILE "${CMAKE_CURRENT_BINARY_DIR}/${PROJECT_NAME}-config-version.cmake")

set(TARGET_FILE "${PROJECT_NAME}-targets.cmake")

install(EXPORT sleqp-targets
  FILE ${TARGET_FILE}
  NAMESPACE
  ${PROJECT_NAMESPACE}
  DESTINATION
  ${CONFIG_DESTINATION}
  COMPONENT "libraries")

write_basic_package_version_file(
  ${CONFIG_VERSION_FILE}
  VERSION ${PROJECT_VERSION}
  COMPATIBILITY SameMajorVersion
)

install(
  FILES
  ${CONFIG_FILE}
  ${CONFIG_VERSION_FILE}
  DESTINATION
  ${CONFIG_DESTINATION}
  COMPONENT "libraries"
)

install(
  FILES
  "${PKGCONFIG_FILE}"
  DESTINATION
  ${PKGCONFIG_DESTINATION}
  COMPONENT "libraries"
)

include(Packaging)

add_feature_info(AMPL SLEQP_ENABLE_AMPL "Interface to AMPL using ASL")
add_feature_info(CUTEst SLEQP_ENABLE_CUTEST "Interface to CUTEst instances")
add_feature_info(Python SLEQP_ENABLE_PYTHON "Python interface")
add_feature_info(MATLAB SLEQP_ENABLE_MATLAB_MEX "MATLAB mex interface")
add_feature_info(Octave SLEQP_ENABLE_OCTAVE_MEX "Octave mex interface")

feature_summary(WHAT ALL)
