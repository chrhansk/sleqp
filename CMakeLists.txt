cmake_minimum_required(VERSION 3.3)

project(sleqp
  VERSION 0.1.0
  LANGUAGES C CXX)

set(PROJECT_DESCRIPTION "Active set-based NLP solver")
set(PROJECT_MAINTAINER_NAME "Christoph Hansknecht")
set(PROJECT_MAINTAINER_EMAIL "c.hansknecht@tu-braunschweig.de")

set(PROJECT_LONG_DESCRIPTION "This package provides sleqp, an
open-source solver for large-scale nonlinear continuous
optimization. sleqp is available as a callable library with interfaces
to C++, C, and python. It uses an active set method based on a
combination of successive linear programming and equality constrained
quadratic programming.")

include(CMakePackageConfigHelpers)
include(FindPackageHandleStandardArgs)
include(GenerateExportHeader)
include(CheckCCompilerFlag)

# Get the current working branch
execute_process(
  COMMAND git rev-parse --abbrev-ref HEAD
  WORKING_DIRECTORY ${CMAKE_SOURCE_DIR}
  OUTPUT_VARIABLE "PROJECT_GIT_BRANCH"
  OUTPUT_STRIP_TRAILING_WHITESPACE
)

# Get the latest abbreviated commit hash of the working branch
execute_process(
  COMMAND git log -1 --format=%h
  WORKING_DIRECTORY ${CMAKE_SOURCE_DIR}
  OUTPUT_VARIABLE "PROJECT_GIT_COMMIT_HASH"
  OUTPUT_STRIP_TRAILING_WHITESPACE
)

set(PROJECT_NAMESPACE "${PROJECT_NAME}::")

string(TOUPPER "${PROJECT_NAME}" PROJECT_PREFIX)

list(APPEND CMAKE_MODULE_PATH "${CMAKE_SOURCE_DIR}/cmake/")

# Options / user variables

option(SLEQP_ENABLE_UNIT_TESTS "Whether or not to enable unit tests" ON)

option(SLEQP_ENABLE_NUM_ASSERTS "Whether or not to numerical assertions" OFF)

option(SLEQP_ENABLE_CUTEST "Whether or not to include the CUTest suite" OFF)

option(SLEQP_ENABLE_PYTHON "Whether or not to enable python bindings" ON)

option(SLEQP_DEBUG "Whether or not to enable debug messages" OFF)

option(SLEQP_FORMAT_CODES "Whether or not to enable ANSI format codes" ON)

set(SLEQP_LPS ""
  CACHE STRING "The LP solver used as a backend")

# Find all packages

include(SearchLPS)

include(SearchFact)

find_package(trlib REQUIRED)

find_package(Check)
find_package(Threads)

add_library(trlib INTERFACE IMPORTED)

set(SLEQP_DEPENDENCIES
  trlib::trlib
  ${SLEQP_LPS_LIBRARIES}
  ${SLEQP_FACT_LIBRARIES}
  m)

if(UNIX AND NOT APPLE)
  list(APPEND SLEQP_DEPENDENCIES rt)
endif()

# Documentation

find_package(Doxygen)

# Config files

configure_file("${CMAKE_SOURCE_DIR}/src/main/sleqp_defs.h.in"
  "${CMAKE_BINARY_DIR}/sleqp_defs.h")

set(PKGCONFIG_FILE "${CMAKE_BINARY_DIR}/sleqp.pc")

configure_file("${CMAKE_SOURCE_DIR}/cmake/sleqp.pc.in"
  "${PKGCONFIG_FILE}"
  @ONLY)

include_directories("${CMAKE_BINARY_DIR}")

# Testing

enable_testing()

if(CHECK_FOUND)
  include_directories(${CHECK_INCLUDE_DIR})
else()
  message(STATUS "Could not find Check, unit tests will be disabled")
  set(SLEQP_ENABLE_UNIT_TESTS OFF)
  set(SLEQP_ENABLE_CUTEST OFF)
endif()

if(SLEQP_ENABLE_CUTEST)
  set(SLEQP_ENABLE_UNIT_TESTS ON)
endif()

# Compiler options

set(CMAKE_POSITION_INDEPENDENT_CODE ON)

check_c_compiler_flag("-fmacro-prefix-map=.=." HAVE_MACRO_PREFIX_MAP)
check_c_compiler_flag("-ffile-prefix-map=.=." HAVE_FILE_PREFIX_MAP)

check_c_source_compiles(
    "
        __attribute__((warn_unused_result)) int f() {return 0;}
        int main(void) {return 0;}
    "
    SLEQP_HAVE_WARN_UNUSED_RESULT
)

# Subdirectories

add_subdirectory(src)

if(SLEQP_ENABLE_PYTHON)
  add_subdirectory(bindings)
endif()

# Miscellaneous installs

set(CONFIG_DESTINATION "lib/cmake/sleqp")

set(PKGCONFIG_DESTINATION "lib/pkgconfig")

set(CONFIG_FILE "${CMAKE_SOURCE_DIR}/cmake/${PROJECT_NAME}-config.cmake")
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
