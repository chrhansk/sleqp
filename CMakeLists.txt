cmake_minimum_required(VERSION 3.3)

project(sleqp
  VERSION 0.1.0
  LANGUAGES C CXX)

set(PROJECT_DESCRIPTION "An active set-based NLP solver")

include(CMakePackageConfigHelpers)
include(FindPackageHandleStandardArgs)

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

option(SLEQP_ENABLE_CUTEST "Whether or not to include the CUTest suite" OFF)

option(SLEQP_ENABLE_PYTHON "Whether or not to enable python bindings" OFF)

option(SLEQP_DEBUG "Whether or not to enable debug messages" OFF)

option(SLEQP_FORMAT_CODES "Whether or not to enable ANSI format codes" ON)

set(SLEQP_LPS ""
  CACHE STRING "The LP solver used as a backend")

# Find all packages

include(SearchLPS)

find_package(Umfpack REQUIRED)

find_package(trlib REQUIRED)

find_package(Check)
find_package(Threads)

add_library(trlib INTERFACE IMPORTED)

find_package(Doxygen)

configure_file("${CMAKE_SOURCE_DIR}/src/main/sleqp_defs.h.in"
  "${CMAKE_BINARY_DIR}/sleqp_defs.h")

include_directories("${CMAKE_BINARY_DIR}")

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

set(CMAKE_POSITION_INDEPENDENT_CODE ON)

add_subdirectory(src)

if(SLEQP_ENABLE_PYTHON)
  add_subdirectory(bindings)
endif()

# Miscellaneous installs

set(CONFIG_DESTINATION "lib/cmake/sleqp")

set(CONFIG_FILE "${CMAKE_SOURCE_DIR}/cmake/${PROJECT_NAME}-config.cmake")
set(CONFIG_VERSION_FILE "${CMAKE_CURRENT_BINARY_DIR}/${PROJECT_NAME}-config-version.cmake")

set(TARGET_FILE "${PROJECT_NAME}-targets.cmake")

install(EXPORT sleqp-targets
  FILE ${TARGET_FILE}
  NAMESPACE
  ${PROJECT_NAMESPACE}
  DESTINATION
  ${CONFIG_DESTINATION})

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
)
