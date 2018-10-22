cmake_minimum_required(VERSION 3.3)

project(sleqp
  VERSION 0.1.0
  LANGUAGES C CXX)

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

set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} "${CMAKE_SOURCE_DIR}/cmake/")

find_package(Umfpack REQUIRED)

find_package(SoPlex REQUIRED)

find_package(trlib REQUIRED)

include_directories(${UMFPACK_INCLUDE_DIRS} ${SOPLEX_INCLUDE_DIRS})

find_package(Check)
find_package(Threads)

add_library(trlib INTERFACE IMPORTED)

find_package(Doxygen)

if(DOXYGEN_FOUND)
  set(DOXYGEN_OUTPUT_DIRECTORY "${CMAKE_BINARY_DIR}/doc")
endif()

enable_testing()
set(ENABLE_UNIT_TESTS TRUE)

if(CHECK_FOUND)
  include_directories(${CHECK_INCLUDE_DIR})
else()
  message(STATUS "Could not find Check, unit tests will be disabled")
  set(ENABLE_UNIT_TESTS FALSE)
endif()

set(CMAKE_POSITION_INDEPENDENT_CODE ON)

add_subdirectory(src)
