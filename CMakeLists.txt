cmake_minimum_required(VERSION 3.1)

project(sleqp C CXX)

set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} "${CMAKE_SOURCE_DIR}/cmake/")

find_package(Umfpack REQUIRED)

find_package(SoPlex REQUIRED)

include_directories(${UMFPACK_INCLUDE_DIRS} ${SOPLEX_INCLUDE_DIRS})

find_package(Check)
find_package(Threads)

enable_testing()
set(ENABLE_UNIT_TESTS TRUE)

if(CHECK_FOUND)
  include_directories(${CHECK_INCLUDE_DIR})
else()
  message(STATUS "Could not find Check, unit tests will be disabled")
  set(ENABLE_UNIT_TESTS FALSE)
endif()

add_subdirectory(src)
