cmake_minimum_required(VERSION 3.1)

project(sleqp C CXX)

set(CMAKE_POSITION_INDEPENDENT_CODE ON)

set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} "${CMAKE_SOURCE_DIR}/cmake/")

find_package(GMP REQUIRED)
find_package(SoPlex REQUIRED)
find_package(ZLIB REQUIRED)

find_package(Check)
find_package(Threads)

enable_testing()
set(ENABLE_UNIT_TESTS TRUE)

if(CHECK_FOUND)
  include_directories(${CHECK_INCLUDE_DIR})
else()
  message(WARNING "Could not find Check, unit tests will be disabled")
  set(ENABLE_UNIT_TESTS FALSE)
endif()

add_subdirectory(src)
