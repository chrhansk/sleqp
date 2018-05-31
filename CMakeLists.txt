cmake_minimum_required(VERSION 3.1)

project(sleqp)

set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} "${CMAKE_SOURCE_DIR}/cmake/")

find_package(GMP REQUIRED)
find_package(SoPlex REQUIRED)
find_package(ZLIB REQUIRED)

find_package(Check)

enable_testing()
set(ENABLE_UNIT_TESTS TRUE)

if(NOT CHECK_FOUND)
  message(WARNING "Could not find Check, unit tests will be disabled")
  set(ENABLE_UNIT_TESTS FALSE)
endif()

add_subdirectory(src)
