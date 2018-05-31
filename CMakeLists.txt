cmake_minimum_required(VERSION 3.1)

project(sleqp)

#set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${CMAKE_ENABLE_C99}")

set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} "${CMAKE_SOURCE_DIR}/cmake/")

find_package(GMP REQUIRED)
find_package(SoPlex REQUIRED)
find_package(ZLIB REQUIRED)

find_package(Check)

add_subdirectory(src)
