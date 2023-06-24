cmake_minimum_required(VERSION 3.3)

project(sleqp_example
  LANGUAGES C)

add_executable (sleqp_example example.c)

find_package(sleqp REQUIRED)

target_link_libraries(sleqp_example sleqp::sleqp)

enable_testing()

add_test(NAME sleqp_example COMMAND sleqp_example)
