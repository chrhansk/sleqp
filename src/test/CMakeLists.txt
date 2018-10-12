include_directories(${CHECK_INCLUDE_DIR})

include_directories("${CMAKE_SOURCE_DIR}/src/main")
include_directories("${CMAKE_SOURCE_DIR}/src/test")

include_directories(${CHECK_INCLUDE_DIRS})

add_custom_target(build_tests)

set(TEST_COMMON_SRC
  rosenbrock_fixture.c)

add_library(test_common ${TEST_COMMON_SRC})

add_custom_target(run_tests
  COMMAND ${CMAKE_CTEST_COMMAND}
  DEPENDS build_tests)

function(add_unit_test TEST_NAME)
  get_filename_component(BASE_NAME ${TEST_NAME} NAME)
  add_executable(${BASE_NAME}
    EXCLUDE_FROM_ALL
    "${TEST_NAME}.c")

  target_link_libraries(${BASE_NAME}
    test_common
    ${CHECK_LIBRARIES}
    ${CMAKE_THREAD_LIBS_INIT}
    sleqp)

  add_dependencies(build_tests ${BASE_NAME})

  add_test(NAME ${BASE_NAME} COMMAND ${BASE_NAME})
endfunction()

add_unit_test(sparse/sleqp_sparse_matrix_test)

add_unit_test(lp/soplex_lpi_test)

add_unit_test(sleqp_cauchy_test)
add_unit_test(sleqp_dual_estimation_test)

add_unit_test(unconstrained_test)
