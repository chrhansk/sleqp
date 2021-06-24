include_directories(${CHECK_INCLUDE_DIR})

include_directories("${CMAKE_SOURCE_DIR}/src/main")
include_directories("${CMAKE_SOURCE_DIR}/src/test")

include_directories(${CHECK_INCLUDE_DIRS})

set(TEST_COMMON_SRC
  test_common.c
  quadcons_fixture.c
  quadfunc_fixture.c
  rosenbrock_fixture.c)

add_library(test_common
  EXCLUDE_FROM_ALL
  ${TEST_COMMON_SRC})

target_link_libraries(test_common
  ${SLEQP_DEPENDENCIES}
  ${CHECK_LIBRARIES}
  m
  sleqp_objects)

add_custom_target(run_tests
  COMMAND ${CMAKE_CTEST_COMMAND}
  DEPENDS build_tests)

function(add_unit_test TEST_NAME)
  get_filename_component(BASE_NAME ${TEST_NAME} NAME)

  add_executable(${BASE_NAME}
    EXCLUDE_FROM_ALL
    "${TEST_NAME}.c")

  target_link_libraries(${BASE_NAME}
    test_common)

  target_compile_options(${BASE_NAME} PRIVATE ${CHECK_CFLAGS})

  add_dependencies(build_tests ${BASE_NAME})

  add_test(NAME ${BASE_NAME} COMMAND ${BASE_NAME} WORKING_DIRECTORY ${CMAKE_BINARY_DIR})
endfunction()

add_unit_test(lp/lpi_test)
add_unit_test(sparse/sleqp_sparse_matrix_test)

add_unit_test(preprocessor/fixed_var_func_test)
add_unit_test(preprocessor/preprocessor_test)

add_unit_test(callback_test)
add_unit_test(cauchy_test)
add_unit_test(constrained_cauchy_test)
add_unit_test(constrained_newton_test)
add_unit_test(constrained_test)
add_unit_test(dual_estimation_test)
add_unit_test(lsq_test)
add_unit_test(mem_test)
add_unit_test(problem_scaling_test)
add_unit_test(scale_test)
add_unit_test(second_order_test)
add_unit_test(unconstrained_newton_test)
add_unit_test(unconstrained_test)
add_unit_test(working_set_var_test)

if(SLEQP_ENABLE_CUTEST)
  add_subdirectory(cutest)
endif()
