add_subdirectory(install)

include_directories(
  "${PROJECT_SOURCE_DIR}/src/main"
  "${PROJECT_SOURCE_DIR}/src/test"
  "${PROJECT_BINARY_DIR}/sleqp")

set(TEST_COMMON_SRC
  test_common.c
  constrained_fixture.c
  dyn_constrained_fixture.c
  dyn_rosenbrock_fixture.c
  quadcons_fixture.c
  quadfunc_fixture.c
  linear_lsq.c
  log_rosenbrock_fixture.c
  rosenbrock_fixture.c
  rosenbrock_lsq_fixture.c
  wachbieg_fixture.c
  zero_func.c)

add_library(test_common
  ${TEST_COMMON_SRC})

target_include_directories(test_common
    PUBLIC
    ${CHECK_INCLUDE_DIRS}
    ${SLEQP_LOCAL_HEADER_DIR})

target_compile_options(test_common
  PUBLIC
  ${CHECK_CFLAGS})

target_link_libraries(test_common
  PUBLIC
  ${SLEQP_DEP_LIBRARIES}
  ${CHECK_LIBRARIES}
  m
  sleqp_objects)

target_link_directories(test_common
  PUBLIC
  ${SLEQP_DEP_LIBRARY_DIRS}
  ${CHECK_LIBRARY_DIRS})

target_compile_definitions(test_common
  PRIVATE
  ${SLEQP_DEP_DEFINITIONS})

add_custom_target(run_tests
  COMMAND ${CMAKE_CTEST_COMMAND}
  DEPENDS build_tests)

function(add_unit_test TEST_NAME)
  get_filename_component(BASE_NAME ${TEST_NAME} NAME)

  add_executable(${BASE_NAME}
    "${TEST_NAME}.c")

  target_link_libraries(${BASE_NAME}
    test_common)

  target_compile_options(${BASE_NAME} PRIVATE ${CHECK_CFLAGS})

  target_include_directories(${BASE_NAME}
    PRIVATE
    ${CHECK_INCLUDE_DIRS}
    ${SLEQP_LOCAL_HEADER_DIR})

  add_dependencies(build_tests ${BASE_NAME})

  add_dependencies(${BASE_NAME} sleqp sleqp_local_headers)

  add_test(NAME ${BASE_NAME} COMMAND ${BASE_NAME} WORKING_DIRECTORY ${PROJECT_BINARY_DIR})
endfunction()

add_unit_test(lp/lpi_test)
add_unit_test(sparse/sleqp_sparse_matrix_test)

add_unit_test(preprocessor/fixed_var_func_test)
add_unit_test(preprocessor/fixed_var_lsq_func_test)
add_unit_test(preprocessor/fixed_var_dyn_func_test)

add_unit_test(preprocessor/preprocessor_test)

add_unit_test(step/step_rule_test)

add_unit_test(box_constrained_cauchy_test)
add_unit_test(callback_test)
add_unit_test(cauchy_test)
add_unit_test(constrained_cauchy_test)
add_unit_test(constrained_newton_test)
add_unit_test(constrained_test)
add_unit_test(degraded_cons_test)
add_unit_test(dual_estimation_test)
add_unit_test(dyn_test)
add_unit_test(dyn_constrained_test)
add_unit_test(gauss_newton_test)
add_unit_test(log_test)
add_unit_test(lsq_test)
add_unit_test(mem_test)
add_unit_test(polish_test)
add_unit_test(problem_scaling_test)
add_unit_test(restoration_test)
add_unit_test(restoration_solver_test)
add_unit_test(scale_test)
add_unit_test(second_order_test)
add_unit_test(settings_test)
add_unit_test(solver_state_test)
add_unit_test(time_limit_test)
add_unit_test(unconstrained_cauchy_test)
add_unit_test(unconstrained_newton_test)
add_unit_test(unconstrained_test)
add_unit_test(unbounded_test)
add_unit_test(working_set_var_test)
add_unit_test(thread_test)

set(EXAMPLE_PARAM_FILE
  "sleqp_example.opt")

file(COPY ${EXAMPLE_PARAM_FILE}
  DESTINATION ${PROJECT_BINARY_DIR})

if(SLEQP_ENABLE_CUTEST)
  add_subdirectory(cutest)
endif()
