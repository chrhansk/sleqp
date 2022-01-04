set(MEX_SOURCES
  mex_func.c
  mex_func_common.c
  mex_lsq.c
  mex_main.c
  mex_output.c
  mex_problem.c
  mex_solve.c
  mex_solve_common.c
  mex_solve_lsq.c)

set(MATLAB_OUTPUT_DIR "${CMAKE_BINARY_DIR}/matlab")
set(SLEQP_MATLAB_OUTPUT_DIR "${MATLAB_OUTPUT_DIR}/+sleqp")

set(MATLAB_SOURCES
  sleqp/solve.m
  sleqp/solve_lsq.m)

set(MATLAB_TESTS
  tests/test_solve.m
  tests/test_solve_lsq.m
  tests/test_solve_quasi_newton.m
  tests/test_working_set.m)

add_custom_target(sleqp_mex ALL)

add_executable(mex_enum_generator
  gen_enums.c)

target_include_directories(mex_enum_generator
  PRIVATE
  ${SLEQP_INCLUDE_DIR}
  ${SLEQP_LOCAL_HEADER_DIR})

add_dependencies(mex_enum_generator sleqp_local_headers)

target_link_libraries(mex_enum_generator
  PRIVATE
  sleqp_objects)

add_custom_command(TARGET sleqp_mex
    POST_BUILD
    COMMAND mex_enum_generator
    WORKING_DIRECTORY ${SLEQP_MATLAB_OUTPUT_DIR})

foreach(MATLAB_SOURCE ${MATLAB_SOURCES})
  get_filename_component(MATLAB_SOURCE_ABS ${MATLAB_SOURCE} ABSOLUTE)

  add_custom_command(TARGET sleqp_mex
    PRE_BUILD
    COMMAND ${CMAKE_COMMAND} -E
    copy ${MATLAB_SOURCE_ABS} ${SLEQP_MATLAB_OUTPUT_DIR})

endforeach()

if(SLEQP_ENABLE_OCTAVE_MEX)

  add_library(sleqp_ext SHARED
    ${MEX_SOURCES})

  set_target_properties(sleqp_ext
    PROPERTIES
    PREFIX ""
    OUTPUT_NAME "_extension"
    SUFFIX ".mex"
    LIBRARY_OUTPUT_DIRECTORY ${SLEQP_MATLAB_OUTPUT_DIR})

  target_link_libraries(sleqp_ext
    PRIVATE
    sleqp)

  target_include_directories(sleqp_ext
    PRIVATE
    ${Octave_INCLUDE_DIR})

  target_include_directories(sleqp_ext
    PRIVATE
    ${SLEQP_LOCAL_HEADER_DIR})

  add_dependencies(sleqp_ext sleqp_local_headers)

  add_dependencies(sleqp_mex sleqp_ext)

  add_dependencies(build_tests sleqp_mex)

  foreach(MATLAB_TEST ${MATLAB_TESTS})
    get_filename_component(MATLAB_TEST_ABS ${MATLAB_TEST} ABSOLUTE)
    get_filename_component(TEST_NAME ${MATLAB_TEST} NAME_WE)

    add_test(
      NAME "matlab_${TEST_NAME}"
      WORKING_DIRECTORY ${MATLAB_OUTPUT_DIR}
      COMMAND ${Octave_EXECUTABLE} ${MATLAB_TEST_ABS})

  endforeach()

endif()
