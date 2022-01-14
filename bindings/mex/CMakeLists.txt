set(MEX_SOURCES
  mex_func.c
  mex_func_common.c
  mex_info.c
  mex_lsq.c
  mex_main.c
  mex_output.c
  mex_problem.c
  mex_solve.c
  mex_solve_common.c
  mex_solve_lsq.c)

set(MATLAB_OUTPUT_DIR "${CMAKE_BINARY_DIR}/matlab")
set(SLEQP_MATLAB_OUTPUT_DIR "${MATLAB_OUTPUT_DIR}/+sleqp")
set(MATLAB_TEST_DIR "${CMAKE_CURRENT_SOURCE_DIR}/tests")

set(MATLAB_SOURCES
  sleqp/info.m
  sleqp/solve.m
  sleqp/solve_lsq.m)

set(MATLAB_TESTS
  tests/test_callback.m
  tests/test_hess_prod.m
  tests/test_solve.m
  tests/test_solve_lsq.m
  tests/test_solve_quasi_newton.m
  tests/test_working_set.m)

add_custom_target(sleqp_mex ALL)

add_custom_command(TARGET sleqp_mex
  PRE_BUILD COMMAND ${CMAKE_COMMAND} -E make_directory ${SLEQP_MATLAB_OUTPUT_DIR})

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
  add_library(sleqp_octave_ext SHARED
    ${MEX_SOURCES})

  set_target_properties(sleqp_octave_ext
    PROPERTIES
    PREFIX ""
    OUTPUT_NAME "extension"
    SUFFIX ".mex"
    LIBRARY_OUTPUT_DIRECTORY ${SLEQP_MATLAB_OUTPUT_DIR})

  target_link_libraries(sleqp_octave_ext
    PRIVATE
    sleqp)

  target_include_directories(sleqp_octave_ext
    PRIVATE
    ${Octave_INCLUDE_DIR})

  target_include_directories(sleqp_octave_ext
    PRIVATE
    ${SLEQP_LOCAL_HEADER_DIR})

  add_dependencies(sleqp_octave_ext sleqp_local_headers)

  add_dependencies(sleqp_mex sleqp_octave_ext)

  add_dependencies(build_tests sleqp_mex)

  foreach(MATLAB_TEST ${MATLAB_TESTS})
    get_filename_component(TEST_NAME ${MATLAB_TEST} NAME_WE)

    add_test(
      NAME "octave_${TEST_NAME}"
      WORKING_DIRECTORY ${MATLAB_TEST_DIR}
      COMMAND ${CMAKE_COMMAND} -E env OCTAVE_PATH=${MATLAB_OUTPUT_DIR}
      ${Octave_EXECUTABLE} --eval ${TEST_NAME})
  endforeach()
endif()

if(SLEQP_ENABLE_MATLAB_MEX)
  matlab_add_mex(
    NAME sleqp_matlab_ext
    SRC ${MEX_SOURCES}
    OUTPUT_NAME "extension"
    LINK_TO sleqp)

  set_target_properties(sleqp_matlab_ext
    PROPERTIES
    LIBRARY_OUTPUT_DIRECTORY ${SLEQP_MATLAB_OUTPUT_DIR})

  foreach(MATLAB_TEST ${MATLAB_TESTS})
    get_filename_component(TEST_NAME ${MATLAB_TEST} NAME_WE)
    add_test(
      NAME "matlab_${TEST_NAME}"
      WORKING_DIRECTORY ${MATLAB_TEST_DIR}
      COMMAND ${CMAKE_COMMAND} -E env MATLABPATH=${MATLAB_OUTPUT_DIR}
      ${Matlab_MAIN_PROGRAM} -nojvm -batch ${TEST_NAME})
  endforeach()

  add_dependencies(sleqp_matlab_ext sleqp_local_headers)

  add_dependencies(sleqp_mex sleqp_matlab_ext)

  add_dependencies(build_tests sleqp_mex)
endif()
