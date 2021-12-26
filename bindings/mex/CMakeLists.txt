set(MEX_SOURCES
  mex_func.c
  mex_main.c
  mex_output.c
  mex_solve.c)

set(MATLAB_OUTPUT_DIR "${CMAKE_BINARY_DIR}/matlab")
set(SLEQP_MATLAB_OUTPUT_DIR "${MATLAB_OUTPUT_DIR}/+sleqp")

set(MATLAB_SOURCES
  sleqp/solve.m)

set(MATLAB_TESTS
  tests/test_solve.m)

add_custom_target(sleqp_mex ALL)

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

  foreach(MATLAB_TEST ${MATLAB_TESTS})
    get_filename_component(MATLAB_TEST_ABS ${MATLAB_TEST} ABSOLUTE)
    get_filename_component(TEST_NAME ${MATLAB_TEST} NAME_WE)

    add_test(
      NAME "matlab_${TEST_NAME}"
      WORKING_DIRECTORY ${MATLAB_OUTPUT_DIR}
      COMMAND ${Octave_EXECUTABLE} ${MATLAB_TEST_ABS})

  endforeach()

endif()
