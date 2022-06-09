include(CheckCCompilerFlag)

if(NOT CMAKE_BUILD_TYPE STREQUAL "Debug")
    message(WARNING "Code coverage results with an optimised (non-Debug) build may be misleading")
endif()

find_program( GCOVR gcovr )

if(NOT GCOVR)
  message(FATAL_ERROR "gcovr not found! Aborting...")
endif()

set(COVERAGE_FLAGS "--coverage")

set(CMAKE_REQUIRED_FLAGS ${COVERAGE_FLAGS})
set(CMAKE_REQUIRED_LIBRARIES ${COVERAGE_FLAGS})

check_c_source_compiles("int main(void) { return 0; }" COMPILER_HAS_GCOV)

if(COMPILER_HAS_GCOV)
  message(STATUS "Generating coverage information")
else()
  message(FATAL_ERROR "Compiler does not support coverage generation")
endif()

function(enable_coverage)
  add_compile_options(${COVERAGE_FLAGS})
  add_link_options(${COVERAGE_FLAGS})
endfunction()

function(add_test_coverage_target)

  set(options "")
  set(oneValueArgs BASE_DIRECTORY NAME OUTPUT)
  set(multiValueArgs EXCLUDE DEPENDS)

  cmake_parse_arguments(Coverage "${options}" "${oneValueArgs}" "${multiValueArgs}" ${ARGN})

  if(DEFINED Coverage_BASE_DIRECTORY)
    get_filename_component(BASEDIR ${Coverage_BASE_DIRECTORY} ABSOLUTE)
  else()
    set(BASEDIR ${PROJECT_SOURCE_DIR})
  endif()

  set(GCOVR_EXCLUDE_ARGS "")

  foreach(EXCLUDE ${Coverage_EXCLUDE})
    list(APPEND GCOVR_EXCLUDE_ARGS "--exclude")
    list(APPEND GCOVR_EXCLUDE_ARGS "${EXCLUDE}")
  endforeach()

  set(GCOVR_OUTPUT_ARGS "")

  if(Coverage_OUTPUT)
    list(APPEND GCOVR_OUTPUT_ARGS "--output")
    list(APPEND GCOVR_OUTPUT_ARGS "${Coverage_OUTPUT}")
  endif()

  message(STATUS "${GCOVR} --print-summary --root ${BASEDIR} ${GCOVR_EXCLUDE_ARGS} ${GCOVR_OUTPUT_ARGS}")

  add_custom_target(${Coverage_NAME}
    COMMAND ctest
    COMMAND ${GCOVR} --xml --print-summary --root ${BASEDIR} ${GCOVR_EXCLUDE_ARGS} ${GCOVR_OUTPUT_ARGS}
    WORKING_DIRECTORY ${PROJECT_BINARY_DIR}
    DEPENDS ${Coverage_DEPENDS})
endfunction()
