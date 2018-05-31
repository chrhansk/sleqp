include_directories(${CHECK_INCLUDE_DIR})

include_directories("${CMAKE_SOURCE_DIR}/src/main")
include_directories("${CMAKE_SOURCE_DIR}/src/test")

add_custom_target(build_tests)

add_custom_target(run_tests
  COMMAND ${CMAKE_CTEST_COMMAND}
  DEPENDS build_tests)

function(add_unit_test TEST_NAME)
  get_filename_component(BASE_NAME ${TEST_NAME} NAME)
  add_executable(${BASE_NAME}
    EXCLUDE_FROM_ALL
    "${TEST_NAME}.c")

  target_link_libraries(${BASE_NAME}
    ${CHECK_LIBRARIES}
    sleqp)

  add_dependencies(build_tests ${BASE_NAME})

  add_test(NAME ${BASE_NAME} COMMAND ${BASE_NAME})
endfunction()


add_unit_test(lp/soplex_lpi_test)
