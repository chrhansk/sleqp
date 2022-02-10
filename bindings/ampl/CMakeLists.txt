add_executable(sleqp_ampl
  ampl_main.c
  ampl_data.c
  ampl_keywords.c
  ampl_output.c
  ampl_problem.c
  ampl_suffix.c
  ampl_util.c)

target_link_directories(sleqp_ampl
  PRIVATE
  ${COINASL_LIBRARY_DIRS})

target_include_directories(sleqp_ampl
  PRIVATE
  ${COINASL_INCLUDE_DIRS}
  ${SLEQP_LOCAL_HEADER_DIR})

add_dependencies(sleqp_ampl
  sleqp_local_headers
  sleqp)

target_link_libraries(sleqp_ampl
  coinasl
  sleqp)

set_target_properties(sleqp_ampl
  PROPERTIES
  OUTPUT_NAME sleqp)

set(TEST_INSTANCES
  instances/constrained.nl)

if(SLEQP_ENABLE_UNIT_TESTS)
  foreach(TEST_INSTANCE ${TEST_INSTANCES})
    get_filename_component(BASE_NAME ${TEST_INSTANCE} NAME_WE)

    add_test(NAME "ampl_${BASE_NAME}"
      COMMAND sleqp_ampl ${TEST_INSTANCE}
      WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR})
  endforeach()

  add_dependencies(build_tests sleqp_ampl)
endif()

install(TARGETS sleqp_ampl
  DESTINATION ${CMAKE_INSTALL_BINDIR}
  COMPONENT "ampl")
