add_executable(sleqp_ampl
  ampl_main.c
  ampl_data.c
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
