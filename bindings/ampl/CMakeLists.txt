pkg_check_modules(COINASL coinasl)

if (COINASL_FOUND)

  add_executable(sleqp_ampl
    ampl_main.c
    ampl_data.c
    ampl_output.c
    ampl_problem.c
    ampl_suffix.c
    ampl_util.c)

  link_directories(${COINASL_LIBRARY_DIRS})
  include_directories(${COINASL_INCLUDE_DIRS})

  add_dependencies(sleqp_ampl
    sleqp_local_headers
    sleqp)

  target_link_libraries(sleqp_ampl
    coinasl
    sleqp)

  set_target_properties(sleqp_ampl
    PROPERTIES
    OUTPUT_NAME sleqp)

endif (COINASL_FOUND)
