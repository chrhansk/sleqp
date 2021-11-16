pkg_check_modules(COINASL coinasl)

if (COINASL_FOUND)
    add_executable(sleqp_ampl
      sleqp_ampl_main.c
      sleqp_ampl_data.c
      sleqp_ampl_problem.c)

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

