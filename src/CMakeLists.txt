add_subdirectory(main)

if(SLEQP_ENABLE_UNIT_TESTS)
  add_subdirectory(test)
endif()
