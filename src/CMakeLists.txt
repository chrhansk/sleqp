add_subdirectory(main)

if(CHECK_FOUND)
  add_subdirectory(test)
  enable_testing()
else()
  message(WARNING "Could not find check. Unit tests will be disabled")
endif()
