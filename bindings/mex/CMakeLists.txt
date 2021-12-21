set( MEX_SOURCES
  sleqp.c)

set(MEX_FILE "${CMAKE_CURRENT_BINARY_DIR}/sleqp.mex")
set(M_FILE "${CMAKE_CURRENT_BINARY_DIR}/sleqp.m")

if(SLEQP_ENABLE_OCTAVE_MEX)

  add_library(sleqp_mex SHARED
    ${MEX_SOURCES})

  set_target_properties(sleqp_mex
    PROPERTIES
    PREFIX ""
    OUTPUT_NAME "sleqp"
    SUFFIX ".mex")

  target_link_libraries(sleqp_mex
    PRIVATE
    sleqp)

  target_include_directories(sleqp_mex
    PRIVATE
    ${Octave_INCLUDE_DIR})

  install(FILES ${MEX_FILE} DESTINATION ${Octave_OCTFILEDIR})
  install(FILES ${M_FILE} DESTINATION ${Octave_OCTFILEDIR})
endif()
