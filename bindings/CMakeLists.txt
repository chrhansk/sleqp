if(SLEQP_ENABLE_AMPL)
  add_subdirectory(ampl)
endif()

if(SLEQP_ENABLE_MEX)
  add_subdirectory(mex)
endif()

if(SLEQP_ENABLE_PYTHON)
  add_subdirectory(python)
endif()
