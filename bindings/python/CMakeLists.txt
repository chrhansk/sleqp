list(APPEND CMAKE_MODULE_PATH "${CMAKE_CURRENT_SOURCE_DIR}/cmake")

find_package(PythonInterp "3")
find_package(PythonLibs "3")

if(PythonInterp_FOUND AND PythonLibs_FOUND)
  find_program(TOX "tox")

  include(PythonHelpers)

  add_python_project(
    PROJECT_NAME "sleqp"
    PROJECT_COMPONENT "python")

  add_dependencies(python_sleqp sleqp)

else()
  message(STATUS "Could not find python, extension will not be built")
endif()
