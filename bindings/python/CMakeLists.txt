list(APPEND CMAKE_MODULE_PATH "${CMAKE_CURRENT_SOURCE_DIR}/cmake")

file(COPY ${LICENSE_FILE}
  DESTINATION ${CMAKE_CURRENT_BINARY_DIR})

file(RENAME
  "${CMAKE_CURRENT_BINARY_DIR}/${LICENSE_FILE_NAME}"
  "${CMAKE_CURRENT_SOURCE_DIR}/LICENSE.txt")

find_package(Python 3.0 COMPONENTS Interpreter)

if(Python_FOUND)
  find_program(TOX "tox")

  include(PythonHelpers)

  add_python_project(
    PROJECT_NAME "sleqp"
    EXTRA_ARGS "SLEQP_LOCAL_BUILD=1"
    EXTRA_INCDIR "${SLEQP_LOCAL_HEADER_DIR}"
    EXTRA_LIBDIR "${SLEQP_LIBRARY_DIR}"
    CONFIG_FILES
    "setup.py.in"
    "tox.ini.in"
    "docs/conf.py.in"
    PROJECT_COMPONENT "python")

  add_dependencies(sleqp_python sleqp_local_headers sleqp)

  if(SLEQP_ENABLE_UNIT_TESTS)
    add_dependencies(build_tests sleqp_python)
  endif()

  add_dependencies(doc sleqp_python_doc)

else()
  message(STATUS "Could not find python, extension will not be built")
endif()
