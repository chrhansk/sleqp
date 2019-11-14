list(APPEND CMAKE_MODULE_PATH "${CMAKE_CURRENT_SOURCE_DIR}/cmake")

find_program(TOX "tox")

include(PythonHelpers)

add_python_project(
  PROJECT_NAME sleqp
  PROJECT_COMPONENT "python")
