find_package(PythonInterp 3.0 REQUIRED)
find_package(PythonLibs 3.0 REQUIRED)

execute_process(COMMAND
  ${PYTHON_EXECUTABLE} -c "import sys, numpy; sys.stdout.write(numpy.get_include())"
  OUTPUT_VARIABLE NUMPY_INCLUDE_DIRS
  RESULT_VARIABLE NUMPY_NOT_FOUND)

if(NUMPY_NOT_FOUND)
  message(FATAL_ERROR "NumPy headers not found")
endif()

find_program(CYTHON "cython")

if(NOT CYTHON)
  message(FATAL_ERROR "Cython executable not found")
endif()

set(CYTHON_SRC "${CMAKE_CURRENT_SOURCE_DIR}/src/sleqp.pyx")
set(CYTHON_OUTPUT "${CMAKE_CURRENT_BINARY_DIR}/sleqp.c")

file(GLOB_RECURSE CYTHON_DEPS
  "src/*.pxi"
  "src/*.pxd"
  "src/*.pyx")

add_custom_command(OUTPUT ${CYTHON_OUTPUT}
  COMMAND ${CYTHON} -o ${CYTHON_OUTPUT} ${CYTHON_SRC}
  MAIN_DEPENDENCY ${CYTHON_SRC}
  DEPENDS ${CYTHON_DEPS})

add_library(sleqp_python SHARED ${CYTHON_OUTPUT})

set_target_properties(
  sleqp_python
  PROPERTIES
  PREFIX ""
  OUTPUT_NAME "sleqp"
  LINKER_LANGUAGE C
  )

target_include_directories(sleqp_python
  PRIVATE
  ${PYTHON_INCLUDE_DIRS}
  ${NUMPY_INCLUDE_DIRS})

target_link_libraries(sleqp_python sleqp)
