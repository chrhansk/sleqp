find_package(PythonInterp 3.0 REQUIRED)
find_package(PythonLibs 3.0 REQUIRED)

execute_process(COMMAND
  ${PYTHON_EXECUTABLE} -c "import sys, numpy; sys.stdout.write(numpy.get_include())"
  OUTPUT_VARIABLE NUMPY_INCLUDE_DIRS
  RESULT_VARIABLE NUMPY_NOT_FOUND)

execute_process(COMMAND
  ${PYTHON_EXECUTABLE} -c "import sys,sysconfig; sys.stdout.write(sysconfig.get_config_var('EXT_SUFFIX'))"
  OUTPUT_VARIABLE PYTHON_EXT_SUFFIX
  RESULT_VARIABLE EXT_SUFFIX_NOT_FOUND)

if(NUMPY_NOT_FOUND)
  message(FATAL_ERROR "NumPy headers not found")
endif()

if(EXT_SUFFIX_NOT_FOUND)
  message(FATAL_ERROR "Extension suffix not found")
endif()

find_program(CYTHON "cython")

if(NOT CYTHON)
  message(FATAL_ERROR "Cython executable not found")
endif()

set(CYTHON_SRC "${CMAKE_CURRENT_SOURCE_DIR}/src/sleqp_python.pyx")
set(CYTHON_OUTPUT "${CMAKE_CURRENT_BINARY_DIR}/sleqp_python.c")

set(PYTHON_MOD_NAME "sleqp")

set(CYTHON_DEPS
  src/csleqp.pxd
  src/call.pxi
  src/cmp.pxi
  src/func.pxi
  src/params.pxi
  src/problem.pxi
  src/solver.pxi
  src/sparse.pxi
  src/types.pxi
  src/version.pxi
  )

add_custom_command(OUTPUT ${CYTHON_OUTPUT}
  COMMAND ${CYTHON} -o ${CYTHON_OUTPUT} ${CYTHON_SRC}
  MAIN_DEPENDENCY ${CYTHON_SRC}
  DEPENDS ${CYTHON_DEPS})

add_library(sleqp_python SHARED ${CYTHON_OUTPUT})

set(PYTHON_LIB_NAME "sleqp_python${PYTHON_EXT_SUFFIX}")

set_target_properties(
  sleqp_python
  PROPERTIES
  PREFIX ""
  SUFFIX ""
  LIBRARY_OUTPUT_DIRECTORY "${CMAKE_CURRENT_SOURCE_DIR}/sleqp"
  OUTPUT_NAME "${PYTHON_LIB_NAME}"
  LINKER_LANGUAGE C
  )

target_include_directories(sleqp_python
  PRIVATE
  ${PYTHON_INCLUDE_DIRS}
  ${NUMPY_INCLUDE_DIRS})

target_link_libraries(sleqp_python sleqp)

set(PYTHON_LIB_OUTPUT "${CMAKE_CURRENT_BINARY_DIR}/${PYTHON_LIB_NAME}")

configure_file("${CMAKE_CURRENT_SOURCE_DIR}/setup.py.in"
  "${CMAKE_CURRENT_SOURCE_DIR}/setup.py")


add_test(NAME python_tests
  COMMAND ${PYTHON_EXECUTABLE} "setup.py" "test"
  WORKING_DIRECTORY "${CMAKE_CURRENT_SOURCE_DIR}")

add_dependencies(build_tests sleqp_python)
