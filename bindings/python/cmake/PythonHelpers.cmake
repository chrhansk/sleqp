find_package(PythonInterp "3" REQUIRED)
find_package(PythonLibs "3" REQUIRED)

set(ENABLE_UNIT_TESTS On)

if(NOT TOX)
  message(WARNING "Could not find tox executable, python unit tests will be disabled")
  set(ENABLE_UNIT_TESTS Off)
endif()

function(add_python_project)

  cmake_parse_arguments(
    PARSE_ARGV 0 ARGS
    ""
    "PROJECT_NAME;PROJECT_COMPONENT"
    ""
    )

  set(PROJECT_NAME "${ARGS_PROJECT_NAME}")

  set(TARGET_NAME "python_${PROJECT_NAME}")

  set(PROJECT_DIR "${CMAKE_CURRENT_SOURCE_DIR}")

  set(SETUP_PY_IN "${CMAKE_CURRENT_SOURCE_DIR}/setup.py.in")
  set(SETUP_PY    "${CMAKE_CURRENT_SOURCE_DIR}/setup.py")

  set(TOX_INI_IN "${CMAKE_CURRENT_SOURCE_DIR}/tox.ini.in")
  set(TOX_INI    "${CMAKE_CURRENT_SOURCE_DIR}/tox.ini")

  set(PYTHON_CFLAGS "")

  if(CMAKE_C_FLAGS)
    set(PYTHON_CFLAGS "${CMAKE_C_FLAGS}")
  elseif(CMAKE_BUILD_TYPE)
    set(PYTHON_CFLAGS "${CMAKE_C_FLAGS_${CMAKE_BUILD_TYPE}}")
  endif()

  configure_file(${SETUP_PY_IN}
    ${SETUP_PY}
    @ONLY)

  configure_file(${TOX_INI_IN}
    ${TOX_INI}
    @ONLY)

  add_custom_target(${TARGET_NAME} ALL)

  add_custom_command(
    TARGET ${TARGET_NAME}
    COMMAND ${CMAKE_COMMAND} -E env CFLAGS=${PYTHON_CFLAGS} ${PYTHON_EXECUTABLE} ${SETUP_PY} build build_ext --inplace
    WORKING_DIRECTORY ${PROJECT_DIR})

  add_custom_target("${TARGET_NAME}_sdist")

  add_custom_command(
    TARGET "${TARGET_NAME}_sdist"
    COMMAND ${CMAKE_COMMAND} -E env CFLAGS=${PYTHON_CFLAGS} ${PYTHON_EXECUTABLE} ${SETUP_PY} sdist --formats=gztar
    WORKING_DIRECTORY ${PROJECT_DIR})

  add_custom_target("${TARGET_NAME}_bdist")

  add_custom_command(
    TARGET "${TARGET_NAME}_bdist"
    COMMAND ${CMAKE_COMMAND} -E env CFLAGS=${PYTHON_CFLAGS} ${PYTHON_EXECUTABLE} ${SETUP_PY} bdist
    WORKING_DIRECTORY ${PROJECT_DIR})

  add_custom_target("${TARGET_NAME}_bdist_wheel")

  add_custom_command(
    TARGET "${TARGET_NAME}_bdist_wheel"
    COMMAND ${CMAKE_COMMAND} -E env CFLAGS=${PYTHON_CFLAGS} ${PYTHON_EXECUTABLE} ${SETUP_PY} bdist_wheel
    WORKING_DIRECTORY ${PROJECT_DIR})

  if(ENABLE_UNIT_TESTS)
    add_test(NAME "${TARGET_NAME}_tests"
      COMMAND ${CMAKE_COMMAND} -E env CFLAGS=${PYTHON_CFLAGS} "${TOX}"
      WORKING_DIRECTORY "${PROJECT_DIR}")
  endif()

  set(INSTALL_SCRIPT_NAME "${CMAKE_CURRENT_BINARY_DIR}/python_install_${PROJECT_NAME}.cmake")

  configure_file("${CMAKE_CURRENT_SOURCE_DIR}/cmake/python_install.cmake.in"
    "${INSTALL_SCRIPT_NAME}"
    @ONLY)

  install(
    SCRIPT "${INSTALL_SCRIPT_NAME}"
    COMPONENT "${ARGS_PROJECT_COMPONENT}")

endfunction()
