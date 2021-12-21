include( FindPackageHandleStandardArgs )

find_program(Octave_CONFIG_EXECUTABLE
  NAMES octave-config
  HINTS ${_hint}
  PATHS ${_path}
  PATH_SUFFIXES ${_suff})

if(Octave_CONFIG_EXECUTABLE)

  execute_process(COMMAND ${Octave_CONFIG_EXECUTABLE} -p BINDIR
    OUTPUT_VARIABLE Octave_BINARY_DIR
    ERROR_QUIET
    OUTPUT_STRIP_TRAILING_WHITESPACE)

  execute_process(COMMAND ${Octave_CONFIG_EXECUTABLE} -p OCTINCLUDEDIR
    OUTPUT_VARIABLE Octave_INCLUDE_DIR
    ERROR_QUIET
    OUTPUT_STRIP_TRAILING_WHITESPACE)

  execute_process(COMMAND ${Octave_CONFIG_EXECUTABLE} -p OCTLIBDIR
    OUTPUT_VARIABLE Octave_OCTLIBDIR
    ERROR_QUIET
    OUTPUT_STRIP_TRAILING_WHITESPACE)

  execute_process(COMMAND ${Octave_CONFIG_EXECUTABLE} -p LIBDIR
    OUTPUT_VARIABLE Octave_LIBDIR
    ERROR_QUIET
    OUTPUT_STRIP_TRAILING_WHITESPACE)

  execute_process(COMMAND ${Octave_CONFIG_EXECUTABLE} -p OCTFILEDIR
    OUTPUT_VARIABLE Octave_OCTFILEDIR
    ERROR_QUIET
    OUTPUT_STRIP_TRAILING_WHITESPACE)

endif()

find_program(Octave_EXECUTABLE
  NAMES octave-cli octave
  HINTS ${Octave_BINARY_DIR})

if(Octave_EXECUTABLE)

  execute_process(COMMAND ${Octave_EXECUTABLE} --version
    OUTPUT_VARIABLE Octave_VERSION
    ERROR_QUIET
    OUTPUT_STRIP_TRAILING_WHITESPACE)

  if(Octave_VERSION MATCHES "GNU Octave, version [0-9]+\\.[0-9]+\\.[0-9]+.*")
    string(REGEX REPLACE "GNU Octave, version ([0-9]+)\\.[0-9]+\\.[0-9]+.*" "\\1" Octave_VERSION_MAJOR ${Octave_VERSION})
    string(REGEX REPLACE "GNU Octave, version [0-9]+\\.([0-9]+)\\.[0-9]+.*" "\\1" Octave_VERSION_MINOR ${Octave_VERSION})
    string(REGEX REPLACE "GNU Octave, version [0-9]+\\.[0-9]+\\.([0-9]+).*" "\\1" Octave_VERSION_PATCH ${Octave_VERSION})

    set(Octave_VERSION ${Octave_VERSION_MAJOR}.${Octave_VERSION_MINOR}.${Octave_VERSION_PATCH})
  endif()

  set(Octave_Interpreter_FOUND true)

endif(Octave_EXECUTABLE)

find_program(Octave_MKOCTFILE
  NAMES mkoctfile
  HINTS ${Octave_BINARY_DIR})

find_package_handle_standard_args(Octave
  VERSION_VAR Octave_VERSION
  REQUIRED_VARS Octave_EXECUTABLE Octave_CONFIG_EXECUTABLE Octave_MKOCTFILE)
