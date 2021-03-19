set(FACTORIZATIONS "")

set(SLEQP_FACT_FOUND FALSE)
set(SLEQP_FACT_INCLUDE_DIRS "")
set(SLEQP_FACT_LIBRARIES "")
set(SLEQP_FACT_SOURCES "")

set(SLEQP_FACT_DEPS_DEBIAN "")

macro(add_factorization)

  cmake_parse_arguments(
    ARGS                  # prefix of output variables
    ""                    # list of names of the boolean arguments (only defined ones will be true)
    "NAME"                # list of names of mono-valued arguments
    "SOURCES;DEPS_DEBIAN" # list of names of multi-valued arguments (output variables are lists)
    ${ARGN}               # arguments of the function to parse, here we take the all original ones
    )

  string(TOUPPER "${ARGS_NAME}" RESULT_NAME)

  set(RESULT_SOURCE_NAME "${RESULT_NAME}_SOURCES")

  set("${RESULT_NAME}_DEPS_DEBIAN" "${ARGS_DEPS_DEBIAN}")

  set("${RESULT_SOURCE_NAME}" "")

  foreach(SOURCE ${ARGS_SOURCES})
    list(APPEND "${RESULT_SOURCE_NAME}" "sparse/${SOURCE}")
  endforeach()

  list(APPEND FACTORIZATIONS ${ARGS_NAME})
endmacro()

add_factorization(
  NAME "Umfpack"
  SOURCES sleqp_sparse_factorization_umfpack.c
  DEPS_DEBIAN "libmumps-5.1.2")

add_factorization(
  NAME "MUMPS"
  SOURCES sleqp_sparse_factorization_mumps.c
  DEPS_DEBIAN "libumfpack5 (>= 5.2)")

add_factorization(
  NAME "MA27"
  SOURCES
  sleqp_sparse_factorization_ma27.c
  hsl_matrix.c)

add_factorization(
  NAME "MA57"
  SOURCES
  sleqp_sparse_factorization_ma57.c
  hsl_matrix.c)

add_factorization(
  NAME "MA86"
  SOURCES
  sleqp_sparse_factorization_ma86.c)

add_factorization(
  NAME "MA97"
  SOURCES
  sleqp_sparse_factorization_ma97.c)

set(_SLEQP_FACTORIZATION_VALUES "")

foreach(FACTORIZATION ${FACTORIZATIONS})
  string(TOUPPER "${FACTORIZATION}" RESULT_NAME)
  list(APPEND _SLEQP_FACTORIZATION_VALUES "SLEQP_SOLVER_${RESULT_NAME}")
endforeach()

string(REPLACE ";" ", " SLEQP_FACTORIZATION_VALUES "${_SLEQP_FACTORIZATION_VALUES}")

macro(find_factorization)
  cmake_parse_arguments(
    ARGS       # prefix of output variables
    "REQUIRED" # list of names of the boolean arguments (only defined ones will be true)
    "NAME"     # list of names of mono-valued arguments
    ""         # list of names of multi-valued arguments (output variables are lists)
    ${ARGN}    # arguments of the function to parse, here we take the all original ones
    )

  message(STATUS "Finding factorization library ${ARGS_NAME}")

  include("Search${ARGS_NAME}")

  string(TOUPPER "${ARGS_NAME}" RESULT_NAME)

  if("${${RESULT_NAME}_FOUND}")
    set(SLEQP_FACT_FOUND TRUE)
    set(SLEQP_FACT_INCLUDE_DIRS "${${RESULT_NAME}_INCLUDE_DIRS}")
    set(SLEQP_FACT_LIBRARIES "${${RESULT_NAME}_LIBRARIES}")
    set(SLEQP_FACT_SOURCES "${${RESULT_NAME}_SOURCES}")
    set(SLEQP_FACT_DEPS_DEBIAN "${${RESULT_NAME}_DEPS_DEBIAN}")
  elseif(${ARGS_REQUIRED})
    message(FATAL_ERROR "Could not find factorization library ${ARGS_NAME}")
  endif()

endmacro()

if(SLEQP_FACT)
  find_factorization(NAME ${SLEQP_FACT} REQUIRED)
else()
  message(STATUS "Searching for factorization libraries")

  foreach(FACTORIZATION ${FACTORIZATIONS})

    find_factorization(NAME ${FACTORIZATION})

    if(${SLEQP_FACT_FOUND})
      set(SLEQP_FACT ${FACTORIZATION} CACHE STRING "The factorization library used as a backend" FORCE)
      break()
    endif()

  endforeach()
endif()

if("${SLEQP_FACT_FOUND}")
  message(STATUS "Using ${SLEQP_FACT} as factorization library")
else()
  message(FATAL_ERROR "Failed to find factorization libraries")
endif()

string(TOUPPER "${SLEQP_FACT}" SLEQP_FACT_NAME)
