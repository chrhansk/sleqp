set(FACTS "")

set(SLEQP_FACT_FOUND FALSE)
set(SLEQP_FACT_INCLUDE_DIRS "")
set(SLEQP_FACT_LIBRARIES "")
set(SLEQP_FACT_SOURCES "")
set(SLEQP_FACT_VERSION "")

set(SLEQP_FACT_DEPS_DEBIAN "")

macro(add_fact)

  cmake_parse_arguments(
    ARGS                  # prefix of output variables
    "QR"                  # list of names of the boolean arguments (only defined ones will be true)
    "NAME"                # list of names of mono-valued arguments
    "SOURCES;DEPS_DEBIAN" # list of names of multi-valued arguments (output variables are lists)
    ${ARGN}               # arguments of the function to parse, here we take the all original ones
    )

  string(TOUPPER "${ARGS_NAME}" RESULT_NAME)

  set(RESULT_SOURCES "${RESULT_NAME}_SOURCES")

  set("${RESULT_NAME}_DEPS_DEBIAN" "${ARGS_DEPS_DEBIAN}")

  set("${RESULT_NAME}_QR" "${ARGS_QR}")

  set("${RESULT_SOURCES}" "${ARGS_SOURCES}")

  list(APPEND FACTS ${ARGS_NAME})
endmacro()

add_fact(
  NAME "Umfpack"
  SOURCES fact/fact_umfpack.c
  DEPS_DEBIAN "libumfpack5 (>= 5.2)")

add_fact(
  NAME "SPQR"
  QR
  SOURCES
  fact/cholmod_helpers.c
  fact/fact_spqr.c
  DEPS_DEBIAN "libspqr2")

add_fact(
  NAME "CHOLMOD"
  SOURCES
  fact/cholmod_helpers.c
  fact/fact_cholmod.c)

add_fact(
  NAME "MUMPS"
  SOURCES fact/fact_mumps.c
  DEPS_DEBIAN "libmumps-5.1.2")

add_fact(
  NAME "MA27"
  SOURCES
  fact/fact_ma27.c
  fact/hsl_matrix.c)

add_fact(
  NAME "MA57"
  SOURCES
  fact/fact_ma57.c
  fact/hsl_matrix.c)

add_fact(
  NAME "MA86"
  SOURCES
  fact/fact_ma86.c)

add_fact(
  NAME "MA97"
  SOURCES
  fact/fact_ma97.c)

add_fact(
  NAME "LAPACK"
  SOURCES
  fact/fact_lapack.c)

set(_SLEQP_FACT_VALUES "")

foreach(FACT ${FACTS})
  string(TOUPPER "${FACT}" RESULT_NAME)
  list(APPEND _SLEQP_FACT_VALUES "SLEQP_SOLVER_${RESULT_NAME}")
endforeach()

string(REPLACE ";" ", " SLEQP_FACT_VALUES "${_SLEQP_FACT_VALUES}")

macro(find_fact)
  cmake_parse_arguments(
    ARGS       # prefix of output variables
    "REQUIRED" # list of names of the boolean arguments (only defined ones will be true)
    "NAME"     # list of names of mono-valued arguments
    ""         # list of names of multi-valued arguments (output variables are lists)
    ${ARGN}    # arguments of the function to parse, here we take the all original ones
    )

  message(STATUS "Finding factorization library ${ARGS_NAME}")

  include("SearchFact${ARGS_NAME}")

  string(TOUPPER "${ARGS_NAME}" RESULT_NAME)

  if("${${RESULT_NAME}_FOUND}")
    set(SLEQP_FACT_FOUND TRUE)
    set(SLEQP_FACT_INCLUDE_DIRS "${${RESULT_NAME}_INCLUDE_DIRS}")
    set(SLEQP_FACT_LIBRARIES "${${RESULT_NAME}_LIBRARIES}")
    set(SLEQP_FACT_SOURCES "${${RESULT_NAME}_SOURCES}")
    set(SLEQP_FACT_VERSION "${${RESULT_NAME}_VERSION}")
    set(SLEQP_FACT_DEPS_DEBIAN "${${RESULT_NAME}_DEPS_DEBIAN}")
  elseif(${ARGS_REQUIRED})
    message(FATAL_ERROR "Could not find factorization library ${ARGS_NAME}")
  endif()

endmacro()

if(SLEQP_FACT)
  find_fact(NAME ${SLEQP_FACT} REQUIRED)
else()
  message(STATUS "Searching for factorization libraries")

  foreach(FACT ${FACTS})

    find_fact(NAME ${FACT})

    if(${SLEQP_FACT_FOUND})
      set(SLEQP_FACT ${FACT} CACHE STRING "The factorization library used as a backend" FORCE)
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

set(SLEQP_FACT_PRETTY_NAME "${SLEQP_FACT}")

if("${${SLEQP_FACT}_QR}")
  message(STATUS "Factorization is QR")
  set(SLEQP_HAVE_QR_FACT On)
endif()
