set(LP_SOLVERS "")

set(SLEQP_LPS_FOUND FALSE)
set(SLEQP_LPS_INCLUDE_DIRS "")
set(SLEQP_LPS_LIBRARIES "")
set(SLEQP_LPS_SOURCES "")
set(SLEQP_LPS_VERSION "")

set(SLEQP_LPS_DEPS_DEBIAN "")

macro(add_lp_solver)

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
    list(APPEND "${RESULT_SOURCE_NAME}" "lp/${SOURCE}")
  endforeach()

  list(APPEND LP_SOLVERS ${ARGS_NAME})
endmacro()

add_lp_solver(
  NAME "Gurobi"
  SOURCES "lpi_gurobi.c")

add_lp_solver(
  NAME "SoPlex"
  SOURCES "lpi_soplex.cc"
  DEPS_DEBIAN "scipoptsuite (>= 7.0.0)")

set(_SLEQP_LP_SOLVER_VALUES "")

foreach(LP_SOLVER ${LP_SOLVERS})
  string(TOUPPER "${LP_SOLVER}" RESULT_NAME)
  list(APPEND _SLEQP_LP_SOLVER_VALUES "SLEQP_SOLVER_${RESULT_NAME}")
endforeach()

string(REPLACE ";" ", " SLEQP_LP_SOLVER_VALUES "${_SLEQP_LP_SOLVER_VALUES}")

macro(find_lp_solver)
  cmake_parse_arguments(
    ARGS       # prefix of output variables
    "REQUIRED" # list of names of the boolean arguments (only defined ones will be true)
    "NAME"     # list of names of mono-valued arguments
    ""         # list of names of multi-valued arguments (output variables are lists)
    ${ARGN}    # arguments of the function to parse, here we take the all original ones
    )

  message(STATUS "Finding LP solver ${ARGS_NAME}")

  include("Search${ARGS_NAME}")

  string(TOUPPER "${ARGS_NAME}" RESULT_NAME)

  if("${${RESULT_NAME}_FOUND}")
    set(SLEQP_LPS_FOUND TRUE)
    set(SLEQP_LPS_INCLUDE_DIRS "${${RESULT_NAME}_INCLUDE_DIRS}")
    set(SLEQP_LPS_LIBRARIES "${${RESULT_NAME}_LIBRARIES}")
    set(SLEQP_LPS_SOURCES "${${RESULT_NAME}_SOURCES}")
    set(SLEQP_LPS_VERSION "${${RESULT_NAME}_VERSION}")
    set(SLEQP_LPS_DEPS_DEBIAN "${${RESULT_NAME}_DEPS_DEBIAN}")
  elseif(${ARGS_REQUIRED})
    message(FATAL_ERROR "Could not find LP solver ${ARGS_NAME}")
  endif()

endmacro()

if(SLEQP_LPS)
  find_lp_solver(NAME ${SLEQP_LPS} REQUIRED)
else()
  message(STATUS "Searching for LP solvers")

  foreach(LP_SOLVER ${LP_SOLVERS})

    find_lp_solver(NAME ${LP_SOLVER})

    if(${SLEQP_LPS_FOUND})
      set(SLEQP_LPS ${LP_SOLVER} CACHE STRING "The LP solver used as a backend" FORCE)
      break()
    endif()

  endforeach()
endif()

if("${SLEQP_LPS_FOUND}")
  message(STATUS "Using ${SLEQP_LPS} as LP solver")
else()
  message(FATAL_ERROR "Failed to find LP solver")
endif()

string(TOUPPER "${SLEQP_LPS}" SLEQP_LPS_NAME)

set(SLEQP_LPS_PRETTY_NAME "${SLEQP_LPS}")
