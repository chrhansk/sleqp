set(LP_SOLVERS "")

set(SLEQP_LPS_FOUND FALSE)
set(SLEQP_LPS_INCLUDE_DIRS "")
set(SLEQP_LPS_LIBRARIES "")
set(SLEQP_LPS_SOURCES "")

macro(add_lp_solver)

  cmake_parse_arguments(
    ARGS # prefix of output variables
    "" # list of names of the boolean arguments (only defined ones will be true)
    "NAME" # list of names of mono-valued arguments
    "SOURCES" # list of names of multi-valued arguments (output variables are lists)
    ${ARGN} # arguments of the function to parse, here we take the all original ones
    )

  string(TOUPPER "${ARGS_NAME}" RESULT_NAME)

  set(RESULT_SOURCE_NAME "${RESULT_NAME}_SOURCES")

  set("${RESULT_SOURCE_NAME}" "")

  foreach(SOURCE ${ARGS_SOURCES})
    list(APPEND "${RESULT_SOURCE_NAME}" "lp/${SOURCE}")
  endforeach()

  list(APPEND LP_SOLVERS ${ARGS_NAME})
endmacro()

add_lp_solver(
  NAME "Gurobi"
  SOURCES sleqp_lpi_gurobi.c)

add_lp_solver(
  NAME "SoPlex"
  SOURCES sleqp_lpi_soplex.cc)

set(_SLEQP_LP_SOLVER_VALUES "")

foreach(LP_SOLVER ${LP_SOLVERS})
  string(TOUPPER "${LP_SOLVER}" RESULT_NAME)
  list(APPEND _SLEQP_LP_SOLVER_VALUES "SLEQP_SOLVER_${RESULT_NAME}")
endforeach()

string(REPLACE ";" ", " SLEQP_LP_SOLVER_VALUES "${_SLEQP_LP_SOLVER_VALUES}")

macro(find_lp_solver)
  cmake_parse_arguments(
    ARGS # prefix of output variables
    "REQUIRED" # list of names of the boolean arguments (only defined ones will be true)
    "NAME" # list of names of mono-valued arguments
    "" # list of names of multi-valued arguments (output variables are lists)
    ${ARGN} # arguments of the function to parse, here we take the all original ones
    )

  if(${ARGS_REQUIRED})
    find_package(${ARGS_NAME} REQUIRED)
  else()
    find_package(${ARGS_NAME})
  endif()

  string(TOUPPER "${ARGS_NAME}" RESULT_NAME)

  if("${${RESULT_NAME}_FOUND}")
    set(SLEQP_LPS_FOUND TRUE)
    set(SLEQP_LPS_INCLUDE_DIRS "${${RESULT_NAME}_INCLUDE_DIRS}")
    set(SLEQP_LPS_LIBRARIES "${${RESULT_NAME}_LIBRARIES}")
    set(SLEQP_LPS_SOURCES "${${RESULT_NAME}_SOURCES}")
  endif()

endmacro()

if(SLEQP_LPS)
  find_lp_solver(NAME ${SLEQP_LPS} REQUIRED)
else()
  message(STATUS "Searching for LP solvers")

  foreach(LP_SOLVER ${LP_SOLVERS})

    find_lp_solver(NAME ${LP_SOLVER})

    if(${SLEQP_LPS_FOUND})
      message(STATUS "Using ${LP_SOLVER} as LP solver")
      set(SLEQP_LPS ${LP_SOLVER} CACHE STRING "The LP solver used as a backend" FORCE)
      break()
    endif()

  endforeach()
endif()

if(NOT "${SLEQP_LPS_FOUND}")
  message(FATAL_ERROR "Failed to find LP solver")
endif()

string(TOUPPER "${SLEQP_LPS}" SLEQP_LPS_NAME)
