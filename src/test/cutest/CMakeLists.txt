include(FetchContent)

FetchContent_Declare(CUTEst
  GIT_REPOSITORY https://github.com/optimizers/cutest-mirror
  )

FetchContent_Declare(sifdecode
  GIT_REPOSITORY https://github.com/optimizers/sifdecode-mirror
  )

FetchContent_Declare(mastsif
  GIT_REPOSITORY https://github.com/optimizers/mastsif-mirror.git
  )

FetchContent_MakeAvailable(CUTEst sifdecode mastsif)

set(CUTEST_BASE ${cutest_SOURCE_DIR})
set(SIFDECODE_BASE ${sifdecode_SOURCE_DIR})
set(MASTSIF_BASE ${mastsif_SOURCE_DIR})

enable_language(Fortran)

add_library(sleqp_cutest
  STATIC
  sleqp_cutest_types.c
  sleqp_cutest_driver.c
  sleqp_cutest_data.c
  sleqp_cutest_constrained.c
  sleqp_cutest_unconstrained.c)

target_link_libraries(sleqp_cutest sleqp)

include_directories("${CUTEST_BASE}/include")

file(MAKE_DIRECTORY "${CMAKE_CURRENT_BINARY_DIR}/sifdecoded")

add_custom_command(
 OUTPUT "${CMAKE_CURRENT_BINARY_DIR}/sifdecoded/sifdecoded.f90"
 COMMAND sed -f "${SIFDECODE_BASE}/seds/double.sed" "${SIFDECODE_BASE}/src/decode/sifdecode.f90" > "${CMAKE_CURRENT_BINARY_DIR}/sifdecoded/sifdecoded.f90")

add_custom_command(
 OUTPUT "${CMAKE_CURRENT_BINARY_DIR}/sifdecoded/sifdecode_maind.f90"
 COMMAND sed -f "${SIFDECODE_BASE}/seds/double.sed" "${SIFDECODE_BASE}/src/decode/sifdecode_main.f90" > "${CMAKE_CURRENT_BINARY_DIR}/sifdecoded/sifdecode_maind.f90")

add_executable(sifdecode
  "${CMAKE_CURRENT_BINARY_DIR}/sifdecoded/sifdecode_maind.f90"
  "${CMAKE_CURRENT_BINARY_DIR}/sifdecoded/sifdecoded.f90")

file(MAKE_DIRECTORY "${CMAKE_CURRENT_BINARY_DIR}/cutestd")

add_custom_command(
  OUTPUT "${CMAKE_CURRENT_BINARY_DIR}/cutestd/fortran_ops_mod.f90"
  COMMAND sed "s/CHARACTER \\( LEN = 64 \\)/CHARACTER \\( LEN = 255 \\)/" "${CUTEST_BASE}/src/tools/fortran_ops.f90" > "${CMAKE_CURRENT_BINARY_DIR}/cutestd/fortran_ops_mod.f90")

set(CUTESTSRC
  ccfg
  ccfsg
  cchprods
  ccifg
  ccifsg
  cdhc
  cdh
  cdimchp
  cdimen
  cdimse
  cdimsh
  cdimsj
  ceh
  cfn
  cgrdh
  cgr
  chcprod
  chprod
  cidh
  cish
  cjprod
  clfg
  cnames
  cofg
  cofsg
  connames
  creport
  csetup
  csgreh
  csgr
  csgrsh
  cshc
  cshcprod
  csh
  cshp
  cshprod
  csjprod
  cstats
  cterminate
  cutest
  cvartype
  fortran_ops
  interface
  lqp
  newthread
  pname
  problem
  probname
  readin
  timings
  ubandh
  udh
  udimen
  udimse
  udimsh
  ueh
  ufn
  ugrdh
  ugreh
  ugr
  ugrsh
  uhprod
  unames
  uofg
  ureport
  usetup
  ush
  ushp
  ushprod
  uterminate
  uvartype
  varnames)

foreach(s ${CUTESTSRC})
  add_custom_command(
    OUTPUT "${CMAKE_CURRENT_BINARY_DIR}/cutestd/${s}d.f90"
    COMMAND sed -f "${CUTEST_BASE}/seds/double.sed" "${CUTEST_BASE}/src/tools/${s}.f90" > "${CMAKE_CURRENT_BINARY_DIR}/cutestd/${s}d.f90")
  add_custom_command(
    OUTPUT "${CMAKE_CURRENT_BINARY_DIR}/cutestd/${s}sized.f90"
    COMMAND sed -f "${CUTEST_BASE}/seds/cutersize.sed" "${CMAKE_CURRENT_BINARY_DIR}/cutestd/${s}d.f90" > "${CMAKE_CURRENT_BINARY_DIR}/cutestd/${s}sized.f90"
    DEPENDS "${CMAKE_CURRENT_BINARY_DIR}/cutestd/${s}d.f90")
  if(NOT (${s} STREQUAL "readin" OR ${s} STREQUAL "fortran_ops"))
    set(cutestsrc ${cutestsrc} "${CMAKE_CURRENT_BINARY_DIR}/cutestd/${s}sized.f90")
  endif(NOT (${s} STREQUAL "readin" OR ${s} STREQUAL "fortran_ops"))
endforeach(s)

set(cutestsrc ${cutestsrc} "${CMAKE_CURRENT_BINARY_DIR}/cutestd/fortran_ops_mod.f90")
add_library(cutest OBJECT ${cutestsrc})
set_property(TARGET cutest PROPERTY POSITION_INDEPENDENT_CODE TRUE)

file(GLOB SIFPROBLEMS_EXT RELATIVE ${MASTSIF_BASE} ${MASTSIF_BASE}/*.SIF)

set(SIFPROBLEMS "")

find_package(GETOPT REQUIRED)

add_custom_target(cutest_solvers)
add_custom_target(cutest_libraries)

set(WITH_OPENMP ${OpenMP_FOUND})

configure_file(cutest_defs.h.in
  "${CMAKE_CURRENT_BINARY_DIR}/cutest_defs.h"
  @ONLY)

foreach(s ${SIFPROBLEMS_EXT})
  string(REPLACE ".SIF" "" PROBLEM ${s})

  list(APPEND SIFPROBLEMS "${PROBLEM}")

  file(MAKE_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/PROBLEMS/${PROBLEM})

  add_custom_command(
    OUTPUT "${CMAKE_CURRENT_BINARY_DIR}/PROBLEMS/${PROBLEM}/A.SIF"
    COMMAND cp "${MASTSIF_BASE}/${PROBLEM}.SIF" "${CMAKE_CURRENT_BINARY_DIR}/PROBLEMS/${PROBLEM}/A.SIF")

  configure_file("${CMAKE_CURRENT_SOURCE_DIR}/AUTOMAT.d.in"
    "${CMAKE_CURRENT_BINARY_DIR}/PROBLEMS/${PROBLEM}/AUTOMAT.d")

  configure_file("${CMAKE_CURRENT_SOURCE_DIR}/SIFDECODE.CNF.in"
    "${CMAKE_CURRENT_BINARY_DIR}/PROBLEMS/${PROBLEM}/SIFDECODE.CNF")

  add_custom_command(OUTPUT
    "${CMAKE_CURRENT_BINARY_DIR}/PROBLEMS/${PROBLEM}/ELFUN.f"
    "${CMAKE_CURRENT_BINARY_DIR}/PROBLEMS/${PROBLEM}/EXTER.f"
    "${CMAKE_CURRENT_BINARY_DIR}/PROBLEMS/${PROBLEM}/GROUP.f"
    "${CMAKE_CURRENT_BINARY_DIR}/PROBLEMS/${PROBLEM}/OUTSDIF.d"
    "${CMAKE_CURRENT_BINARY_DIR}/PROBLEMS/${PROBLEM}/RANGE.f"
    COMMAND sifdecode
    WORKING_DIRECTORY "PROBLEMS/${PROBLEM}/"
    DEPENDS
    "${CMAKE_CURRENT_BINARY_DIR}/PROBLEMS/${PROBLEM}/A.SIF"
    "${CMAKE_CURRENT_BINARY_DIR}/PROBLEMS/${PROBLEM}/SIFDECODE.CNF")

  add_library(${PROBLEM} OBJECT
    "${CMAKE_CURRENT_BINARY_DIR}/PROBLEMS/${PROBLEM}/ELFUN.f"
    "${CMAKE_CURRENT_BINARY_DIR}/PROBLEMS/${PROBLEM}/EXTER.f"
    "${CMAKE_CURRENT_BINARY_DIR}/PROBLEMS/${PROBLEM}/GROUP.f"
    "${CMAKE_CURRENT_BINARY_DIR}/PROBLEMS/${PROBLEM}/RANGE.f")

  set_property(TARGET ${PROBLEM} PROPERTY POSITION_INDEPENDENT_CODE TRUE)

  add_library("cutest_${PROBLEM}"
    SHARED
    EXCLUDE_FROM_ALL
    $<TARGET_OBJECTS:cutest>
    $<TARGET_OBJECTS:${PROBLEM}>)

  set_target_properties("cutest_${PROBLEM}" PROPERTIES
    ARCHIVE_OUTPUT_DIRECTORY "${CMAKE_CURRENT_BINARY_DIR}/PROBLEMS/${PROBLEM}"
    LIBRARY_OUTPUT_DIRECTORY "${CMAKE_CURRENT_BINARY_DIR}/PROBLEMS/${PROBLEM}")

  add_dependencies(cutest_libraries "cutest_${PROBLEM}")

  set(SOLVER_NAME "cutest_${PROBLEM}_solver")

  add_executable("${SOLVER_NAME}"
    EXCLUDE_FROM_ALL
    sleqp_cutest_main.c)

  target_include_directories("${SOLVER_NAME}"
    PRIVATE "${GETOPT_INCLUDE_DIR}" "${CMAKE_CURRENT_BINARY_DIR}")

  target_compile_definitions("${SOLVER_NAME}"
    PRIVATE
    PROBLEM_OUTSDIF_D="${CMAKE_CURRENT_BINARY_DIR}/PROBLEMS/${PROBLEM}/OUTSDIF.d"
    PROBLEM_NAME="${PROBLEM}")

  target_link_libraries("${SOLVER_NAME}"
    sleqp_objects
    sleqp_cutest
    m
    ${CMAKE_DL_LIBS}
    "cutest_${PROBLEM}"
    gfortran
    cutest)

  if(OpenMP_FOUND)
    target_link_libraries("${SOLVER_NAME}"
      OpenMP::OpenMP_C)
  endif()

  add_test(NAME "${SOLVER_NAME}" COMMAND "${SOLVER_NAME}")

  add_dependencies(cutest_solvers "${SOLVER_NAME}")

endforeach(s)

add_dependencies(build_tests cutest_solvers)

list(LENGTH SIFPROBLEMS NUM_PROBLEMS)

message(STATUS "Suite contains ${NUM_PROBLEMS} problems")
