include(ExternalProject)
ExternalProject_Add(CUTEst_clone
  GIT_REPOSITORY https://github.com/optimizers/cutest-mirror
  CONFIGURE_COMMAND ""
  BUILD_COMMAND ""
  INSTALL_COMMAND ""
  )
ExternalProject_Get_Property(CUTEst_clone SOURCE_DIR)
set(CUTEST_BASE ${SOURCE_DIR} PARENT_SCOPE)
set(CUTEST_BASE ${SOURCE_DIR})

ExternalProject_Add(sifdecode_clone
  GIT_REPOSITORY https://github.com/optimizers/sifdecode-mirror
  CONFIGURE_COMMAND ""
  BUILD_COMMAND ""
  INSTALL_COMMAND ""
  )

ExternalProject_Add(mastsif_clone
  GIT_REPOSITORY https://github.com/optimizers/mastsif-mirror.git
  CONFIGURE_COMMAND ""
  BUILD_COMMAND ""
  INSTALL_COMMAND ""
  )

ExternalProject_Get_Property(sifdecode_clone SOURCE_DIR)
set(SIFDEC_BASE ${SOURCE_DIR})

ExternalProject_Get_Property(mastsif_clone SOURCE_DIR)
set(MASTSIF_BASE ${SOURCE_DIR})

enable_language (Fortran)

add_library(sleqp_cutest
  STATIC
  sleqp_cutest_defs.c
  sleqp_cutest_driver.c
  sleqp_cutest_constrained.c
  sleqp_cutest_unconstrained.c)

target_link_libraries(sleqp_cutest sleqp)

include_directories("${CUTEST_BASE}/include")

file(MAKE_DIRECTORY "${CMAKE_CURRENT_BINARY_DIR}/sifdecoded")

add_custom_command(
 OUTPUT "${CMAKE_CURRENT_BINARY_DIR}/sifdecoded/sifdecoded.f90"
 COMMAND sed -f "${SIFDEC_BASE}/seds/double.sed" "${SIFDEC_BASE}/src/decode/sifdecode.f90" > "${CMAKE_CURRENT_BINARY_DIR}/sifdecoded/sifdecoded.f90"
 DEPENDS sifdecode_clone)

add_custom_command(
 OUTPUT "${CMAKE_CURRENT_BINARY_DIR}/sifdecoded/sifdecode_maind.f90"
 COMMAND sed -f "${SIFDEC_BASE}/seds/double.sed" "${SIFDEC_BASE}/src/decode/sifdecode_main.f90" > "${CMAKE_CURRENT_BINARY_DIR}/sifdecoded/sifdecode_maind.f90"
 DEPENDS sifdecode_clone)

add_executable(sifdecode "${CMAKE_CURRENT_BINARY_DIR}/sifdecoded/sifdecode_maind.f90" "${CMAKE_CURRENT_BINARY_DIR}/sifdecoded/sifdecoded.f90")

file(MAKE_DIRECTORY "${CMAKE_CURRENT_BINARY_DIR}/cutestd")

add_custom_command(
  OUTPUT "${CMAKE_CURRENT_BINARY_DIR}/cutestd/fortran_ops_mod.f90"
  COMMAND sed "s/CHARACTER \\( LEN = 64 \\)/CHARACTER \\( LEN = 255 \\)/" "${CUTEST_BASE}/src/tools/fortran_ops.f90" > "${CMAKE_CURRENT_BINARY_DIR}/cutestd/fortran_ops_mod.f90"
  DEPENDS CUTEst_clone)

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
    COMMAND sed -f "${CUTEST_BASE}/seds/double.sed" "${CUTEST_BASE}/src/tools/${s}.f90" > "${CMAKE_CURRENT_BINARY_DIR}/cutestd/${s}d.f90"
    DEPENDS CUTEst_clone)
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

file(GLOB SIFPROBLEMS_EXT RELATIVE ${MASTSIF_BASE} ${MASTSIF_BASE}/*.SIF DEPENDS mastsif_clone)

set(SIFPROBLEMS "")

add_custom_target(cutest_solvers)
add_custom_target(cutest_libraries)

foreach(s ${SIFPROBLEMS_EXT})
  string(REPLACE ".SIF" "" PROBLEM ${s})

  list(APPEND SIFPROBLEMS "${PROBLEM}")

  file(MAKE_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/PROBLEMS/${PROBLEM})

  add_custom_command(
    OUTPUT "${CMAKE_CURRENT_BINARY_DIR}/PROBLEMS/${PROBLEM}/A.SIF"
    COMMAND cp "${MASTSIF_BASE}/${PROBLEM}.SIF" "${CMAKE_CURRENT_BINARY_DIR}/PROBLEMS/${PROBLEM}/A.SIF"
    DEPENDS mastsif_clone)

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

  target_compile_definitions("${SOLVER_NAME}"
    PRIVATE
    PROBLEM_OUTSDIF_D="${CMAKE_CURRENT_BINARY_DIR}/PROBLEMS/${PROBLEM}/OUTSDIF.d"
    PROBLEM_NAME="${PROBLEM}")

  target_link_libraries("${SOLVER_NAME}"
    sleqp_cutest
    m
    ${CMAKE_DL_LIBS}
    "cutest_${PROBLEM}"
    gfortran
    cutest)

  add_test(NAME "${SOLVER_NAME}" COMMAND "${SOLVER_NAME}")

  add_dependencies(cutest_solvers "${SOLVER_NAME}")

endforeach(s)

add_dependencies(build_tests cutest_solvers)

list(LENGTH SIFPROBLEMS NUM_PROBLEMS)

message(STATUS "Suite contains ${NUM_PROBLEMS} problems")
