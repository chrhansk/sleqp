set(COMMON_SRC
  aug_jac/aug_jac.c
  aug_jac/box_constrained_aug_jac.c
  aug_jac/direct_aug_jac.c
  aug_jac/reduced_aug_jac.c
  aug_jac/standard_aug_jac.c
  aug_jac/unconstrained_aug_jac.c
  callback_handler.c
  cauchy/box_constrained_cauchy.c
  cauchy/cauchy.c
  cauchy/standard_cauchy.c
  cauchy/unconstrained_cauchy.c
  cmp.c
  deriv_check.c
  direction.c
  dual_estimation/dual_estimation.c
  dual_estimation/dual_estimation_lp.c
  dual_estimation/dual_estimation_lsq.c
  dual_estimation/dual_estimation_mixed.c
  dyn.c
  enum.c
  eqp.c
  error.c
  fact/fact.c
  fact/fact_qr.c
  feas.c
  func.c
  gauss_newton.c
  hess_struct.c
  iterate.c
  linesearch.c
  log.c
  lp/lpi.c
  lsq.c
  measure.c
  merit.c
  newton.c
  parametric.c
  penalty.c
  polish.c
  preprocessor/fixed_var_func.c
  preprocessor/preprocessing.c
  preprocessor/preprocessing_state.c
  preprocessor/preprocessor.c
  preprocessor/restore.c
  preprocessor/transform.c
  problem.c
  problem_scaling.c
  problem_solver.c
  problem_solver/callback.c
  problem_solver/iteration.c
  problem_solver/print.c
  problem_solver/solve.c
  problem_solver/state.c
  problem_solver/step.c
  problem_solver/trust_radius.c
  quasi_newton/bfgs.c
  quasi_newton/quasi_newton.c
  quasi_newton/sr1.c
  restoration.c
  scale.c
  settings.c
  soc.c
  solver.c
  solver/callback.c
  solver/phase.c
  solver/print.c
  solver/solve.c
  solver/state.c
  sparse/mat.c
  sparse/vec.c
  step/step_rule.c
  step/step_rule_direct.c
  step/step_rule_minstep.c
  step/step_rule_window.c
  timer.c
  tr/lsqr.c
  tr/steihaug_solver.c
  tr/tr_solver.c
  tr/tr_util.c
  tr/trlib_solver.c
  trial_point.c
  trial_point/cauchy_step.c
  trial_point/print.c
  trial_point/time.c
  types.c
  util.c
  working_set.c
  working_step.c)

set(MPI_SRC
  mpi_utils.c)

set(SOURCES
  ${COMMON_SRC}
  ${SLEQP_LPS_SOURCES}
  ${SLEQP_FACT_SOURCES})

if(SLEQP_WITH_MPI)
  set(SOURCES
    ${SOURCES}
    ${MPI_SRC})
endif()

add_library(sleqp_objects OBJECT
  ${SOURCES})

set_property(TARGET sleqp_objects PROPERTY C_STANDARD 11)

if(HAVE_MACRO_PREFIX_MAP)
  target_compile_options(sleqp_objects PRIVATE "-fmacro-prefix-map=${CMAKE_CURRENT_SOURCE_DIR}=.")
endif()

if(HAVE_FILE_PREFIX_MAP)
  target_compile_options(sleqp_objects PRIVATE "-ffile-prefix-map=${CMAKE_CURRENT_SOURCE_DIR}=.")
endif()

target_link_libraries(sleqp_objects
  PRIVATE
  ${SLEQP_DEP_LIBRARIES})

target_compile_definitions(sleqp_objects
  PRIVATE
  ${SLEQP_DEP_DEFINITIONS})

set_target_properties(sleqp_objects
  PROPERTIES
  C_VISIBILITY_PRESET hidden
  VISIBILITY_INLINES_HIDDEN ON)

set(SLEQP_INCLUDE_DIR "${CMAKE_CURRENT_SOURCE_DIR}" CACHE STRING "" FORCE)
set(SLEQP_LIBRARY_DIR "${CMAKE_CURRENT_BINARY_DIR}" CACHE STRING "" FORCE)

set(SLEQP_LOCAL_HEADER_DIR "${PROJECT_BINARY_DIR}" CACHE STRING "" FORCE)

target_include_directories(sleqp_objects
  PRIVATE
  ${SLEQP_LOCAL_HEADER_DIR}
  "${SLEQP_LOCAL_HEADER_DIR}/sleqp"
  ${SLEQP_DEP_INCLUDE_DIRS})

target_include_directories(sleqp_objects
  PUBLIC
  $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}>)

add_library(sleqp SHARED
  sleqp.c)

target_include_directories(sleqp
  PUBLIC
  $<BUILD_INTERFACE:${SLEQP_LOCAL_HEADER_DIR}>
  $<INSTALL_INTERFACE:${CMAKE_INSTALL_INCLUDEDIR}>)

target_link_libraries(sleqp
  PRIVATE
  sleqp_objects
  ${SLEQP_DEP_LIBRARIES})

target_link_directories(sleqp
  PRIVATE
  ${SLEQP_DEP_LIBRARY_DIRS})

set_target_properties(sleqp
  PROPERTIES
  VERSION ${PROJECT_VERSION})

set(SLEQP_EXPORT_FILE_NAME "${PROJECT_BINARY_DIR}/sleqp/export.h")
set(SLEQP_DEFS_FILE_NAME "${PROJECT_BINARY_DIR}/sleqp/defs.h")

add_library(${PROJECT_NAMESPACE}sleqp ALIAS sleqp)

# Generate include list

file(GLOB_RECURSE PUB_HEADERS RELATIVE ${CMAKE_CURRENT_SOURCE_DIR} "pub_*.h")

set(HEADER_INCLUDES "")

add_custom_target(sleqp_local_headers)

foreach(PUB_HEADER ${PUB_HEADERS})

  get_filename_component(HEADER_DIR ${PUB_HEADER} DIRECTORY)
  get_filename_component(HEADER_NAME ${PUB_HEADER} NAME)

  install(FILES ${PUB_HEADER}
    DESTINATION "${CMAKE_INSTALL_INCLUDEDIR}/sleqp/${HEADER_DIR}"
    COMPONENT "headers")

  set(HEADER_FULL "${CMAKE_CURRENT_SOURCE_DIR}/${PUB_HEADER}")
  set(HEADER_DEST "${PROJECT_BINARY_DIR}/sleqp/${HEADER_DIR}/${HEADER_NAME}")

  add_custom_command(
    TARGET sleqp_local_headers
    COMMAND ${CMAKE_COMMAND} -E copy_if_different
    ${HEADER_FULL}
    ${HEADER_DEST}
    DEPENDS ${HEADER_FULL})

  set(HEADER_INCLUDE "#include \"sleqp/${PUB_HEADER}\"")
  list(APPEND HEADER_INCLUDES ${HEADER_INCLUDE})
endforeach()

add_dependencies(sleqp_objects sleqp_local_headers)

string(JOIN "\n" PUB_HEADER_INCLUDES ${HEADER_INCLUDES})

configure_file("${CMAKE_CURRENT_SOURCE_DIR}/sleqp.h.in"
  "${PROJECT_BINARY_DIR}/sleqp.h"
  @ONLY)

configure_file("${PROJECT_SOURCE_DIR}/src/main/defs.h.in"
  "${SLEQP_DEFS_FILE_NAME}")

generate_export_header(sleqp
  EXPORT_FILE_NAME ${SLEQP_EXPORT_FILE_NAME})

# Make sure we are exporting symbols when compiling sleqp_objects
# NOTE: implicitly added `target_compile_definitions(sleqp PRIVATE sleqp_EXPORTS)` is not enough
target_compile_definitions(sleqp_objects PRIVATE sleqp_EXPORTS)

# Install for main library target (runtime destination required for
# windows, "COMPONENT" has to be repeated everywhere)
install(TARGETS sleqp
  EXPORT sleqp-targets
  ARCHIVE DESTINATION ${CMAKE_INSTALL_LIBDIR}
  COMPONENT "libraries"
  LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR}
  COMPONENT "libraries"
  RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR}
  COMPONENT "libraries")


install(
  FILES
  "${SLEQP_DEFS_FILE_NAME}"
  "${SLEQP_EXPORT_FILE_NAME}"
  DESTINATION "${CMAKE_INSTALL_INCLUDEDIR}/sleqp"
  COMPONENT "headers")

install(
  FILES
  "${PROJECT_BINARY_DIR}/sleqp.h"
  DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}
  COMPONENT "headers")
