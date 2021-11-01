set(COMMON_SRC
  aug_jac/aug_jac.c
  aug_jac/box_constrained_aug_jac.c
  aug_jac/standard_aug_jac.c
  aug_jac/unconstrained_aug_jac.c
  cauchy/cauchy.c
  cauchy/box_constrained_cauchy.c
  cauchy/standard_cauchy.c
  cauchy/unconstrained_cauchy.c
  sparse/sparse_factorization.c
  sparse/sparse_matrix.c
  sparse/sparse_vec.c
  tr/lsqr.c
  tr/tr_solver.c
  tr/tr_util.c
  tr/trlib_solver.c
  tr/steihaug_solver.c
  preprocessor/fixed_var_func.c
  preprocessor/preprocessing.c
  preprocessor/preprocessing_state.c
  preprocessor/preprocessor.c
  preprocessor/restore.c
  preprocessor/transform.c
  quasi_newton/quasi_newton.c
  quasi_newton/bfgs.c
  quasi_newton/sr1.c
  solver/callback.c
  #solver/cauchy_step.c
  solver/iteration.c
  solver/print.c
  solver/state.c
  solver/step.c
  solver/time.c
  solver/trust_radius.c
  assert.c
  hess_struct.c
  callback_handler.c
  cmp.c
  deriv_check.c
  dual_estimation.c
  dyn.c
  eqp.c
  feas.c
  func.c
  gauss_newton.c
  iterate.c
  params.c
  linesearch.c
  log.c
  lsq.c
  options.c
  parametric.c
  penalty.c
  polish.c
  problem_scaling.c
  restoration.c
  newton.c
  merit.c
  problem.c
  scale.c
  soc.c
  solver.c
  timer.c
  types.c
  trial_point.c
  trial_point/cauchy_step.c
  trial_point/print.c
  trial_point/time.c
  util.c
  working_set.c
  working_step.c
  lp/lpi.c
  step/step_rule.c
  step/step_rule_direct.c
  step/step_rule_minstep.c
  step/step_rule_window.c)


add_library(sleqp_objects OBJECT
  ${COMMON_SRC}
  ${SLEQP_LPS_SOURCES}
  ${SLEQP_FACT_SOURCES})

set_property(TARGET sleqp_objects PROPERTY C_STANDARD 11)

if(HAVE_MACRO_PREFIX_MAP)
  target_compile_options(sleqp_objects PRIVATE "-fmacro-prefix-map=${CMAKE_CURRENT_SOURCE_DIR}=.")
endif()

if(HAVE_FILE_PREFIX_MAP)
  target_compile_options(sleqp_objects PRIVATE "-ffile-prefix-map=${CMAKE_CURRENT_SOURCE_DIR}=.")
endif()

target_link_libraries(sleqp_objects
  PRIVATE
  ${SLEQP_DEPENDENCIES})

set_target_properties(sleqp_objects
  PROPERTIES
  C_VISIBILITY_PRESET hidden
  VISIBILITY_INLINES_HIDDEN ON)

if(DOXYGEN_FOUND)
  set(DOXYGEN_OUTPUT_DIRECTORY "${CMAKE_BINARY_DIR}/doc")
  set(DOXYGEN_USE_MATHJAX YES)
  set(DOXYGEN_FILE_PATTERNS "*.h")
  doxygen_add_docs(c_doc ${CMAKE_CURRENT_SOURCE_DIR})

  add_dependencies(doc c_doc)
endif()

set(SLEQP_INCLUDE_DIR "${CMAKE_CURRENT_SOURCE_DIR}" CACHE STRING "" FORCE)
set(SLEQP_LIBRARY_DIR "${CMAKE_CURRENT_BINARY_DIR}" CACHE STRING "" FORCE)

set(SLEQP_LOCAL_HEADER_DIR "${CMAKE_BINARY_DIR}" CACHE STRING "" FORCE)

target_include_directories(sleqp_objects
  PRIVATE
  ${SLEQP_LOCAL_HEADER_DIR}
  "${SLEQP_LOCAL_HEADER_DIR}/sleqp"
  ${SLEQP_FACT_INCLUDE_DIRS}
  ${SLEQP_LPS_INCLUDE_DIRS})

target_include_directories(sleqp_objects
  PUBLIC
  $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}>
  $<INSTALL_INTERFACE:"include/sleqp">)

add_library(sleqp SHARED
  sleqp.c)

target_link_libraries(sleqp
  PRIVATE
  sleqp_objects
  ${SLEQP_DEPENDENCIES})

set_target_properties(sleqp
    PROPERTIES
    VERSION ${PROJECT_VERSION})

set(SLEQP_EXPORT_FILE_NAME "${CMAKE_BINARY_DIR}/sleqp/export.h")
set(SLEQP_DEFS_FILE_NAME "${CMAKE_BINARY_DIR}/sleqp/defs.h")

# Generate include list

file(GLOB_RECURSE PUB_HEADERS RELATIVE ${CMAKE_CURRENT_SOURCE_DIR} "pub_*.h")

set(HEADER_INCLUDES "")

add_custom_target(sleqp_local_headers)

foreach(PUB_HEADER ${PUB_HEADERS})

  get_filename_component(HEADER_DIR ${PUB_HEADER} DIRECTORY)
  get_filename_component(HEADER_NAME ${PUB_HEADER} NAME)

  install(FILES ${PUB_HEADER}
    DESTINATION "include/sleqp/${HEADER_DIR}"
    COMPONENT "headers")

  set(HEADER_FULL "${CMAKE_CURRENT_SOURCE_DIR}/${PUB_HEADER}")
  set(HEADER_DEST "${CMAKE_BINARY_DIR}/sleqp/${HEADER_DIR}/${HEADER_NAME}")

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
  "${CMAKE_BINARY_DIR}/sleqp.h"
  @ONLY)

configure_file("${CMAKE_SOURCE_DIR}/src/main/defs.h.in"
  "${SLEQP_DEFS_FILE_NAME}")

generate_export_header(sleqp
  EXPORT_FILE_NAME ${SLEQP_EXPORT_FILE_NAME})

install(TARGETS sleqp
  EXPORT sleqp-targets
  ARCHIVE DESTINATION lib
  LIBRARY DESTINATION lib
  COMPONENT "libraries")


install(
  FILES
  "${SLEQP_DEFS_FILE_NAME}"
  "${SLEQP_EXPORT_FILE_NAME}"
  DESTINATION include/sleqp
  COMPONENT "headers")

install(
  FILES
  "${CMAKE_BINARY_DIR}/sleqp.h"
  DESTINATION include
  COMPONENT "headers")
