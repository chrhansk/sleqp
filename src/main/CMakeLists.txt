set(COMMON_SRC
  sparse/sleqp_sparse_factorization.c
  sparse/sleqp_sparse_factorization_umfpack.c
  sparse/sleqp_sparse_matrix.c
  sparse/sleqp_sparse_vec.c
  sleqp_assert.c
  sleqp_aug_jacobian.c
  sleqp_bfgs.c
  sleqp_hess_struct.c
  sleqp_cauchy.c
  sleqp_cmp.c
  sleqp_deriv_check.c
  sleqp_dual_estimation.c
  sleqp_feas.c
  sleqp_func.c
  sleqp_iterate.c
  sleqp_params.c
  sleqp_linesearch.c
  sleqp_log.c
  sleqp_lsq.c
  sleqp_options.c
  sleqp_problem_scaling.c
  sleqp_newton.c
  sleqp_merit.c
  sleqp_problem.c
  sleqp_scale.c
  sleqp_sr1.c
  sleqp_soc.c
  sleqp_solver.c
  sleqp_timer.c
  sleqp_tr_solver.c
  sleqp_types.c
  sleqp_util.c
  sleqp_working_set.c
  lp/sleqp_lpi.c)


add_library(sleqp_objects OBJECT
  ${COMMON_SRC}
  ${SLEQP_LPS_SOURCES})

set_property(TARGET sleqp_objects PROPERTY C_STANDARD 11)

target_link_libraries(sleqp_objects
  PRIVATE
  trlib::trlib
  ${SLEQP_LPS_LIBRARIES}
  ${UMFPACK_LIBRARIES}
  m)

if(UNIX AND NOT APPLE)
  target_link_libraries(sleqp_objects
    PRIVATE
    rt)
endif()

set_target_properties(sleqp_objects
  PROPERTIES
  C_VISIBILITY_PRESET hidden
  VISIBILITY_INLINES_HIDDEN ON)

if(DOXYGEN_FOUND)
  set(DOXYGEN_OUTPUT_DIRECTORY "${CMAKE_BINARY_DIR}/doc")
  set(DOXYGEN_USE_MATHJAX YES)
  set(DOXYGEN_FILE_PATTERNS "*.h")
  doxygen_add_docs(doc ${CMAKE_CURRENT_SOURCE_DIR})
endif()

target_include_directories(sleqp_objects
  PRIVATE
  ${CMAKE_BINARY_DIR}
  ${UMFPACK_INCLUDE_DIRS}
  ${SLEQP_LPS_INCLUDE_DIRS})

set(SLEQP_INCLUDE_DIR "${CMAKE_CURRENT_SOURCE_DIR}" CACHE STRING "" FORCE)
set(SLEQP_LIBRARY_DIR "${CMAKE_CURRENT_BINARY_DIR}" CACHE STRING "" FORCE)

target_include_directories(sleqp_objects
  PUBLIC
  $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}>
  $<INSTALL_INTERFACE:"include/sleqp">)

add_library(sleqp SHARED
  sleqp.c)

target_link_libraries(sleqp
  PRIVATE
  sleqp_objects)

set_target_properties(sleqp
    PROPERTIES
    VERSION ${PROJECT_VERSION})

set(SLEQP_EXPORT_FILE_NAME "${CMAKE_BINARY_DIR}/sleqp_export.h")

generate_export_header(sleqp
  EXPORT_FILE_NAME SLEQP_EXPORT_FILE_NAME)

install(TARGETS sleqp
  EXPORT sleqp-targets
  ARCHIVE DESTINATION lib
  LIBRARY DESTINATION lib
  COMPONENT "libraries")

install(DIRECTORY "${CMAKE_SOURCE_DIR}/src/main/"
  DESTINATION include/sleqp
  COMPONENT "headers"
  FILES_MATCHING PATTERN "sleqp*.h")

install(
  FILES
  "${CMAKE_BINARY_DIR}/sleqp_defs.h"
  "${SLEQP_EXPORT_FILE_NAME}"
  DESTINATION include/sleqp
  COMPONENT "headers")
