set(COMMON_SRC
  sparse/sleqp_sparse_factorization.c
  sparse/sleqp_sparse_matrix.c
  sparse/sleqp_sparse_vec.c
  sleqp_active_set.c
  sleqp_aug_jacobian.c
  sleqp_bfgs.c
  sleqp_cauchy.c
  sleqp_cmp.c
  sleqp_deriv_check.c
  sleqp_dual_estimation.c
  sleqp_func.c
  sleqp_iterate.c
  sleqp_params.c
  sleqp_log.c
  sleqp_newton.c
  sleqp_merit.c
  sleqp_problem.c
  sleqp_scale.c
  sleqp_sr1.c
  sleqp_soc.c
  sleqp_solver.c
  sleqp_timer.c
  sleqp_util.c
  lp/sleqp_lpi.c)

add_library(sleqp SHARED
  ${COMMON_SRC}
  ${SLEQP_LPS_SOURCES})

set_property(TARGET sleqp PROPERTY C_STANDARD 11)

add_executable(sleqp_bin
  sleqp_main.c)

target_link_libraries(sleqp
  trlib::trlib
  ${SLEQP_LPS_LIBRARIES}
  ${UMFPACK_LIBRARIES}
  m)

if(DOXYGEN_FOUND)
  set(DOXYGEN_OUTPUT_DIRECTORY "${CMAKE_BINARY_DIR}/doc")
  set(DOXYGEN_USE_MATHJAX YES)
  doxygen_add_docs(doc ${CMAKE_CURRENT_SOURCE_DIR})
endif()

set_target_properties(sleqp_bin
  PROPERTIES OUTPUT_NAME sleqp)

include_directories(sleqp
  ${UMFPACK_INCLUDE_DIRS}
  ${SLEQP_LPS_INCLUDE_DIRS})

target_include_directories(sleqp
  PUBLIC
  $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}>
  $<INSTALL_INTERFACE:"include/sleqp">)

target_link_libraries(sleqp_bin sleqp)

install(TARGETS sleqp
  EXPORT sleqp-targets
  ARCHIVE DESTINATION lib
  LIBRARY DESTINATION lib)

install(DIRECTORY "${CMAKE_SOURCE_DIR}/src/main/"
  DESTINATION include/sleqp
  FILES_MATCHING PATTERN "sleqp*.h")

install(FILES "${CMAKE_BINARY_DIR}/sleqp_defs.h"
  DESTINATION include/sleqp)
