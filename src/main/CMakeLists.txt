configure_file("${CMAKE_CURRENT_SOURCE_DIR}/sleqp_defs.h.in"
  "${CMAKE_BINARY_DIR}/sleqp_defs.h")

include_directories("${CMAKE_BINARY_DIR}")

set(SOPLEX_LPI_SRC
  lp/sleqp_lpi_soplex.cc)

set(COMMON_SRC
  sparse/sleqp_sparse_factorization.c
  sparse/sleqp_sparse_matrix.c
  sparse/sleqp_sparse_vec.c
  sleqp_active_set.c
  sleqp_aug_jacobian.c
  sleqp_cauchy.c
  sleqp_cmp.c
  sleqp_dual_estimation.c
  sleqp_func.c
  sleqp_iterate.c
  sleqp_params.c
  sleqp_log.c
  sleqp_newton.c
  sleqp_merit.c
  sleqp_problem.c
  sleqp_soc.c
  sleqp_solver.c
  sleqp.c
  lp/sleqp_lpi.c)

add_library(sleqp
  ${SOPLEX_LPI_SRC}
  ${COMMON_SRC})

set_property(TARGET sleqp PROPERTY C_STANDARD 11)

add_executable(sleqp_bin
  sleqp_main.c)

target_link_libraries(sleqp
  trlib::trlib
  ${UMFPACK_LIBRARIES}
  ${SOPLEX_LIBRARIES}
  m)

if(DOXYGEN_FOUND)
  set(DOXYGEN_OUTPUT_DIRECTORY "${CMAKE_BINARY_DIR}/doc")
  set(DOXYGEN_USE_MATHJAX YES)
  doxygen_add_docs(doc ${CMAKE_CURRENT_SOURCE_DIR})
endif()

set_target_properties(sleqp_bin
  PROPERTIES OUTPUT_NAME sleqp)

target_include_directories(sleqp
  PUBLIC
  $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}>
  $<INSTALL_INTERFACE:"include/sleqp">)

target_link_libraries(sleqp_bin sleqp)

install(TARGETS sleqp
  EXPORT sleqp-config
  ARCHIVE DESTINATION lib
  LIBRARY DESTINATION lib)

install(EXPORT sleqp-config
  NAMESPACE "${PROJECT_NAMESPACE}"
  DESTINATION "lib/cmake/sleqp")

install(DIRECTORY "${CMAKE_SOURCE_DIR}/src/main/"
  DESTINATION include/sleqp
  FILES_MATCHING PATTERN "sleqp*.h")

install(FILES "${CMAKE_BINARY_DIR}/sleqp_defs.h"
  DESTINATION include/sleqp)
