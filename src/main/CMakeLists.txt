set(SOPLEX_LPI_SRC
  lp/sleqp_lpi_soplex.cc)

set(COMMON_SRC
  sparse/sleqp_sparse.c
  sleqp_cauchy.c
  sleqp_cmp.c
  sleqp_func.c
  sleqp_iterate.c
  sleqp_log.c
  sleqp.c
  lp/sleqp_lpi.c)

add_library(sleqp
  ${SOPLEX_LPI_SRC}
  ${COMMON_SRC})

add_executable(sleqp_bin
  sleqp_main.c)

target_link_libraries(sleqp
  ${UMFPACK_LIBRARIES}
  ${SOPLEX_LIBRARIES})

set_target_properties(sleqp_bin
  PROPERTIES OUTPUT_NAME sleqp)

include_directories("${CMAKE_CURRENT_SOURCE_DIR}")

target_link_libraries(sleqp_bin sleqp)
