set(SOPLEX_LPI_SRC
  lp/sleqp_lpi_soplex.cc)

set(COMMON_SRC
  sparse/sleqp_sparse.c
  sleqp_cmp.c
  lp/sleqp_lpi.c)

add_library(sleqp STATIC
  ${SOPLEX_LPI_SRC}
  ${COMMON_SRC})

add_executable(sleqp_bin
  sleqp_main.c)

target_link_libraries(sleqp
  ${SOPLEX_LIBRARIES}
  ${GMP_LIBRARIES}
  ${ZLIB_LIBRARIES})

set_target_properties(sleqp_bin
  PROPERTIES OUTPUT_NAME sleqp)

include_directories("${CMAKE_CURRENT_SOURCE_DIR}")

target_link_libraries(sleqp_bin sleqp)
