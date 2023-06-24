set(EXAMPLE_FILE example.c)

configure_file(cmake_fetch/CMakeLists.txt.in
  ${CMAKE_CURRENT_SOURCE_DIR}/cmake_fetch/CMakeLists.txt
  @ONLY)

file(COPY ${EXAMPLE_FILE} DESTINATION ${CMAKE_CURRENT_SOURCE_DIR}/make)
file(COPY ${EXAMPLE_FILE} DESTINATION ${CMAKE_CURRENT_SOURCE_DIR}/cmake_install)
file(COPY ${EXAMPLE_FILE} DESTINATION ${CMAKE_CURRENT_SOURCE_DIR}/cmake_fetch)
