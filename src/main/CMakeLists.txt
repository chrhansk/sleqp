set(COMMON_SRC
  sparse/sparse_factorization.c
  sparse/sparse_matrix.c
  sparse/sparse_vec.c
  tr/tr_solver.c
  tr/trlib_solver.c
  tr/steihaug_solver.c
  preprocessor/fixed_var_func.c
  preprocessor/preprocessing.c
  preprocessor/preprocessing_state.c
  preprocessor/preprocessor.c
  preprocessor/restore.c
  preprocessor/transform.c
  assert.c
  aug_jacobian.c
  bfgs.c
  hess_struct.c
  callback_handler.c
  cauchy.c
  cmp.c
  deriv_check.c
  dual_estimation.c
  feas.c
  func.c
  iterate.c
  params.c
  linesearch.c
  log.c
  lsq.c
  options.c
  parametric.c
  penalty.c
  problem_scaling.c
  newton.c
  merit.c
  problem.c
  scale.c
  sr1.c
  soc.c
  solver.c
  timer.c
  types.c
  util.c
  working_set.c
  lp/lpi.c)


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

target_include_directories(sleqp_objects
  PRIVATE
  ${CMAKE_BINARY_DIR}
  ${CMAKE_BINARY_DIR}/sleqp
  ${SLEQP_FACT_INCLUDE_DIRS}
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

foreach(PUB_HEADER ${PUB_HEADERS})

  get_filename_component(HEADER_DIR ${PUB_HEADER} DIRECTORY)
  get_filename_component(HEADER_NAME ${PUB_HEADER} NAME)

  install(FILES ${PUB_HEADER}
    DESTINATION "include/sleqp/${HEADER_DIR}"
    COMPONENT "headers")

  add_custom_command(TARGET sleqp POST_BUILD
    COMMAND ${CMAKE_COMMAND} -E copy
    "${CMAKE_CURRENT_SOURCE_DIR}/${PUB_HEADER}"
    "${CMAKE_BINARY_DIR}/sleqp/${HEADER_DIR}/${HEADER_NAME}")

  set(HEADER_INCLUDE "#include \"sleqp/${PUB_HEADER}\"")
  list(APPEND HEADER_INCLUDES ${HEADER_INCLUDE})
endforeach()

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
