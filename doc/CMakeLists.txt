set(DOXYGEN_OUTPUT_DIRECTORY "${CMAKE_BINARY_DIR}/doc")
set(DOXYGEN_USE_MATHJAX YES)
set(DOXYGEN_FILE_PATTERNS "pub*.h")
# doxygen_add_docs(c_doc ${CMAKE_CURRENT_SOURCE_DIR})

set(DOXYGEN_PROJECT_VERSION ${PROJECT_VERSION})
string(TOUPPER ${PROJECT_NAME} DOXYGEN_PROJECT_NAME)
set(DOXYGEN_INPUT ${SLEQP_INCLUDE_DIR})

set(DOXYGEN_FORMULA_MACROFILE "${CMAKE_CURRENT_SOURCE_DIR}/macros.tex")

configure_file(Doxyfile.in
  Doxyfile
  @ONLY)

add_custom_target(doxygen_doc)

add_custom_command(TARGET doxygen_doc
  POST_BUILD
  COMMAND ${DOXYGEN}
  WORKING_DIRECTORY ${DOXYGEN_OUTPUT_DIRECTORY})

if(MAKE)
  add_custom_target(doxygen_latex)

  add_custom_command(
    TARGET doxygen_latex
    COMMAND ${MAKE}
    WORKING_DIRECTORY "${DOXYGEN_OUTPUT_DIRECTORY}/latex")

  add_dependencies(doc doxygen_latex)
  add_dependencies(doxygen_latex doxygen_doc)
endif()

add_dependencies(doc doxygen_doc)
