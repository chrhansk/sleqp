find_path(GETOPT_INCLUDE_DIR getopt.h)

set(GETOPT_INCLUDE_DIRS ${GETOPT_INCLUDE_DIR})

include(FindPackageHandleStandardArgs)

find_package_handle_standard_args(GETOPT DEFAULT_MSG GETOPT_INCLUDE_DIR)

mark_as_advanced(GETOPT_INCLUDE_DIR)
