# - Try to find the CHECK libraries
#  Once done this will define
#
#  CHECK_FOUND - system has check
#  CHECK_INCLUDE_DIRS - the check include directory
#  CHECK_LIBRARIES - check library
#

find_package(PkgConfig)

# Take care about check.pc settings

if(PKG_CONFIG_FOUND)
  pkg_search_module(CHECK check)
endif()

if(NOT CHECK_FOUND)
  find_path(CHECK_INCLUDE_DIR NAMES check.h)
  find_library(CHECK_LIBRARY NAMES check)

  find_package_handle_standard_args(CHECK
    REQUIRED_VARS CHECK_INCLUDE_DIR CHECK_LIBRARY)

  set(CHECK_INCLUDE_DIRS ${CHECK_INCLUDE_DIR})
  set(CHECK_LIBRARIES ${CHECK_LIBRARY})
endif()

# hide advanced variables from CMake GUIs
mark_as_advanced(CHECK_INCLUDE_DIRS CHECK_LIBRARIES)
