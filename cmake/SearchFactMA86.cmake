# Once done, this will define
#
#  MA86_LIBRARIES      - List of libraries when using MA86.
#  MA86_FOUND          - True if MA86 found.

include(SearchCoinHSL)

set(MA86_INCLUDE_DIRS "")

if(METIS_FOUND)
  if(CoinHSL_FOUND)
    set(MA86_LIBRARIES ${CoinHSL_LIBRARIES} ${MEITS_LIBRARY})
    set(MA86_FOUND 1)
  else()
    find_library(MA86_LIBRARY ma86 PATHS $ENV{HSLDIR})

    include(SearchMETIS)

    include(FindPackageHandleStandardArgs)

    find_package_handle_standard_args(MA86
      REQUIRED_VARS MA86_LIBRARY METIS_LIBRARIES)

    set(MA86_LIBRARIES ${MA86_LIBRARY} ${METIS_LIBRARIES})
  endif()
endif()
