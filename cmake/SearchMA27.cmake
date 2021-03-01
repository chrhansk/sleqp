# Once done, this will define
#
#  MA27_LIBRARIES      - List of libraries when using MA27.
#  MA27_FOUND          - True if MA27 found.

include(SearchCoinHSL)

set(MA27_INCLUDE_DIRS "")

if(METIS_FOUND)
  if(CoinHSL_FOUND)
    set(MA27_LIBRARIES ${CoinHSL_LIBRARIES} ${MEITS_LIBRARY})
    set(MA27_FOUND 1)
  else()
    find_library(MA27_LIBRARY ma27 PATHS $ENV{HSLDIR})

    include(SearchMETIS)

    include(FindPackageHandleStandardArgs)

    find_package_handle_standard_args(MA27
      REQUIRED_VARS MA27_LIBRARY METIS_LIBRARIES)

    set(MA27_LIBRARIES ${MA27_LIBRARY} ${METIS_LIBRARIES})
  endif()
endif()
