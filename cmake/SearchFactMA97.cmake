# Once done, this will define
#
#  MA97_LIBRARIES      - List of libraries when using MA97.
#  MA97_FOUND          - True if MA97 found.

include(SearchCoinHSL)

set(MA97_INCLUDE_DIRS "")

if(METIS_FOUND)
  if(CoinHSL_FOUND)
    set(MA97_LIBRARIES ${CoinHSL_LIBRARIES} ${MEITS_LIBRARY})
    set(MA97_FOUND 1)
  else()
    find_library(MA97_LIBRARY ma97 PATHS $ENV{HSLDIR})

    include(SearchMETIS)

    include(FindPackageHandleStandardArgs)

    find_package_handle_standard_args(MA97
      REQUIRED_VARS MA97_LIBRARY METIS_LIBRARIES)

    set(MA97_LIBRARIES ${MA97_LIBRARY} ${METIS_LIBRARIES})
  endif()
endif()
