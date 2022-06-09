# Once done, this will define
#
#  MA57_LIBRARIES      - List of libraries when using MA57.
#  MA57_FOUND          - True if MA57 found.

include(SearchCoinHSL)

set(MA57_INCLUDE_DIRS "")

if(METIS_FOUND)
  if(CoinHSL_FOUND)
    set(MA57_LIBRARIES ${CoinHSL_LIBRARIES} ${MEITS_LIBRARY})
    set(MA57_FOUND 1)
  else()
    find_library(MA57_LIBRARY ma57 PATHS $ENV{HSLDIR})

    include(SearchMETIS)

    include(FindPackageHandleStandardArgs)

    find_package_handle_standard_args(MA57
      REQUIRED_VARS MA57_LIBRARY METIS_LIBRARIES)

    set(MA57_LIBRARIES ${MA57_LIBRARY} ${METIS_LIBRARIES})
  endif()
endif()
