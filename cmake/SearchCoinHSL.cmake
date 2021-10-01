# Once done, this will define
#
#  CoinHSL_INCLUDE_DIRS   - where to find umfpack.h, etc.
#  CoinHSL_LIBRARIES      - List of libraries when using Umfpack.
#  CoinHSL_FOUND          - True if CoinHSL found.

find_package(PkgConfig)

if(PKG_CONFIG_FOUND)
  pkg_check_modules(CoinHSL coinhsl)
endif()

find_package(BLAS)
find_package(LAPACK)

if(CoinHSL_FOUND)
  include(SearchMETIS)
  if(METIS_FOUND)
    set(CoinHSL_LIBRARIES ${CoinHSL_LIBRARIES} ${METIS_LIBRARIES})
  else()
    unset(CoinHSL_FOUND)
  endif()

  if(BLAS_FOUND)
    set(CoinHSL_LIBRARIES ${CoinHSL_LIBRARIES} ${BLAS_LIBRARIES})
  else()
    unset(CoinHSL_FOUND)
  endif()

  if(LAPACK_FOUND)
    set(CoinHSL_LIBRARIES ${CoinHSL_LIBRARIES} ${LAPACK_LIBRARIES})
  else()
    unset(CoinHSL_FOUND)
  endif()

  if(UNIX AND NOT APPLE)
    set(CoinHSL_LIBRARIES ${CoinHSL_LIBRARIES} ${CMAKE_DL_LIBS})
  endif()
endif()
