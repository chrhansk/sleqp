set(CPACK_DEBIAN_PACKAGE_DEPENDS "cssrobopec,libqt4-xml,libqt4-network,libqtgui4,treeupdatablereeti")
set(CPACK_PACKAGE_DESCRIPTION_SUMMARY "${PROJECT_DESCRIPTION}")
set(CPACK_PACKAGE_CONTACT "Christoph Hansknecht <c.hansknecht@tu-braunschweig.de>")
set(CPACK_DEBIAN_PACKAGE_SECTION "science")

set(CPACK_DEBIAN_PACKAGE_DEPENDS "libsuitesparse-dev (>= 5.2), trlib (>= 0.2)")

if(SLEQP_LPS_DEPS_DEBIAN)
  set(CPACK_DEBIAN_PACKAGE_DEPENDS "${CPACK_DEBIAN_PACKAGE_DEPENDS}, ${SLEQP_LPS_DEPS_DEBIAN}")
endif()

set(CPACK_GENERATOR "DEB;TBZ2")

include(CPack)
