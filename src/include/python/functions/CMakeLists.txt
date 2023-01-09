set(HEADERS
  discretefunction.hh
  globalbasis.hh
  hierarchicvectorwrapper.hh
  interpolate.hh
  tree.hh
)

install(FILES ${HEADERS} DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}/dune/python/functions)
