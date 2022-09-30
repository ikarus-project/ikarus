# install headers
install(
  FILES controlLogger.hh
        controlVTKWriter.hh
        genericControlObserver.hh
        gridDrawerObserver.hh
        loadControlObserver.hh
        nonLinearSolverLogger.hh
        observer.hh
        observerMessages.hh
  DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}/ikarus/utils/observer
)
