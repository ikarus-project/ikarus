# SPDX-FileCopyrightText: 2022 The Ikarus Developers mueller@ibb.uni-stuttgart.de
# SPDX-License-Identifier: LGPL-3.0-or-later

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
