# SPDX-FileCopyrightText: 2021-2025 The Ikarus Developers ikarus@ibb.uni-stuttgart.de
# SPDX-License-Identifier: LGPL-3.0-or-later
add_subdirectory(strainenhancements)
add_subdirectory(assumedstress)
add_subdirectory(loads)
add_subdirectory(materials)
# install headers
install(
  FILES nonlinearelastic.hh
        assumedstress.hh
        enhancedassumedstrains.hh
        linearelastic.hh
        materials.hh
        kirchhoffloveshell.hh
        membranestrains.hh
        loads.hh
        truss.hh
  DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}/ikarus/finiteelements/mechanics
)
