# SPDX-FileCopyrightText: 2021-2025 The Ikarus Developers ikarus@ibb.uni-stuttgart.de
# SPDX-License-Identifier: LGPL-3.0-or-later

# install headers
install(
  FILES newtonraphson.hh
        trustregion.hh
        newtonraphsonwithscalarsubsidiaryfunction.hh
        nonlinearsolverfactory.hh
        solverinfos.hh
        nonlinearsolverbase.hh
        nonlinearsolverstate.hh
  DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}/ikarus/solver/nonlinearsolver
)
