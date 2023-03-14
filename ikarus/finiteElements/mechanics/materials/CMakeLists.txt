# SPDX-FileCopyrightText: 2022 The Ikarus Developers mueller@ibb.uni-stuttgart.de
# SPDX-License-Identifier: LGPL-3.0-or-later

# install headers
install(
  FILES interface.hh
        linearElasticity.hh
        neohooke.hh
        strainConversions.hh
        svk.hh
        tags.hh
        vanishingStress.hh
  DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}/ikarus/finiteElements/mechanics/materials
)
