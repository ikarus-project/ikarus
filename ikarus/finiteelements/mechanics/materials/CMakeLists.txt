# SPDX-FileCopyrightText: 2021-2025 The Ikarus Developers ikarus@ibb.uni-stuttgart.de
# SPDX-License-Identifier: LGPL-3.0-or-later

# install headers
install(
  FILES interface.hh
        linearelasticity.hh
        strainconversions.hh
        stressconversions.hh
        svk.hh
        tags.hh
        vanishingstress.hh
        vanishingstrain.hh
        materialhelpers.hh
        numericalmaterialinversion.hh
  DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}/ikarus/finiteelements/mechanics/materials
)
add_subdirectory(hyperelastic)
add_subdirectory(muesli)
