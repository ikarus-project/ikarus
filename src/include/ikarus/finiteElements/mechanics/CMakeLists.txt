# SPDX-FileCopyrightText: 2022 The Ikarus Developers mueller@ibb.uni-stuttgart.de
# SPDX-License-Identifier: LGPL-2.1-or-later

# install headers
install(FILES nonLinearElasticityFE.hh enhancedAssumedStrains.hh linearElastic.hh
        DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}/ikarus/finiteElements/mechanics
)
