# SPDX-FileCopyrightText: 2022 The Ikarus Developers mueller@ibb.uni-stuttgart.de
# SPDX-License-Identifier: LGPL-3.0-or-later

# install headers
install(FILES feRequirements.hh feTraits.hh fEparameter.hh physicsHelper.hh
        DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}/ikarus/finiteElements
)

add_subdirectory(feBases)
add_subdirectory(mechanics)
