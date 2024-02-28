# SPDX-FileCopyrightText: 2021-2024 The Ikarus Developers mueller@ibb.uni-stuttgart.de
# SPDX-License-Identifier: LGPL-3.0-or-later

# install headers
install(FILES ferequirements.hh fetraits.hh fehelper.hh physicshelper.hh febase.hh
        DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}/ikarus/finiteelements
)

add_subdirectory(autodiff)
add_subdirectory(mechanics)
