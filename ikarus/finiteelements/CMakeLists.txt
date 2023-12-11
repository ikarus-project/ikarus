# SPDX-FileCopyrightText: 2021-2023 The Ikarus Developers mueller@ibb.uni-stuttgart.de
# SPDX-License-Identifier: LGPL-3.0-or-later

# install headers
install(FILES ferequirements.hh fetraits.hh feparameter.hh physicshelper.hh
        DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}/ikarus/finiteelements
)

add_subdirectory(febases)
add_subdirectory(mechanics)
