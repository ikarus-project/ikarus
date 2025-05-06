# SPDX-FileCopyrightText: 2021-2025 The Ikarus Developers ikarus@ibb.uni-stuttgart.de
# SPDX-License-Identifier: LGPL-3.0-or-later
add_subdirectory(easfunctions)
add_subdirectory(easvariants)
# install headers
install(FILES easvariants.hh easfunctions.hh
        DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}/ikarus/finiteelements/mechanics/strainenhancements
)
