# SPDX-FileCopyrightText: 2021-2024 The Ikarus Developers mueller@ibb.uni-stuttgart.de
# SPDX-License-Identifier: LGPL-3.0-or-later

# install headers
install(FILES simpleassemblers.inl simpleassemblers.hh dirichletbcenforcement.hh
              assemblermanipulatorfuser.hh interface.hh assemblermanipulatorbuildingblocks.hh
        DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}/ikarus/assembler
)
