# SPDX-FileCopyrightText: 2022 The Ikarus Developers mueller@ibb.uni-stuttgart.de
# SPDX-License-Identifier: LGPL-3.0-or-later

set(modules IkarusMacros.cmake AddAutoDiffFlags.cmake AddEigenFlags.cmake AddMatplotppFlags.cmake
            AddSpdlogFlags.cmake
)
install(FILES ${modules} DESTINATION ${DUNE_INSTALL_MODULEDIR})
