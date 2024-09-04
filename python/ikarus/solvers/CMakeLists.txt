# SPDX-FileCopyrightText: 2021-2024 The Ikarus Developers mueller@ibb.uni-stuttgart.de
# SPDX-License-Identifier: LGPL-3.0-or-later

add_python_targets(solvers __init__)

install(FILES registersolverinfo.hh
        DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}/python/ikarus/solvers
)