# SPDX-FileCopyrightText: 2021-2024 The Ikarus Developers mueller@ibb.uni-stuttgart.de
# SPDX-License-Identifier: LGPL-3.0-or-later

set(HEADERS material.hh)

install(FILES ${HEADERS}
        DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}/ikarus/python/finiteelements/materials
)
