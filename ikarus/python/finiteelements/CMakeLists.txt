# SPDX-FileCopyrightText: 2021-2025 The Ikarus Developers ikarus@ibb.uni-stuttgart.de
# SPDX-License-Identifier: LGPL-3.0-or-later

set(HEADERS fe.hh registerpreelement.hh material.hh registerferequirements.hh registerpreelement.hh
            scalarwrapper.hh
)

install(FILES ${HEADERS} DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}/ikarus/python/finiteelements)
