# SPDX-FileCopyrightText: 2022 The Ikarus Developers mueller@ibb.uni-stuttgart.de
# SPDX-License-Identifier: LGPL-2.1-or-later

# install headers
install(
  FILES constant.hh
        binaryExpr.hh
        dotProductExpr.hh
        exprChecks.hh
        linearStrainsExpr.hh
        negateExpr.hh
        normSquaredExpr.hh
        rebind.hh
        scaleExpr.hh
        sumExpr.hh
        unaryExpr.hh
  DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}/ikarus/localFunctions/expressions
)

add_subdirectory(scalarUnaryExpressions)
