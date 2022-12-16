# SPDX-FileCopyrightText: 2022 The Ikarus Developers mueller@ibb.uni-stuttgart.de
# SPDX-License-Identifier: LGPL-3.0-or-later

# install headers
install(FILES nonLinearOperator.hh truncatedConjugateGradient.hh dirichletValues.hh
        DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}/ikarus/linearAlgebra
)
