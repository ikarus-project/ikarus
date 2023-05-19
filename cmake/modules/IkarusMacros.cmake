# SPDX-FileCopyrightText: 2022 The Ikarus Developers mueller@ibb.uni-stuttgart.de
# SPDX-License-Identifier: LGPL-3.0-or-later

include(CMakeFindDependencyMacro)

find_package(spdlog)
include(AddSpdlogFlags)

find_package(Matplot++)
include(AddMatplotppFlags)

find_package(Eigen3 3.3.9 REQUIRED)
include(AddEigenFlags)

find_package(autodiff REQUIRED)
include(AddAutoDiffFlags)
