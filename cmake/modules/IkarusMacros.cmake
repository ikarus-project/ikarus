# SPDX-FileCopyrightText: 2021-2024 The Ikarus Developers mueller@ibb.uni-stuttgart.de
# SPDX-License-Identifier: LGPL-3.0-or-later

include(CMakeFindDependencyMacro)

find_package(spdlog)
include(AddSpdlogFlags)

find_package(Matplot++)
include(AddMatplotppFlags)

find_package(Eigen3 3.4.90 REQUIRED)
include(AddEigenFlags)

find_package(autodiff REQUIRED)
include(AddAutoDiffFlags)

find_package(Spectra REQUIRED)
include(AddSpectraFlags)


set(MUESLI_LIB_DIR "/muesli/lib/lib")

find_library(MUESLI liblibmuesli_linux_ia64.a.a MUESLI_LIB_DIR REQUIRED)
include(AddMuesliFlags)

