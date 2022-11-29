# SPDX-FileCopyrightText: 2022 The Ikarus Developers mueller@ibb.uni-stuttgart.de
# SPDX-License-Identifier: LGPL-2.1-or-later

get_filename_component(IKARUS_CMAKE_DIR "IkarusMacros.cmake" PATH)
include(CMakeFindDependencyMacro)

find_dependency(spdlog)
find_dependency(Matplot++)
find_dependency(Eigen3)

set(IKARUS_lIBRARIES ikarus)
