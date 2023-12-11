# SPDX-FileCopyrightText: 2021-2022 The Ikarus Developers mueller@ibb.uni-stuttgart.de
# SPDX-License-Identifier: LGPL-3.0-or-later

# install headers
install(
  FILES algorithms.hh
        autodiffhelper.hh
        basis.hh
        concepts.hh
        defaultfunctions.hh
        duneutilities.hh
        eigendunetransformations.hh
        eigensparseaddon.hh
        findlinesegment.hh
        flatprebasis.hh
        functionsanitychecks.hh
        init.hh
        linearalgebrahelper.hh
        makeenum.hh
        math.hh
        pathfollowingfunctions.hh
        polyfit.hh
        tensorutils.hh
        traits.hh
  DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}/ikarus/utils
)

add_subdirectory(drawing)
add_subdirectory(observer)
