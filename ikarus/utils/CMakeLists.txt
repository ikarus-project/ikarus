# SPDX-FileCopyrightText: 2021-2025 The Ikarus Developers ikarus@ibb.uni-stuttgart.de
# SPDX-License-Identifier: LGPL-3.0-or-later

# install headers
install(
  FILES algorithms.hh
        autodiffhelper.hh
        basis.hh
        concepts.hh
        defaultfunctions.hh
        derivativetraits.hh
        differentiablefunction.hh
        dirichletvalues.hh
        eigendunetransformations.hh
        eigensparseaddon.hh
        findlinesegment.hh
        flatprebasis.hh
        functionsanitychecks.hh
        functionhelper.hh
        init.hh
        lambertw.hh
        linearalgebrahelper.hh
        makeenum.hh
        math.hh
        differentiablefunction.hh
        differentiablefunctionfactory.hh
        polyfit.hh
        pythonautodiffdefinitions.hh
        tensorutils.hh
        tensorproductquadrule.hh
        traits.hh
        traversal.hh
  DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}/ikarus/utils
)

add_subdirectory(drawing)
add_subdirectory(listener)
add_subdirectory(broadcaster)
