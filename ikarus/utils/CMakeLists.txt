# SPDX-FileCopyrightText: 2022 The Ikarus Developers mueller@ibb.uni-stuttgart.de
# SPDX-License-Identifier: LGPL-3.0-or-later

# install headers
install(
  FILES algorithms.hh
        autodiffHelper.hh
        basis.hh
        concepts.hh
        defaultFunctions.hh
        duneUtilities.hh
        eigenDuneTransformations.hh
        eigenSparseAddon.hh
        findLineSegment.hh
        flatPreBasis.hh
        functionSanityChecks.hh
        init.hh
        linearAlgebraHelper.hh
        makeEnum.hh
        pathFollowingFunctions.hh
        polyfit.hh
        tensorUtils.hh
        traits.hh
  DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}/ikarus/utils
)

add_subdirectory(drawing)
add_subdirectory(observer)
