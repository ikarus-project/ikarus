# SPDX-FileCopyrightText: 2022 The Ikarus Developers mueller@ibb.uni-stuttgart.de
# SPDX-License-Identifier: LGPL-2.1-or-later

# install headers
install(
  FILES algorithms.hh
        autodiffHelper.hh
        concepts.hh
        duneUtilities.hh
        eigenDuneTransformations.hh
        eigenSparseAddon.hh
        findLineSegment.hh
        functionSanityChecks.hh
        linearAlgebraHelper.hh
        pathFollowingFunctions.hh
        polyfit.hh
        traits.hh
  DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}/ikarus/utils
)

add_subdirectory(drawing)
add_subdirectory(observer)
