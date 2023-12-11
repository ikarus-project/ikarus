// SPDX-FileCopyrightText: 2021-2022 The Ikarus Developers mueller@ibb.uni-stuttgart.de
// SPDX-License-Identifier: LGPL-3.0-or-later
#pragma once

#include <ikarus/utils/makeenum.hh>

namespace Ikarus {
  MAKE_ENUM(StrainTags, linear, deformationGradient, displacementGradient, greenLagrangian, rightCauchyGreenTensor)
  MAKE_ENUM(StressTags, linear, PK2, PK1, Cauchy, Kirchhoff)
  MAKE_ENUM(TangentModuliTags, Material, Spatial, TwoPoint)
}  // namespace Ikarus
