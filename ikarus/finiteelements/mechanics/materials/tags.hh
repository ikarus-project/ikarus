// SPDX-FileCopyrightText: 2021-2025 The Ikarus Developers ikarus@ibb.uni-stuttgart.de
// SPDX-License-Identifier: LGPL-3.0-or-later

/**
 * \file tags.hh
 * \brief Definition of several material related enums
 */

#pragma once

#include <ikarus/utils/makeenum.hh>

namespace Ikarus {
/**
 *
 * \ingroup Materialtags
 * \brief A strongly typed enum class representing the type of the passed strain
 */
MAKE_ENUM(StrainTags, linear, displacementGradient, deformationGradient, rightCauchyGreenTensor, greenLagrangian);

/**
 *
 * \ingroup Materialtags
 * \brief A strongly typed enum class representing the type of the computed stresses
 */
MAKE_ENUM(StressTags, linear, PK2, PK1, Cauchy, Kirchhoff);

/**
 *
 * \ingroup Materialtags
 * \brief A strongly typed enum class representing the type of the computed tangent moduli
 */
MAKE_ENUM(TangentModuliTags, Material, Spatial, TwoPoint);

/**
 *
 * \ingroup StretchTag
 * \brief A strongly typed enum class representing the type of the the principal strechts used in hyperelastic materials
 */
MAKE_ENUM(PrincipalStretchTags, total, deviatoric);

} // namespace Ikarus
