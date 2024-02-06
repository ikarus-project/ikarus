// SPDX-FileCopyrightText: 2021-2024 The Ikarus Developers mueller@ibb.uni-stuttgart.de
// SPDX-License-Identifier: LGPL-3.0-or-later

/**
 * \file easvariants.hh
 * \brief Definition of the EAS variants.
 * \ingroup  mechanics
 */

#pragma once

#include <ikarus/finiteelements/mechanics/eas/eas2d.hh>
#include <ikarus/finiteelements/mechanics/eas/eas3d.hh>
#include <ikarus/utils/tensorutils.hh>

namespace Ikarus::EAS {
/**
 * \brief EASVariants structure holding variants of different Enhanced Assumed Strains (EAS).
 *
 * \details EASVariants holds the different variants of EAS formulations for
 * both 2D and 3D elements.
 *
 * \tparam GEO The geometry type.
 */
template <typename GEO>
struct Variants
{
  static constexpr int worldDim = GEO::coorddimension;
  static_assert((worldDim == 2) or (worldDim == 3), "EAS variants are only available 2D and 3D elements.");

  /** \brief Variant type holding different EAS formulations for 2D elements. */
  using EAS2D = std::variant<std::monostate, Q1E4<GEO>, Q1E5<GEO>, Q1E7<GEO>>;

  /** \brief Variant type holding different EAS formulations for 3D elements. */
  using EAS3D = std::variant<std::monostate, H1E9<GEO>, H1E21<GEO>>;

  /** \brief Type of the EAS variant depending on the grid dimension. */
  using type = std::conditional_t<worldDim == 2, EAS2D, EAS3D>;
};
} // namespace Ikarus::EAS
