// SPDX-FileCopyrightText: 2021-2024 The Ikarus Developers mueller@ibb.uni-stuttgart.de
// SPDX-License-Identifier: LGPL-3.0-or-later

/**
 * \file resultevaluators.hh
 * \brief Ikarus Result Evaluators for Stress Analysis
 * \ingroup resultevaluators
 *
 */

#pragma once

#include <dune/common/math.hh>

#include <ikarus/finiteelements/ferequirements.hh>

namespace Dune {
// Forward declaration
template <typename ScalarType, int size>
class FieldVector;
} // namespace Dune

namespace Ikarus::ResultEvaluators {

/**
 * \brief Struct for calculating von Mises stress
 * \ingroup resultevaluators
 * \details The VonMises struct provides a function call operator to calculate von Mises stress.
 * \tparam dim dimension of stress state
 */
template <int dim>
requires(dim == 2 or dim == 3)
struct VonMises
{
  /**
   * \brief Calculate the result quantity (von Mises stress)
   * \param resultArray EigenMatrix containing the stress state
   * \param comp component of result (not used here)
   * \return von Mises stress
   */
  double operator()(const auto& resultArray, [[maybe_unused]] const int comp) const {
    if constexpr (dim == 2) {
      const auto s_x  = resultArray(0, 0);
      const auto s_y  = resultArray(1, 0);
      const auto s_xy = resultArray(2, 0);

      return std::sqrt(Dune::power(s_x, 2) + Dune::power(s_y, 2) - s_x * s_y + 3 * Dune::power(s_xy, 2));
    } else {
      const auto s_x  = resultArray(0, 0);
      const auto s_y  = resultArray(1, 0);
      const auto s_z  = resultArray(2, 0);
      const auto s_yz = resultArray(4, 0);
      const auto s_xz = resultArray(5, 0);
      const auto s_xy = resultArray(6, 0);

      return std::sqrt(Dune::power(s_x, 2) + Dune::power(s_y, 2) + Dune::power(s_z, 2) - s_x * s_y - s_x * s_z -
                       s_y * s_z + 3 * (Dune::power(s_xy, 2) + Dune::power(s_xz, 2) + Dune::power(s_yz, 2)));
    }
  }

  /**
   * \brief Get the name of the result type (VonMises)
   * \return String representing the name
   */
  static std::string name() { return "VonMises"; }

  /**
   * \brief Get the number of components in the result (always 1 for VonMises)
   * \return Number of components
   */
  static int ncomps() { return 1; }
};

/**
 * \brief Struct for calculating principal stresses
 * \ingroup resultevaluators
 * \details The PrincipalStress struct provides a function call operator to calculate principal stresses.
 * \remark  Only 2D stresses are supported
 */
struct PrincipalStress
{
  /**
   * \brief Calculate the result quantity (principal stress)
   * \param resultArray EigenMatrix containing the stress state
   * \param comp component of result
   * \return principal stress
   */
  double operator()(const auto& resultArray, const int comp) const {
    const auto s_x  = resultArray(0, 0);
    const auto s_y  = resultArray(1, 0);
    const auto s_xy = resultArray(2, 0);

    auto t1 = (s_x + s_y) / 2;
    auto t2 = std::sqrt(Dune::power((s_x - s_y) / 2, 2) + Dune::power(s_xy, 2));

    return comp == 0 ? t1 + t2 : t1 - t2;
  }

  /**
   * \brief Get the name of the result type (PrincipalStress)
   * \return String representing the name
   */
  static std::string name() { return "PrincipalStress"; }

  /**
   * \brief Get the number of components in the result (always 2 for PrincipalStress)
   * \return Number of components
   */
  static int ncomps() { return 2; }
};

} // namespace Ikarus::ResultEvaluators
