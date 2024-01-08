// SPDX-FileCopyrightText: 2021-2024 The Ikarus Developers mueller@ibb.uni-stuttgart.de
// SPDX-License-Identifier: LGPL-3.0-or-later

/**
 * @file resultevaluators.hh
 * @brief Ikarus Result Evaluators for Stress Analysis
 * @ingroup resultevaluators
 *
 */

#pragma once

#include <dune/common/math.hh>

#include <ikarus/finiteelements/ferequirements.hh>

namespace Dune {
  // Forward declaration
  template <typename ScalarType, int size>
  class FieldVector;
}  // namespace Dune

namespace Ikarus::ResultEvaluators {

  /**
   * @brief Struct for calculating von Mises stress
   * @ingroup resultevaluators
   * @details The VonMises struct provides a function call operator to calculate von Mises stress.
   * @remark  Only 2D stresses are supported
   * @tparam ElementType Type of the finite element
   * @tparam FERequirements Type representing the requirements for finite element calculations
   * @tparam size Size of the stress vector
   * @tparam ScalarType Scalar type for numerical calculations
   */
  struct VonMises {
    template <typename ElementType, typename FERequirements, int size, typename ScalarType>
    double operator()(const ElementType& fe, const ResultRequirements<FERequirements>& req,
                      const Dune::FieldVector<ScalarType, size>& pos, [[maybe_unused]] int comp) const
        requires(size == 2) {
      ResultTypeMap<ScalarType> res_;
      fe.calculateAt(req, pos, res_);

      const auto& [resultType, sigma] = res_.getSingleResult();
      assert(resultType == ResultType::cauchyStress or resultType == ResultType::PK2Stress);
      const auto s_x  = sigma(0, 0);
      const auto s_y  = sigma(1, 0);
      const auto s_xy = sigma(2, 0);

      return std::sqrt(std::pow(s_x, 2) + Dune::power(s_y, 2) - s_x * s_y + 3 * Dune::power(s_xy, 2));
    }

    /**
     * @brief Get the name of the result type (VonMises)
     * @return String representing the name
     */
    static std::string name() { return "VonMises"; }

    /**
     * @brief Get the number of components in the result (always 1 for VonMises)
     * @return Number of components
     */
    static int ncomps() { return 1; }
  };

  /**
   * @brief Struct for calculating principal stresses
   * @ingroup resultevaluators
   *  @details The PrincipalStress struct provides a function call operator to calculate principal stresses.
   * @remark  Only 2D stresses are supported
   *
   * @tparam ElementType Type of the finite element
   * @tparam FERequirements Type representing the requirements for finite element calculations
   * @tparam size Size of the stress vector
   * @tparam ScalarType Scalar type for numerical calculations
   */
  struct PrincipalStress {
    template <typename ElementType, typename FERequirements, int size, typename ScalarType>
    double operator()(const ElementType& fe, const ResultRequirements<FERequirements>& req,
                      const Dune::FieldVector<ScalarType, size>& pos, int comp) const requires(size == 2) {
      ResultTypeMap<ScalarType> res_;
      fe.calculateAt(req, pos, res_);

      const auto& [resultType, sigma] = res_.getSingleResult();
      assert(resultType == ResultType::cauchyStress or resultType == ResultType::PK2Stress);
      const auto s_x  = sigma(0, 0);
      const auto s_y  = sigma(1, 0);
      const auto s_xy = sigma(2, 0);

      auto t1 = (s_x + s_y) / 2;
      auto t2 = std::sqrt(Dune::power((s_x - s_y) / 2, 2) + Dune::power(s_xy, 2));

      return comp == 0 ? t1 + t2 : t1 - t2;
    }

    /**
     * @brief Get the name of the result type (PrincipalStress)
     * @return String representing the name
     */
    static std::string name() { return "PrincipalStress"; }

    /**
     * @brief Get the number of components in the result (always 2 for PrincipalStress)
     * @return Number of components
     */
    static int ncomps() { return 2; }
  };

}  // namespace Ikarus::ResultEvaluators
