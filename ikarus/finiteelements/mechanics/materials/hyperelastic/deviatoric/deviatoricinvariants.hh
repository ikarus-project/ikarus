// SPDX-FileCopyrightText: 2021-2025 The Ikarus Developers ikarus@ibb.uni-stuttgart.de
// SPDX-License-Identifier: LGPL-3.0-or-later

/**
 * \file deviatoricinvariants.hh
 * \brief Implementation of the computation of the deviatoric invariants and its derivatives.
 * \ingroup materials
 */

#pragma once

#include <ikarus/finiteelements/mechanics/materials/materialhelpers.hh>
#include <ikarus/finiteelements/mechanics/materials/tags.hh>
#include <ikarus/utils/tensorutils.hh>

namespace Ikarus::Materials {
/**
 * \brief Implementation of the deviatoric invariants and its derivatives.
 * \ingroup materials
 *
 * \details The three invariants are
 * \f[ I_1 = \la_1^2 + \la_2^2 + \la_3^2, \f]
 * \f[ I_2 = \la_1^2 \la_2^2 + \la_2^2 \la_3^2 +
 \la_1^2 \la_3^2, \f]
 * \f[ I_3 = \la_1^2 \la_2^2 \la_3^2. \f]
 *
 * The deviatoric invariants are then defined as
 * \f[ W_1 = I_1 I_3^{-1/3}, \f]
 * \f[ W_2 = I_2 I_3^{-2/3}, \f]
 *
 * \tparam PS Type of the principal stretches.
 */
template <typename PS>
struct DeviatoricInvariants
{
  using PrincipalStretches = PS;
  using ScalarType         = PrincipalStretches::value_type;
  static constexpr int dim = PrincipalStretches::RowsAtCompileTime;
  using Invariants         = PrincipalStretches;

  using FirstDerivative  = Eigen::Vector<ScalarType, dim>;
  using SecondDerivative = Eigen::Matrix<ScalarType, dim, dim>;

  /**
   * \brief Constructor for DeviatoricInvariants
   *
   * \param lambda The principal stretches.
   */
  explicit DeviatoricInvariants(const PrincipalStretches& lambda)
      : lambda_{lambda} {}

  /**
   * \brief Computation of value of the deviatoric invariants
   */
  auto value() const {
    const Invariants& invariants = Impl::invariants(lambda_);
    ScalarType W1                = invariants[0] * pow(invariants[2], -1.0 / 3.0);
    ScalarType W2                = invariants[1] * pow(invariants[2], -2.0 / 3.0);
    return std::array<ScalarType, 2>{W1, W2};
  }

  /**
   * \brief Computation of the first derivatives of the deviatoric invariants w.r.t the total principal stretches.
   */
  auto firstDerivative() const {
    const Invariants& invariants            = Impl::invariants(lambda_);
    auto [I1, I2, I3, I3Pow1by3, I3Pow2by3] = computeInvariants(invariants);

    FirstDerivative dW1dLambda = 2.0 * (3.0 * lambda_.cwisePow(2).array() - I1) / (3.0 * lambda_.array() * I3Pow1by3);
    FirstDerivative dW2dLambda =
        -2.0 * (3.0 * I3 * lambda_.cwisePow(2).cwiseInverse().array() - I2) / (3.0 * lambda_.array() * I3Pow2by3);

    return std::make_pair(dW1dLambda, dW2dLambda);
  }

  /**
   * \brief Computation of the second derivatives of the deviatoric invariants w.r.t the total principal stretches.
   */
  auto secondDerivative() const {
    const Invariants& invariants            = Impl::invariants(lambda_);
    auto ddW1dLambda                        = SecondDerivative::Zero().eval();
    auto ddW2dLambda                        = SecondDerivative::Zero().eval();
    auto [I1, I2, I3, I3Pow1by3, I3Pow2by3] = computeInvariants(invariants);

    auto lam               = lambda_.array();
    ddW1dLambda.diagonal() = (2.0 / 9.0) * (5 * I1 - 3 * lam.square()) / (lam.square() * I3Pow1by3);
    ddW2dLambda.diagonal() = (2.0 / 9.0) * ((15.0 * I3 / lam.square()) - I2) / (lam.square() * I3Pow2by3);
    for (const auto j : dimensionRange()) {
      for (const auto i : dimensionRange())
        if (i != j) {
          ddW1dLambda(i, j) = (4.0 / 9.0) * (I1 - 3 * (pow(lambda_[i], 2.0) + pow(lambda_[j], 2.0))) /
                              (lambda_[i] * lambda_[j] * I3Pow1by3);
          ddW2dLambda(i, j) = (-4.0 / 9.0) * (2.0 * I2 - 3 * pow(lambda_[i], 2.0) * pow(lambda_[j], 2.0)) /
                              (lambda_[i] * lambda_[j] * I3Pow2by3);
        }
    }
    return std::make_pair(ddW1dLambda, ddW2dLambda);
  }

private:
  PrincipalStretches lambda_;

  inline static constexpr auto dimensionRange() { return Dune::range(dim); }

  auto computeInvariants(const Invariants& invariants) const {
    ScalarType I1        = invariants[0];
    ScalarType I2        = invariants[1];
    ScalarType I3        = invariants[2];
    ScalarType I3Pow1by3 = pow(I3, 1.0 / 3.0);
    ScalarType I3Pow2by3 = pow(I3, 2.0 / 3.0);

    return std::make_tuple(I1, I2, I3, I3Pow1by3, I3Pow2by3);
  }
};

} // namespace Ikarus::Materials
