// SPDX-FileCopyrightText: 2021-2025 The Ikarus Developers ikarus@ibb.uni-stuttgart.de
// SPDX-License-Identifier: LGPL-3.0-or-later

/**
 * \file ogden.hh
 * \brief Implementation of the Ogden material model.
 * \ingroup materials
 */

#pragma once

#include <ikarus/finiteelements/mechanics/materials/materialhelpers.hh>
#include <ikarus/finiteelements/mechanics/materials/tags.hh>
#include <ikarus/utils/tensorutils.hh>

namespace Ikarus::Materials {

/**
 * \brief Implementation of the Ogden material model.
 * \ingroup materials
 *
 * \details If total principal stretches are used, the energy is computed via
 * \f[ \hat{\Psi}(\la_1, \la_2, \la_3) = \sum_{p=1}^n{
 \frac{\mu_n}{\al_n} ({\la_1}^{\al_n} + {\la_2}^{\al_n} +
 {\la_3}^{\al_n} - 3) - \mu_n \ln J}, \f]
 * with \f$ J = \la_1 \la_2 \la_3 \f$.
 *
 * If deviatoric principal stretches (\f$ \bar{\la_i} = \la_i J^{-1/3} \f$)
 * are used, the energy is computed via
 * \f[ \hat{\Psi}(\la_1, \la_2, \la_3) = \sum_{p=1}^n{
 \frac{\mu_n}{\al_n}({\bar{\la_1}}^{\al_n} + {\bar{\la_2}}^{\al_n} +
 {\bar{\la_3}}^{\al_n} - 3) }. \f]
 *
 * \remark See \cite bergstromMechanicsSolidPolymers2015 for details on this material.
 * \tparam ST_ The underlying scalar type.
 * \tparam n Number of ogden parameters
 * \tparam tag Type of principal stretch quantity, either total stretches or deviatoric stretches
 */
template <typename ST_, int n, PrincipalStretchTags tag>
struct OgdenT
{
  using ScalarType = ST_;

  template <typename ST = ScalarType>
  using PrincipalStretches = Eigen::Vector<ST, 3>;

  static constexpr PrincipalStretchTags stretchTag = tag;
  static constexpr int numMatParameters            = n;
  static constexpr int dim                         = 3;
  static constexpr bool usesDeviatoricStretches    = stretchTag == PrincipalStretchTags::deviatoric;

  template <typename ST = ScalarType>
  using FirstDerivative = Eigen::Vector<ST, dim>;
  template <typename ST = ScalarType>
  using SecondDerivative = Eigen::Matrix<ST, dim, dim>;

  using MaterialParameters = std::array<double, numMatParameters>;
  using MaterialExponents  = std::array<double, numMatParameters>;

  [[nodiscard]] constexpr static std::string name() noexcept {
    return "Ogden (n = " + std::to_string(numMatParameters) + ", stretch type = " + toString(tag) + ")";
  }

  /**
   * \brief Constructor for OgdenT
   *
   * \param mpt material parameters (array of mu values)
   * \param mex material exponents (array of alpha values)
   */
  explicit OgdenT(const MaterialParameters& mpt, const MaterialExponents& mex)
      : materialParameters_{mpt},
        materialExponents_{mex} {}

  /**
   * \brief Returns the material parameters (mu values) stored in the material
   */
  MaterialParameters materialParametersImpl() const { return materialParameters_; }

  /**
   * \brief Returns the material exponents (alpha values) stored in the material
   */
  const MaterialExponents& materialExponents() const { return materialExponents_; }

  /**
   * \brief Computes the stored energy in the Ogden material model.
   *
   * \param lambda principal stretches
   * \tparam ST the scalartype of the principal stretches
   * \return ScalarType
   */

  template <typename ST>
  auto storedEnergyImpl(const PrincipalStretches<ST>& lambda) const {
    auto& mu    = materialParameters_;
    auto& alpha = materialExponents_;

    ST energy{};

    if constexpr (usesDeviatoricStretches) {
      auto lambdaBar = Impl::deviatoricStretches(lambda);
      for (auto i : parameterRange())
        energy += mu[i] / alpha[i] * (lambdaBar.array().pow(alpha[i]).sum() - 3);
    } else {
      auto J    = lambda[0] * lambda[1] * lambda[2];
      auto logJ = log(J);
      for (auto i : parameterRange())
        energy += mu[i] / alpha[i] * (lambda.array().pow(alpha[i]).sum() - 3) - mu[i] * logJ;
    }
    return energy;
  }

  /**
   * \brief Computes the first derivative of the stored energy function w.r.t. the total principal stretches.
   *
   * \param lambda principal stretches
   * \tparam ST the scalartype of the principal stretches
   * \return ScalarType
   */

  template <typename ST>
  auto firstDerivativeImpl(const PrincipalStretches<ST>& lambda) const {
    auto& mu       = materialParameters_;
    auto& alpha    = materialExponents_;
    auto dWdLambda = FirstDerivative<ST>::Zero().eval();

    if constexpr (usesDeviatoricStretches) {
      auto lambdaBar    = Impl::deviatoricStretches(lambda);
      auto dWdLambdaBar = Eigen::Array<ST, dim, 1>::Zero().eval();

      for (const auto j : parameterRange())
        dWdLambdaBar += mu[j] * lambdaBar.array().pow(alpha[j] - 1);

      const ST sumLambdaBar = (lambdaBar.array() * dWdLambdaBar).sum();
      dWdLambda             = (lambdaBar.array() * dWdLambdaBar - (1.0 / 3.0) * sumLambdaBar) / lambda.array();
    } else
      for (const auto j : parameterRange())
        dWdLambda.array() += (mu[j] * (lambda.array().pow(alpha[j]) - 1)) / lambda.array();

    return dWdLambda;
  }

  /**
   * \brief Computes the second derivatives of the stored energy function w.r.t. the total principal stretches.
   *
   * \param lambda principal stretches
   * \tparam ST the scalartype of the principal stretches
   * \return ScalarType
   */

  template <typename ST>
  auto secondDerivativeImpl(const PrincipalStretches<ST>& lambda) const {
    auto& mu    = materialParameters_;
    auto& alpha = materialExponents_;
    auto dS     = SecondDerivative<ST>::Zero().eval();

    if constexpr (usesDeviatoricStretches) {
      const auto lambdaBar = Impl::deviatoricStretches(lambda);
      const auto dWdLambda = firstDerivativeImpl(lambda);

      auto lambdaBarPowSum = Eigen::Array<ST, dim, 1>::Zero().eval();
      for (const auto p : parameterRange())
        lambdaBarPowSum[p] = lambdaBar.array().pow(alpha[p]).sum();

      for (const auto a : dimensionRange()) {
        for (const auto b : dimensionRange()) {
          if (a == b)
            for (const auto p : parameterRange())
              dS(a, b) += mu[p] * alpha[p] * (1.0 / 3.0 * pow(lambdaBar[a], alpha[p]) + 1.0 / 9.0 * lambdaBarPowSum[p]);

          else
            for (const auto p : parameterRange())
              dS(a, b) += mu[p] * alpha[p] *
                          (-(1.0 / 3.0) * (pow(lambdaBar[a], alpha[p]) + pow(lambdaBar[b], alpha[p])) +
                           1.0 / 9.0 * lambdaBarPowSum[p]);

          dS(a, b) *= 1.0 / (lambda[a] * lambda[b]);
          if (a == b)
            dS(a, b) -= (2.0 / lambda[a]) * dWdLambda[a];
        }
      }
    } else {
      for (const auto j : parameterRange()) {
        dS.diagonal().array() +=
            (-2 * mu[j] * (lambda.array().pow(alpha[j]) - 1) + mu[j] * lambda.array().pow(alpha[j]) * alpha[j]) /
            lambda.array().square();
      }
    }
    return dS;
  }

  /**
   * \brief Rebinds the material to a different scalar type.
   * \tparam STO The target scalar type.
   * \return OgdenT<ScalarTypeOther> The rebound Ogden material.
   */
  template <typename STO>
  auto rebind() const {
    return OgdenT<STO, numMatParameters, stretchTag>(materialParameters_, materialExponents_);
  }

private:
  MaterialParameters materialParameters_;
  MaterialExponents materialExponents_;

  inline static constexpr auto parameterRange() { return Dune::range(numMatParameters); }
  inline static constexpr auto dimensionRange() { return Dune::range(dim); }
};

/**
 * \brief Alias for OgdenT with double as the default scalar type.
 */
template <int n, PrincipalStretchTags tag>
using Ogden = OgdenT<double, n, tag>;

} // namespace Ikarus::Materials
