// SPDX-FileCopyrightText: 2021-2025 The Ikarus Developers ikarus@ibb.uni-stuttgart.de
// SPDX-License-Identifier: LGPL-3.0-or-later

/**
 * \file invariantbased.hh
 * \brief Implementation of the InvariantBased material model.
 * \ingroup materials
 */

#pragma once

#include <ikarus/finiteelements/mechanics/materials/hyperelastic/deviatoric/deviatoricinvariants.hh>
#include <ikarus/finiteelements/mechanics/materials/materialhelpers.hh>
#include <ikarus/finiteelements/mechanics/materials/tags.hh>
#include <ikarus/utils/tensorutils.hh>

namespace Ikarus::Materials {

/**
 * \brief Implementation of the InvariantBased material model.
 * \ingroup materials
 *
 * \details This is a general model based on deviatoric invariants. It can be used to derive other specific material
 * models, for instance Mooney-Rivlin and Yeoh models.
 * The energy is computed as
 * \f[ \hat{\Psi}(\la_1, \la_2, \la_3) = \sum_{p+q=1}^n{C_{pq} (W_1 - 3)^p (W_2 - 3)^q}, \quad p,q \geq 0. \f]
 *
 * \remark See \cite bergstromMechanicsSolidPolymers2015 for details on this material. For information on the deviatoric
 * invariants \f$ W_1, W_2 \f$, see \ref DeviatoricInvariants.
 *
 * \tparam ST_ The underlying scalar type.
 * \tparam n Number of material parameters
 */
template <typename ST_, int n>
struct InvariantBasedT
{
  using ScalarType = ST_;

  template <typename ST = ScalarType>
  using PrincipalStretches = Eigen::Vector<ST, 3>;
  template <typename ST = ScalarType>
  using Invariants = PrincipalStretches<ST>;

  static constexpr int dim              = 3;
  static constexpr int numMatParameters = n;

  using Exponents          = std::array<std::size_t, numMatParameters>;
  using MaterialParameters = std::array<double, numMatParameters>;

  template <typename ST = ScalarType>
  using FirstDerivative = Eigen::Vector<ST, dim>;
  template <typename ST = ScalarType>
  using SecondDerivative = Eigen::Matrix<ST, dim, dim>;

  static constexpr auto stretchTag = PrincipalStretchTags::deviatoric;

  [[nodiscard]] constexpr static std::string name() noexcept {
    return "InvariantBased (n = " + std::to_string(numMatParameters);
  }

  MaterialParameters materialParametersImpl() const { return matParameters_; }

  const Exponents& pExponents() const { return pex_; }

  const Exponents& qExponents() const { return qex_; }

  /**
   * \brief Constructor for InvariantBasedT
   *
   * \param pex Array of exponents related to the first invariant.
   * \param qex Array of exponents related to the second invariant.
   * \param matParameters Array of material parameters.
   */
  explicit InvariantBasedT(const Exponents& pex, const Exponents& qex, const MaterialParameters& matParameters)
      : pex_{pex},
        qex_{qex},
        matParameters_{matParameters} {
    checkExponents();
  }

  /**
   * \brief Computes the stored energy in the InvariantBased material model.
   *
   * \param lambda principal stretches
   * \tparam ST the scalartype of the principal stretches
   * \return ScalarType the energy
   */
  template <typename ST = ScalarType>
  ST storedEnergyImpl(const PrincipalStretches<ST>& lambda) const {
    const Invariants<ST>& invariants = Impl::invariants(lambda);
    ST energy{0.0};
    const auto& devInvariants = DeviatoricInvariants<PrincipalStretches<ST>>(lambda);
    auto [W1, W2]             = devInvariants.value();
    W1 -= 3.0;
    W2 -= 3.0;

    for (const auto i : parameterRange())
      energy += matParameters_[i] * pow(W1, pex_[i]) * pow(W2, qex_[i]);

    return energy;
  }

  /**
   * \brief Computes the first derivative of the stored energy function w.r.t. the total principal stretches.
   *
   * \param lambda principal stretches
   * \tparam ST the scalartype of the principal stretches
   * \return FirstDerivative the first derivative
   */
  template <typename ST = ScalarType>
  FirstDerivative<ST> firstDerivativeImpl(const PrincipalStretches<ST>& lambda) const {
    const Invariants<ST>& invariants = Impl::invariants(lambda);
    auto dWdLambda                   = FirstDerivative<ST>::Zero().eval();

    const auto& mu            = matParameters_;
    const auto& devInvariants = DeviatoricInvariants<PrincipalStretches<ST>>(lambda);
    auto [W1, W2]             = devInvariants.value();
    W1 -= 3.0;
    W2 -= 3.0;
    const auto& [dW1dLambda, dW2dLambda] = devInvariants.firstDerivative();

    for (const auto p : parameterRange())
      for (const auto k : dimensionRange()) {
        auto W1pm1p = powerAndMultiply<ST>(W1, pex_[p] - 1.0, pex_[p]);
        auto W2qm1q = powerAndMultiply<ST>(W2, qex_[p] - 1.0, qex_[p]);
        dWdLambda[k] +=
            mu[p] * ((W1pm1p * pow(W2, qex_[p]) * dW1dLambda[k]) + (pow(W1, pex_[p]) * W2qm1q * dW2dLambda[k]));
      }

    return dWdLambda;
  }

  /**
   * \brief Computes the second derivatives of the stored energy function w.r.t. the total principal stretches.
   *
   * \param lambda principal stretches.
   * \tparam ST the scalartype of the principal stretches
   * \return SecondDerivative the second derivative
   */
  template <typename ST = ScalarType>
  SecondDerivative<ST> secondDerivativeImpl(const PrincipalStretches<ST>& lambda) const {
    const Invariants<ST>& invariants = Impl::invariants(lambda);
    auto dS                          = SecondDerivative<ST>::Zero().eval();

    const auto& mu            = matParameters_;
    const auto& devInvariants = DeviatoricInvariants<PrincipalStretches<ST>>(lambda);
    auto [W1, W2]             = devInvariants.value();
    W1 -= 3.0;
    W2 -= 3.0;
    const auto& [dW1dLambda, dW2dLambda]   = devInvariants.firstDerivative();
    const auto& [ddW1dLambda, ddW2dLambda] = devInvariants.secondDerivative();

    for (const auto p : parameterRange())
      for (const auto i : dimensionRange())
        for (const auto j : dimensionRange()) {
          auto W1pm1p       = powerAndMultiply<ST>(W1, pex_[p] - 1.0, pex_[p]);
          auto W2qm1q       = powerAndMultiply<ST>(W2, qex_[p] - 1.0, qex_[p]);
          auto W1pm2pp      = powerAndMultiply<ST>(W1, pex_[p] - 2.0, pex_[p] * (pex_[p] - 1.0));
          auto W2qm2qq      = powerAndMultiply<ST>(W2, qex_[p] - 2.0, qex_[p] * (qex_[p] - 1.0));
          auto dW1W2dlambda = dW1dLambda[i] * dW2dLambda[j] + dW1dLambda[j] * dW2dLambda[i];
          auto factor1      = (W2qm1q * dW1W2dlambda + pow(W2, qex_[p]) * ddW1dLambda(i, j)) * W1pm1p * mu[p];
          auto factor2      = W2qm1q * ddW2dLambda(i, j) * pow(W1, pex_[p]) * mu[p];
          auto factor3      = W1pm2pp * pow(W2, qex_[p]) * dW1dLambda[i] * dW1dLambda[j] * mu[p];
          auto factor4      = W2qm2qq * pow(W1, pex_[p]) * dW2dLambda[i] * dW2dLambda[j] * mu[p];
          dS(i, j) += factor1 + factor2 + factor3 + factor4;
          if (i == j) {
            auto factor5 =
                mu[p] * (W1pm1p * pow(W2, qex_[p]) * dW1dLambda[i] + pow(W1, pex_[p]) * W2qm1q * dW2dLambda[i]);
            dS(i, j) -= (1.0 / lambda[i]) * factor5;
          }
        }
    return dS;
  }

  /**
   * \brief Rebinds the material to a different scalar type.
   * \tparam STO The target scalar type.
   * \return InvariantBasedT<ScalarTypeOther> The rebound InvariantBased material.
   */
  template <typename STO>
  auto rebind() const {
    return InvariantBasedT<STO, numMatParameters>(pex_, qex_, matParameters_);
  }

private:
  Exponents pex_, qex_;
  MaterialParameters matParameters_;

  inline static constexpr auto parameterRange() { return Dune::range(numMatParameters); }
  inline static constexpr auto dimensionRange() { return Dune::range(dim); }

  /**
   * \brief A function that computes \f$ x^p * m \f$.
   *
   * \details In the undeformed configuration, all principal stretches are equal to one. In this scenario, the
   * deviatoric invariants W1 and W2 are zero. For any positive integer n, 0^{-n} is inf. Multiplying this with 0.0
   * again leads to nan. In order to circumvent this, in such a scenario, zero is returned instead of nan.
   *
   * \remark If inf is multiplied with any other number apart from zero, nan is returned.
   *
   * \param x Base number.
   * \param p Exponent to which the base number is raised.
   * \param m Scalar number to be multiplied to \f$ x^p \f$.
   * \return ScalarType Either 0.0 or \f$ x^p * m \f$.
   */
  template <typename ST>
  ST powerAndMultiply(ST x, int p, int m) const {
    if (m == 0)
      return 0; // anything multiplied with 0 is 0
    if (p == 0)
      return m; // x^0 * m = 1 * m = m
    const double epsilon = std::numeric_limits<double>::epsilon();
    assert(not(std::abs(x) < epsilon and p < 0) && "Raising zero to a negative power results in NaN.");
    return std::pow(x, p) * m;
  }

  /** \brief check that: \f$ \lnot (q_i = 0 \land p_i = 0) \forall i \f$ */
  void checkExponents() const {
    for (const auto i : parameterRange())
      if (pex_[i] == 0 and qex_[i] == 0)
        DUNE_THROW(Dune::InvalidStateException, "The exponent q" + std::to_string(i) + "and the exponent p" +
                                                    std::to_string(i) + "should not be zero at the same time.");
  }
};

/**
 * \brief Alias for InvariantBasedT with double as the default scalar type.
 */
template <int n>
using InvariantBased = InvariantBasedT<double, n>;

} // namespace Ikarus::Materials
