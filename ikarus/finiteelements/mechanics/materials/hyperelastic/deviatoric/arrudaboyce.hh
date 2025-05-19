// SPDX-FileCopyrightText: 2021-2025 The Ikarus Developers ikarus@ibb.uni-stuttgart.de
// SPDX-License-Identifier: LGPL-3.0-or-later

/**
 * \file arrudaboyce.hh
 * \brief Implementation of the ArrudaBoyce material model.
 * \ingroup materials
 */

#pragma once

#include <ikarus/finiteelements/mechanics/materials/hyperelastic/deviatoric/deviatoricinvariants.hh>
#include <ikarus/finiteelements/mechanics/materials/materialhelpers.hh>
#include <ikarus/finiteelements/mechanics/materials/tags.hh>
#include <ikarus/utils/tensorutils.hh>

namespace Ikarus {

///< Structure representing material parameters for the Arrudy-Boyce material model.
struct ArrudaBoyceMatParameters
{
  double mu;      ///< Denotes the shear modulus.
  double lambdaM; ///< Denotes the maximum (fully extended) stretch that a molecule is exposed to.
};
} // namespace Ikarus

namespace Ikarus::Materials {

/**
 * \brief Implementation of the ArrudaBoyce material model (also referred as Eight-Chain model).
 * \ingroup materials
 *
 * \details The energy is computed as
 * \f[ \hat{\Psi}(\la_1, \la_2, \la_3) = \mu  \sum_{p=0}^4{\al_p  \beta^p}  (W_1^{p+1} - 3^{p+1}), \f]
 * with \f$ \beta = \frac{1}{\la_m^2} , \al_0 = \frac{1}{2}, \al_1 = \frac{1}{20}, \al_2 = \frac{11}{1050},
 \al_3 = \frac{19}{7000}, \al_4 = \frac{519}{673750} \f$.
 *
 * The first derivatives w.r.t the total principal stretches are
 * \f[ \fracpt{\Psi}{\la_i} = \mu \sum_{p=0}^4{\al_p  \beta^p  W_1^p (p+1)
 \fracpt{W_1}{\la_i}}. \f]
 *
 * The second derivatives w.r.t the total principal stretches are
 * \f[ \fracpt{^2 \Psi}{\la_i\partial\la_j} = \mu \sum_{p=0}^4{
  \al_p \beta^p (W_1^p (p+1) \fracpt{^2 W_1}{\la_i\partial\la_j} +
 W_1^{p-1} p (p+1) \fracpt{W_1}{\la_i} \fracpt{W_1}{\la_j} - \delta_{ij}
 \frac{1}{\la_i} W_1^p (p+1) \fracpt{W_1}{\la_i}) } \f]
 *
 * \remark See \cite hiermaierStructuresCrashImpact2010 and \cite bergstromMechanicsSolidPolymers2015 for details on
 * this material. For information on the deviatoric invariant \f$ W_1 \f$, see \ref DeviatoricInvariants.
 *
 *\tparam ST_ The underlying scalar type.
 */
template <typename ST_>
struct ArrudaBoyceT
{
  static constexpr int dim = 3;
  using ScalarType         = ST_;

  template <typename ST = ScalarType>
  using PrincipalStretches = Eigen::Vector<ST, 3>;
  template <typename ST = ScalarType>
  using Invariants = PrincipalStretches<ST>;

  template <typename ST = ScalarType>
  using FirstDerivative = Eigen::Vector<ST, dim>;
  template <typename ST = ScalarType>
  using SecondDerivative = Eigen::Matrix<ST, dim, dim>;

  using MaterialParameters      = ArrudaBoyceMatParameters;
  static constexpr int numTerms = 5;

  static constexpr auto stretchTag = PrincipalStretchTags::deviatoric;

  [[nodiscard]] constexpr static std::string name() noexcept { return "ArrudaBoyce"; }

  /**
   * \brief Constructor for ArrudaBoyceT
   *
   * \param matPar The material parameters for the Arruda-Boyce material model.
   */
  explicit ArrudaBoyceT(const MaterialParameters& matPar)
      : matPar_{matPar} {}

  /**
   * \brief Returns the material parameters stored in the material
   */
  MaterialParameters materialParametersImpl() const { return matPar_; }

  /**
   * \brief Computes the stored energy in the ArrudaBoyce material model.
   * \details Using only the first five terms of the inverse Langevin function.
   *
   * \param lambda principal stretches
   * \tparam ST the scalartype of the principal stretches
   * \return ScalarType
   */
  template <typename ST>
  ST storedEnergyImpl(const PrincipalStretches<ST>& lambda) const {
    const Invariants<ST>& invariants = Impl::invariants(lambda);
    ST energy{0.0};
    const auto& devInvariants = DeviatoricInvariants<PrincipalStretches<ST>>(lambda);
    const auto W1             = devInvariants.value()[0];
    const auto mu             = matPar_.mu;
    const auto lambdaM        = matPar_.lambdaM;
    const auto beta           = 1 / pow(lambdaM, 2.0);

    for (const auto i : parameterRange())
      energy += alphas_[i] * pow(beta, i) * (pow(W1, i + 1) - pow(3, i + 1));
    energy *= mu;

    return energy;
  }

  /**
   * \brief Computes the first derivative of the stored energy function w.r.t. the total principal stretches.
   *
   * \param lambda principal stretches
   * \tparam ST the scalartype of the principal stretches
   * \return FirstDerivative The first derivatives of the stored energy function w.r.t. the total principal stretches.
   */
  template <typename ST>
  FirstDerivative<ST> firstDerivativeImpl(const PrincipalStretches<ST>& lambda) const {
    const Invariants<ST>& invariants = Impl::invariants(lambda);
    auto dWdLambda                   = FirstDerivative<ST>::Zero().eval();

    const auto& devInvariants = DeviatoricInvariants<PrincipalStretches<ST>>(lambda);
    const auto W1             = devInvariants.value()[0];
    const auto& dW1dLambda    = devInvariants.firstDerivative().first;
    const auto mu             = matPar_.mu;
    const auto lambdaM        = matPar_.lambdaM;
    const auto beta           = 1 / pow(lambdaM, 2.0);

    for (const auto j : parameterRange())
      for (const auto k : dimensionRange())
        dWdLambda[k] += mu * alphas_[j] * pow(beta, j) * pow(W1, j) * dW1dLambda[k] * (j + 1);

    return dWdLambda;
  }

  /**
   * \brief Computes the second derivatives of the stored energy function w.r.t. the total principal stretches.
   *
   * \param lambda principal stretches
   * \tparam ST the scalartype of the principal stretches
   * \return SecondDerivative The second derivatives of the stored energy function w.r.t. the total principal stretches.
   */
  template <typename ST>
  SecondDerivative<ST> secondDerivativeImpl(const PrincipalStretches<ST>& lambda) const {
    const Invariants<ST>& invariants = Impl::invariants(lambda);
    auto dS                          = SecondDerivative<ST>::Zero().eval();

    const auto& devInvariants = DeviatoricInvariants<PrincipalStretches<ST>>(lambda);
    const auto W1             = devInvariants.value()[0];
    const auto& dW1dLambda    = devInvariants.firstDerivative().first;
    const auto& ddW1dLambda   = devInvariants.secondDerivative().first;
    const auto mu             = matPar_.mu;
    const auto lambdaM        = matPar_.lambdaM;
    const auto beta           = 1 / pow(lambdaM, 2.0);

    const auto dW1dLambdaDyad = dyadic(dW1dLambda, dW1dLambda);
    for (const auto p : parameterRange()) {
      const auto factor1 = mu * alphas_[p] * pow(beta, p);
      const auto factor2 = pow(W1, p) * ddW1dLambda * (p + 1);
      const auto factor3 = pow(W1, p - 1) * dW1dLambdaDyad * p * (p + 1);
      dS += factor1 * (factor2 + factor3);
      dS.diagonal() -= (dW1dLambda.array() / lambda.array() * factor1 * pow(W1, p) * (p + 1)).matrix();
    }

    return dS;
  }

  /**
   * \brief Rebinds the material to a different scalar type.
   * \tparam STO The target scalar type.
   * \return ArrudaBoyceT<ScalarTypeOther> The rebound ArrudaBoyce material.
   */
  template <typename STO>
  auto rebind() const {
    return ArrudaBoyceT<STO>(matPar_);
  }

private:
  MaterialParameters matPar_;
  static constexpr std::array<double, numTerms> alphas_ = {0.5, 1.0 / 20.0, 11.0 / 1050.0, 19.0 / 7000.0,
                                                           519.0 / 673750.0};

  inline static constexpr auto parameterRange() { return Dune::range(numTerms); }
  inline static constexpr auto dimensionRange() { return Dune::range(dim); }
};

/**
 * \brief Alias for ArrudaBoyceT with double as the default scalar type.
 */
using ArrudaBoyce = ArrudaBoyceT<double>;

} // namespace Ikarus::Materials
