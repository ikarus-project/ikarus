// SPDX-FileCopyrightText: 2021-2025 The Ikarus Developers ikarus@ibb.uni-stuttgart.de
// SPDX-License-Identifier: LGPL-3.0-or-later

/**
 * \file gent.hh
 * \brief Implementation of the Gent material model.
 * \ingroup materials
 */

#pragma once

#include <ikarus/finiteelements/mechanics/materials/hyperelastic/deviatoric/deviatoricinvariants.hh>
#include <ikarus/finiteelements/mechanics/materials/materialhelpers.hh>
#include <ikarus/finiteelements/mechanics/materials/tags.hh>
#include <ikarus/utils/tensorutils.hh>

namespace Ikarus {

///< Structure representing material parameters for the Gent material model.
struct GentMatParameters
{
  double mu; ///< Denotes the shear modulus.
  double Jm; ///< Denotes a dimensionless parameter that controls the extensibility of chains
};
} // namespace Ikarus

namespace Ikarus::Materials {

/**
 * \brief Implementation of the Gent material model.
 * \ingroup materials
 *
 * \details The energy is computed as
 * \f[ \hat{\Psi}(\la_1, \la_2, \la_3) = -\frac{\mu}{2} J_m \ln{(1-\frac{W_1 - 3}{J_m})}, \f]
 * with \f$ J_m > W_1-3\f$.
 *
 * The first derivatives w.r.t the total principal stretches are
 * \f[ \fracpt{\Psi}{\la_i} = \frac{\mu J_m}{2(J_m - W_1) + 6}\fracpt{W_1}{\la_i}. \f]
 *
 * The second derivatives w.r.t the total principal stretches are
 * \f[ \fracpt{^2 \Psi}{\la_i\partial\la_j} =
 \frac{\mu}{2\al} (\fracpt{^2 W_1}{\la_i\partial\la_j} +
 \frac{1}{\al J_m} \fracpt{W_1}{\la_i} \fracpt{W_1}{\la_j}) -
 \delta_{ij} \frac{\mu}{2\al\la_i} \fracpt{W_1}{\la_i}, \f]
 * with \f$ \al = 1 - \frac{W_1 - 3}{J_m} \f$.
 *
 * \remark See \cite bergstromMechanicsSolidPolymers2015 for details on this material. For information on the deviatoric
 * invariant \f$ W_1 \f$, see \ref DeviatoricInvariants.
 *
 * \tparam ST_ The underlying scalar type.
 */
template <typename ST_>
struct GentT
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

  using MaterialParameters = GentMatParameters;

  static constexpr auto stretchTag = PrincipalStretchTags::deviatoric;

  [[nodiscard]] constexpr static std::string name() noexcept { return "Gent"; }

  /**
   * \brief Constructor for GentT
   *
   * \param C material constant
   * \param lambdaM maximum stretch at which the polymer chain locks
   */
  explicit GentT(const MaterialParameters& matPar)
      : matPar_{matPar} {}

  /**
   * \brief Returns the material parameters stored in the material
   */
  MaterialParameters materialParametersImpl() const { return matPar_; }

  /**
   * \brief Computes the stored energy in the Gent material model.
   * \details Using only the first five terms of the inverse Langevin function.
   *
   * \param lambda principal stretches
   * \tparam ST the scalartype of the principal stretches
   * \return ScalarType
   */
  template <typename ST = ScalarType>
  ST storedEnergyImpl(const PrincipalStretches<ST>& lambda) const {
    const Invariants<ST>& invariants = Impl::invariants(lambda);
    const auto& devInvariants        = DeviatoricInvariants<PrincipalStretches<ST>>(lambda);
    const auto W1                    = devInvariants.value()[0];
    checkJm(W1);
    return -(matPar_.mu / 2.0) * matPar_.Jm * log(1.0 - ((W1 - 3.0) / matPar_.Jm));
  }

  /**
   * \brief Computes the first derivative of the stored energy function w.r.t. the total principal stretches.
   *
   * \param lambda principal stretches
   * \tparam ST the scalartype of the principal stretches
   * \return ScalarType
   */
  template <typename ST = ScalarType>
  FirstDerivative<ST> firstDerivativeImpl(const PrincipalStretches<ST>& lambda) const {
    const Invariants<ST>& invariants = Impl::invariants(lambda);
    const auto& devInvariants        = DeviatoricInvariants<PrincipalStretches<ST>>(lambda);
    const auto W1                    = devInvariants.value()[0];
    const auto& dW1dLambda           = devInvariants.firstDerivative().first;
    const auto mu                    = matPar_.mu;
    const auto Jm                    = matPar_.Jm;
    checkJm(W1);

    FirstDerivative<ST> dWdLambda = (mu * dW1dLambda * Jm) / (2.0 * (Jm - W1) + 6.0);

    return dWdLambda;
  }

  /**
   * \brief Computes the second derivatives of the stored energy function w.r.t. the total principal stretches.
   *
   * \param lambda principal stretches
   * \tparam ST the scalartype of the principal stretches
   * \return ScalarType
   */
  template <typename ST = ScalarType>
  SecondDerivative<ST> secondDerivativeImpl(const PrincipalStretches<ST>& lambda) const {
    const Invariants<ST>& invariants = Impl::invariants(lambda);
    auto dS                          = SecondDerivative<ST>::Zero().eval();

    const auto& devInvariants = DeviatoricInvariants<PrincipalStretches<ST>>(lambda);
    const auto W1             = devInvariants.value()[0];
    const auto& dW1dLambda    = devInvariants.firstDerivative().first;
    const auto& ddW1dLambda   = devInvariants.secondDerivative().first;
    const auto mu             = matPar_.mu;
    const auto Jm             = matPar_.Jm;
    checkJm(W1);

    for (const auto i : dimensionRange())
      for (const auto j : dimensionRange()) {
        auto factor1 = 1.0 - ((W1 - 3.0) / Jm);
        dS(i, j) += (mu / (2.0 * factor1)) * (ddW1dLambda(i, j) + (dW1dLambda[i] * dW1dLambda[j] / (factor1 * Jm)));
        if (i == j)
          dS(i, j) -= (mu / (2.0 * lambda[i] * factor1)) * dW1dLambda[i];
      }

    return dS;
  }

  /**
   * \brief Rebinds the material to a different scalar type.
   * \tparam STO The target scalar type.
   * \return GentT<ScalarTypeOther> The rebound Gent material.
   */
  template <typename STO>
  auto rebind() const {
    return GentT<STO>(matPar_);
  }

private:
  MaterialParameters matPar_;

  inline static constexpr auto dimensionRange() { return Dune::range(dim); }

  void checkJm(ScalarType W1) const {
    if (Dune::FloatCmp::le(matPar_.Jm, static_cast<double>(W1) - 3.0, 1e-14))
      DUNE_THROW(Dune::InvalidStateException, "The material parameter Jm should be greater than (W1 - 3)");
  }
};

/**
 * \brief Alias for GentT with double as the default scalar type.
 */
using Gent = GentT<double>;

} // namespace Ikarus::Materials
