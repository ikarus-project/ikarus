// SPDX-FileCopyrightText: 2021-2025 The Ikarus Developers ikarus@ibb.uni-stuttgart.de
// SPDX-License-Identifier: LGPL-3.0-or-later

/**
 * \file BlatzKo.hh
 * \brief Implementation of the BlatzKo material model.
 * \ingroup  materials
 */

#pragma once

#include <ikarus/finiteelements/mechanics/materials/interface.hh>
#include <ikarus/finiteelements/mechanics/materials/materialhelpers.hh>
#include <ikarus/finiteelements/mechanics/materials/tags.hh>
#include <ikarus/utils/tensorutils.hh>

namespace Ikarus::Materials {

/**
 * \brief Implementation of the Blatz-Ko material model.
 * \ingroup materials
 *
 * \details The energy is computed as
 * \f[ \hat{\Psi}(\la_1, \la_2, \la_3) = \frac{\mu}{2}(
 \frac{1}{\la_1^2} + \frac{1}{\la_2^2} + \frac{1}{\la_3^2}
 + 2\la_1 \la_2 \la_3 - 5). \f]
 *
 * \remark A special Blatz-Ko material model is implemented here. It assumes material parameters which corresponds to a
 * Poisson's ratio of 0.25.
 *
 * \tparam ST_ The underlying scalar type.
 */
template <typename ST_>
struct BlatzKoT
{
  using ScalarType = ST_;

  template <typename ST = ScalarType>
  using PrincipalStretches = Eigen::Vector<ST, 3>;

  static constexpr int dim         = 3;
  static constexpr auto stretchTag = PrincipalStretchTags::total;

  template <typename ST = ScalarType>
  using FirstDerivative = Eigen::Vector<ST, dim>;
  template <typename ST = ScalarType>
  using SecondDerivative = Eigen::Matrix<ST, dim, dim>;

  using MaterialParameters = double;

  [[nodiscard]] constexpr static std::string name() noexcept { return "BlatzKo"; }

  /**
   * \brief Constructor for BlatzKoT.
   * \param mu The shear modulus.
   */
  explicit BlatzKoT(MaterialParameters mu)
      : mu_{mu} {}

  /**
   * \brief Returns the material parameters stored in the material
   */
  MaterialParameters materialParametersImpl() const { return mu_; }

  /**
   * \brief Computes the stored energy in the BlatzKo material model.
   *
   * \param lambda principal stretches
   * \tparam ST the scalartype of the principal stretches
   * \return ScalarType
   */
  template <typename ST>
  ST storedEnergyImpl(const PrincipalStretches<ST>& lambda) const {
    return mu_ / 2 * (lambda.cwiseInverse().squaredNorm() + 2 * lambda.prod() - 5);
  }

  /**
   * \brief Computes the first derivative of the stored energy function w.r.t. the total principal stretches.
   *
   * \param lambda principal stretches
   * \return ScalarType
   */
  template <typename ST>
  FirstDerivative<ST> firstDerivativeImpl(const PrincipalStretches<ST>& lambda) const {
    auto dWdLambda = FirstDerivative<ST>::Zero().eval();
    const ST J     = Impl::determinantFromPrincipalValues(lambda);

    return mu_ * (-lambda.cwisePow(3).cwiseInverse() + (J * lambda.cwiseInverse()));
  }

  /**
   * \brief Computes the second derivatives of the stored energy function w.r.t. the total principal stretches.
   *
   * \param lambda principal stretches
   * \tparam ST the scalartype of the principal stretches
   * \return ScalarType
   */
  template <typename ST>
  SecondDerivative<ST> secondDerivativeImpl(const PrincipalStretches<ST>& lambda) const {
    auto dS    = SecondDerivative<ST>::Zero().eval();
    const ST J = Impl::determinantFromPrincipalValues(lambda);

    auto lam      = lambda.array();
    dS            = J / (lambda * lambda.transpose()).array();
    dS.diagonal() = lam.square().inverse() * (lam.square().inverse() - J) + (3.0 / lam.pow(4));
    dS.array() *= mu_;

    return dS;
  }

  /**
   * \brief Rebinds the material to a different scalar type.
   * \tparam STO The target scalar type.
   * \return BlatzKoT<ScalarTypeOther> The rebound BlatzKo material.
   */
  template <typename STO>
  auto rebind() const {
    return BlatzKoT<STO>(mu_);
  }

private:
  MaterialParameters mu_;

  inline static constexpr auto dimensionRange() { return Dune::range(dim); }
};

/**
 * \brief Alias for BlatzKoT with double as the default scalar type.
 */
using BlatzKo = BlatzKoT<double>;

} // namespace Ikarus::Materials
