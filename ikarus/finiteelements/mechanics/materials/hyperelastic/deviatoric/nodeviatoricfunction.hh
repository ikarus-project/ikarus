// SPDX-FileCopyrightText: 2021-2025 The Ikarus Developers ikarus@ibb.uni-stuttgart.de
// SPDX-License-Identifier: LGPL-3.0-or-later

/**
 * \file nodeviatoricfunction.hh
 * \brief Implementation of the NoDev material model.
 * \ingroup  materials
 */

#pragma once

#include <ikarus/finiteelements/mechanics/materials/interface.hh>
#include <ikarus/finiteelements/mechanics/materials/materialhelpers.hh>
#include <ikarus/finiteelements/mechanics/materials/tags.hh>
#include <ikarus/utils/tensorutils.hh>

namespace Ikarus::Materials {

/**
 * \brief Dummy class for no deviatoric function.
 * \ingroup materials
 *
 */
template <typename ST_>
struct NoDevT
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

  [[nodiscard]] constexpr static std::string name() noexcept { return "None"; }

  /**
   * \brief Returns the material parameters stored in the material (returns zero)
   */
  MaterialParameters materialParametersImpl() const { return 0.0; }

  /**
   * \brief Computes the stored energy (returns zero).
   *
   * \param lambda principal stretches
   * \tparam ST the scalartype of the principal stretches
   * \return ScalarType
   */
  template <typename ST>
  ST storedEnergyImpl(const PrincipalStretches<ST>& lambda) const {
    return ST{};
  }

  /**
   * \brief Computes the first derivative of the stored energy function w.r.t. the total principal stretches (returns
   * zero vector).
   *
   * \param lambda principal stretches
   * \return ScalarType
   */
  template <typename ST>
  FirstDerivative<ST> firstDerivativeImpl(const PrincipalStretches<ST>& lambda) const {
    return FirstDerivative<ST>::Zero().eval();
  }

  /**
   * \brief Computes the second derivatives of the stored energy function w.r.t. the total principal stretches (returns
   * zero matrix).
   *
   * \param lambda principal stretches
   * \tparam ST the scalartype of the principal stretches
   * \return ScalarType
   */
  template <typename ST>
  SecondDerivative<ST> secondDerivativeImpl(const PrincipalStretches<ST>& lambda) const {
    return SecondDerivative<ST>::Zero().eval();
  }

  /**
   * \brief Rebinds the material to a different scalar type.
   * \tparam STO The target scalar type.
   * \return NoDevT<ScalarTypeOther> The rebound NoDev material.
   */
  template <typename STO>
  auto rebind() const {
    return NoDevT<STO>();
  }

private:
  // inline static constexpr auto dimensionRange() { return Dune::range(dim); }
};

/**
 * \brief Alias for NoDevT with double as the default scalar type.
 */
using NoDev = NoDevT<double>;

} // namespace Ikarus::Materials
