// SPDX-FileCopyrightText: 2021-2025 The Ikarus Developers ikarus@ibb.uni-stuttgart.de
// SPDX-License-Identifier: LGPL-3.0-or-later

/**
 * \file interface.hh
 * \brief Implementation of the interface for the deviatoric part of a hyperelastic material.
 * \ingroup  materials
 */

#pragma once

#include <dune/common/float_cmp.hh>

#include <ikarus/finiteelements/mechanics/materials/hyperelastic/concepts.hh>
#include <ikarus/utils/tensorutils.hh>

namespace Ikarus::Materials {

/**
 * \brief This is the interface implementation for the deviatoric part of a hyperelastic material.
 *    It is intended to be used with the hyperelastic material model.
 * \ingroup materials
 *
 * \details The deviatoric part of the hyperelastic model, i.e., related to
 * \f$ \hat{\Psi}(\la_1, \la_2, \la_3) \f$, is parametrized with a certain
 * deviatoric function (DF) implemented in terms of principal stretches.
 * The three interface functions (energy, streses and tangentModulus) are called with the argument being the principal
 * stretches (\f$ \la_i \f$). The underlying deviatoric function must only implement the energy
 * \f$ \hat{\Psi}(\la_1, \la_2, \la_3) \f$ and its first and second derivatives
 * w.r.t the total principal stretches.
 *
 * \tparam DF Deviatoric function.
 */
template <Concepts::DeviatoricFunction DF>
struct Deviatoric
{
  using ScalarType = typename DF::ScalarType;

  template <typename ST = ScalarType>
  using PrincipalStretches = typename DF::template PrincipalStretches<ST>;
  using MaterialParameters = typename DF::MaterialParameters;

  using DeviatoricFunction = DF;

  template <typename ST = ScalarType>
  using FirstDerivative = typename DF::template FirstDerivative<ST>;
  template <typename ST = ScalarType>
  using SecondDerivative = typename DF::template SecondDerivative<ST>;

  static constexpr int dim = 3;

  template <typename ST = ScalarType>
  using StressMatrix = Eigen::Vector<ST, dim>;

  template <typename ST = ScalarType>
  using MaterialTensor = Eigen::TensorFixedSize<ST, Eigen::Sizes<dim, dim, dim, dim>>;

  [[nodiscard]] constexpr static std::string name() noexcept { return "Deviatoric function: " + DF::name(); }

  Deviatoric(const DF df)
      : deviatoricFunction_{df} {}

  /**
   * \brief Returns the material parameters stored in the deviatoric part of the material.
   */
  const MaterialParameters materialParameters() const { return deviatoricFunction_.materialParametersImpl(); }

  /**
   * \brief Returns the stored energy obtained from the deviatoric function.
   *
   * \param lambdas the principal stretches.
   * \tparam ST the scalartype of the principal stretches
   * \return ScalarType the energy.
   */

  template <typename ST>
  ST storedEnergy(const PrincipalStretches<ST>& lambda) const {
    return deviatoricFunction_.storedEnergyImpl(lambda);
  };

  /**
   * \brief Returns the principal PK2 stresses obtained from the first derivative of the deviatoric function.
   *
   * \param lambda the principal stretches.
   * \tparam ST the scalartype of the principal stretches
   * \return StressMatrix the stresses in principal strains coordinate system.
   */
  template <typename ST>
  StressMatrix<ST> stresses(const PrincipalStretches<ST>& lambda) const {
    auto dWdLambda = deviatoricFunction_.firstDerivativeImpl(lambda);
    return (dWdLambda.array() / lambda.array()).eval();
  }

  /**
   * \brief Returns the material tangent modulus obtained from the second derivative of the deviatoric function.
   *
   * \param lambda the principal stretches.
   * \tparam ST the scalartype of the principal stretches
   * \return MaterialTensor the tangentModuli in principal strains coordinate system.
   */

  template <typename ST>
  MaterialTensor<ST> tangentModuli(const PrincipalStretches<ST>& lambda) const {
    auto S  = stresses(lambda);
    auto dS = deviatoricFunction_.secondDerivativeImpl(lambda);

    auto L = MaterialTensor<ST>{};
    L.setZero();

    for (const auto i : dimensionRange())
      for (const auto k : dimensionRange())
        L(i, i, k, k) = 1.0 / (lambda(i) * lambda(k)) * dS(i, k);

    Eigen::ArrayXd lambdaSquared = lambda.array().square();
    for (const auto i : dimensionRange())
      for (const auto k : dimensionRange())
        if (i != k) {
          if (Dune::FloatCmp::eq(lambda(i), lambda(k), 1e-8))
            L(i, k, i, k) = 0.5 * (L(i, i, i, i) - L(i, i, k, k));
          else
            L(i, k, i, k) += (S(i) - S(k)) / (lambdaSquared(i) - lambdaSquared(k));
        }

    return L;
  };

  /**
   * \brief Rebinds the material to a different scalar type.
   * \tparam STO The target scalar type.
   * \return The rebound deviatoric part.
   */
  template <typename STO>
  auto rebind() const {
    auto reboundDF = deviatoricFunction_.template rebind<STO>();
    return Deviatoric<decltype(reboundDF)>{reboundDF};
  }

private:
  DF deviatoricFunction_;

  inline static constexpr auto dimensionRange() { return Dune::range(dim); }
};
} // namespace Ikarus::Materials