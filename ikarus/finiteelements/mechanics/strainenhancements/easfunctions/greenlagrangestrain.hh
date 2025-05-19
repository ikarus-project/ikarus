// SPDX-FileCopyrightText: 2021-2025 The Ikarus Developers ikarus@ibb.uni-stuttgart.de
// SPDX-License-Identifier: LGPL-3.0-or-later

/**
 * \file greenlagrangestrain.hh
 * \brief Implementation of the EAS function, where the Green-Lagrange strain is enhanced.
 *
 * \ingroup eas
 */

#pragma once

#include <dune/common/fvector.hh>
#include <dune/localfefunctions/eigenDuneTransformations.hh>
#include <dune/localfefunctions/expressions/greenLagrangeStrains.hh>

#include <Eigen/Core>

namespace Ikarus::EAS {
/**
 * \brief A struct computing the value, first and second derivatives of the Green-Lagrange strain tensor (in Voigt
 * notation), where the Green-Lagrange strain itself is enhanced.
 */
struct GreenLagrangeStrain
{
  /**
   * \brief Compute the strain vector at a given integration point or its index.
   *
   * \param geo The geometry object providing the transposed Jacobian.
   * \param uFunction The function representing the displacement field.
   * \param gpPos The position of the integration point.
   * \param easFunction The EAS function.
   * \param alpha The coefficients of the EAS function.
   *
   * \tparam GEO The type of the geometry object.
   * \tparam EAST The type of the EAS function.
   *
   * \return The Green-Lagrange strain vector at the given integration point.
   */
  template <typename GEO, typename EAST>
  static auto value(const GEO& geo, const auto& uFunction, const Dune::FieldVector<double, GEO::mydimension>& gpPos,
                    const EAST& easFunction, const auto& alpha) {
    using namespace Dune::DerivativeDirections;
    using namespace Dune;
    constexpr int myDim        = GEO::mydimension;
    const auto& strainFunction = Dune::greenLagrangeStrains(uFunction);
    const auto E               = (strainFunction.evaluate(gpPos, on(gridElement))).eval();
    const auto Etilde          = (easFunction(gpPos) * alpha).eval();
    return (E + Etilde).eval();
  }

  /**
   * \brief Compute the first derivative of the Green-Lagrange strain w.r.t d or alpha for a given node and integration
   * point.
   *
   * \param geo The geometry object of the finite element.
   * \param uFunction The function representing the displacement field.
   * \param localBasis The local basis object.
   * \param gpPos The position of the integration point.
   * \param gpIndex The index of the integration point.
   * \param node The FE node index (defaults to sNaN).
   * \param easFunction The EAS function.
   * \param alpha The coefficients of the EAS function.
   *
   * \tparam wrtCoeff An integer indicating the coefficient w.r.t which the first derivative of Green-Lagrange strain is
   * to be returned (0 = E,d; 1 = E,a).
   * \tparam GEO The type of the geometry object.
   * \tparam EAST The type of the EAS
   * function.
   *
   * \return The first derivative of the Green-Lagrange strain w.r.t d or alpha.
   */
  template <int wrtCoeff, typename GEO, typename EAST>
  static auto firstDerivative(const GEO& geo, const auto& uFunction, const auto& localBasis, const auto& gpIndex,
                              const Dune::FieldVector<double, GEO::mydimension>& gpPos, const EAST& easFunction,
                              const auto& alpha, const int node = sNaN) {
    if constexpr (wrtCoeff == 0) {
      using namespace Dune::DerivativeDirections;
      using namespace Dune;
      constexpr int myDim        = GEO::mydimension;
      const auto& strainFunction = Dune::greenLagrangeStrains(uFunction);
      const auto bopI            = strainFunction.evaluateDerivative(gpIndex, wrt(coeff(node)), on(gridElement));
      return bopI.eval();
    } else if constexpr (wrtCoeff == 1) {
      return easFunction(gpPos);
    } else
      static_assert(Dune::AlwaysFalse<GEO>::value,
                    "firstDerivative can only be called with wrtCoeff as 0 and 1 indicating first derivative of the "
                    "Green-Lagrange strain w.r.t d and alpha, respectively.");
  }

  /**
   * \brief Compute the second derivative of the Green-Lagrange strain w.r.t d or alpha for a given node and integration
   * point.
   *
   * \param geo The geometry object of the finite element.
   * \param uFunction The function representing the displacement field.
   * \param localBasis The local basis object.
   * \param gpPos The position of the integration point.
   * \param gpIndex The index of the integration point.
   * \param I The FE node index I (defaults to sNaN).
   * \param J The FE node index J (defaults to sNaN).
   * \param S The PK2 stress (in Voigt notation).
   * \param easFunction The EAS function.
   * \param alpha The coefficients of the EAS function.
   *
   * \tparam wrtCoeff An integer indicating the coefficient w.r.t which the second derivative of Green-Lagrange strain
   * is to be returned (0 = E,dd; 1 = E,aa; 2 = E,ad = E,da).
   * \tparam GEO The type of the geometry object.
   * \tparam ST
   * The underlying scalar type.
   * \tparam EAST The type of the EAS function.
   *
   * \return The second derivative of the Green-Lagrange strain w.r.t d, alpha or the mixed derivative
   * (0 = E,dd; 1 = E,aa; 2 = E,ad = E,da).
   */
  template <int wrtCoeff, typename GEO, typename ST, typename EAST>
  static auto secondDerivative(const GEO& geo, const auto& uFunction, const auto& localBasis, const auto& gpIndex,
                               const Dune::FieldVector<double, GEO::mydimension>& gpPos,
                               const Eigen::Vector<ST, GEO::mydimension*(GEO::mydimension + 1) / 2>& S,
                               const EAST& easFunction, const auto& alpha, const int I = sNaN, const int J = sNaN) {
    constexpr int myDim              = GEO::mydimension;
    constexpr int enhancedStrainSize = EAST::enhancedStrainSize;
    if constexpr (wrtCoeff == 0) {
      using namespace Dune::DerivativeDirections;
      using namespace Dune;
      const auto& strainFunction = Dune::greenLagrangeStrains(uFunction);
      const auto kg = strainFunction.evaluateDerivative(gpIndex, wrt(coeff(I, J)), along(S), on(gridElement));
      return kg.eval(); // E_dd
    } else if constexpr (wrtCoeff == 1) {
      const auto kg = Eigen::Matrix<ST, enhancedStrainSize, enhancedStrainSize>::Zero().eval();
      return kg; // E_aa
    } else if constexpr (wrtCoeff == 2) {
      const auto kg = Eigen::Matrix<ST, enhancedStrainSize, myDim>::Zero().eval();
      return kg; // E_ad
    } else
      static_assert(
          Dune::AlwaysFalse<GEO>::value,
          "secondDerivative can only be called with wrtCoeff as 0, 1 and 2 indicating second derivative of the "
          "Green-Lagrange strain w.r.t d and alpha and the mixed derivative, respectively.");
  }

  /** \brief The name of the strain measure enhanced w.r.t EAS method. */
  static constexpr auto name() { return std::string("Green-Lagrange Strain"); }

private:
  static constexpr int sNaN = std::numeric_limits<int>::signaling_NaN();
};

} // namespace Ikarus::EAS
