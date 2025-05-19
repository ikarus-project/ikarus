// SPDX-FileCopyrightText: 2021-2025 The Ikarus Developers ikarus@ibb.uni-stuttgart.de
// SPDX-License-Identifier: LGPL-3.0-or-later

/**
 * \file linearstress.hh
 * \brief Implementation of the PS function, where the linear stress is assumed.
 *
 * \ingroup ps
 */

#pragma once

#include <dune/common/fvector.hh>
#include <dune/localfefunctions/eigenDuneTransformations.hh>
#include <dune/localfefunctions/expressions/linearStrainsExpr.hh>

#include <Eigen/Core>

namespace Ikarus::PS {

struct LinearStress
{
  /**
   * \brief Compute the stress vector at a given integration point or its index.
   *
   * \param geo The geometry object providing the transposed Jacobian.
   * \param gpPos The position of the integration point.
   * \param asFunction The AssumedStress function.
   * \param beta The coefficients of the PS function.
   *
   * \tparam GEO The type of the geometry object.
   * \tparam AST The type of the Assumed Stress function
   *
   * \return The linear stress vector at the given integration point.
   */
  template <typename GEO, typename AST>
  static auto value(const GEO& geo, const Dune::FieldVector<double, GEO::mydimension>& gpPos, const AST& asFunction,
                    const auto& beta) {
    const auto P = asFunction(gpPos);
    return (P * beta).eval();
  }

  /**
   * \brief Compute the first derivative of the linear stress w.r.t beta for a given node and integration point.
   *
   * \param geo The geometry object of the finite element.
   * \param uFunction The function representing the displacement field.
   * \param localBasis The local basis object.
   * \param gpPos The position of the integration point.
   * \param gpIndex The index of the integration point.
   * \param node The FE node index (defaults to sNaN).
   * \param asFunction The AssumedStress function.
   * \param beta The coefficients of the PS function.
   *
   * \tparam GEO The type of the geometry object.
   * \tparam AST The type of the Assumed Stress function
   * function.
   *
   * \return The first derivative of the linear stress w.r.t beta.
   */
  template <typename GEO, typename AST>
  static auto firstDerivative(const GEO& geo, const auto& uFunction, const auto& localBasis, const auto& gpIndex,
                              const Dune::FieldVector<double, GEO::mydimension>& gpPos, const AST& asFunction,
                              const auto& beta, const int node = sNaN) {
    return asFunction(gpPos);
  }

  /**
   * \brief Compute the second derivative of the linear stress w.r.t beta for a given node and integration
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
   * \param asFunction The AssumedStress function.
   * \param beta The coefficients of the PS function.
   *
   * \tparam GEO The type of the geometry object.
   * \tparam ST
   * The underlying scalar type.
   * \tparam AST The type of the Assumed Stress function.
   *
   * \return The second derivative of the linear stress w.r.t beta
   */
  template <typename GEO, typename ST, typename AST>
  static auto secondDerivative(const GEO& geo, const auto& uFunction, const auto& localBasis, const auto& gpIndex,
                               const Dune::FieldVector<double, GEO::mydimension>& gpPos,
                               const Eigen::Vector<ST, GEO::mydimension*(GEO::mydimension + 1) / 2>& S,
                               const AST& asFunction, const auto& beta, const int I = sNaN, const int J = sNaN) {
    constexpr int myDim             = GEO::mydimension;
    constexpr int assumedStressSize = AST::assumedStressSize;
    return Eigen::Matrix<ST, assumedStressSize, assumedStressSize>::Zero().eval();
  }

  /** \brief The name of the assumed stress type. */
  static constexpr auto name() { return std::string("Linear Stress"); }

private:
  static constexpr int sNaN = std::numeric_limits<int>::signaling_NaN();
};

} // namespace Ikarus::PS
