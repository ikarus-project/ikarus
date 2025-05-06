// SPDX-FileCopyrightText: 2021-2025 The Ikarus Developers ikarus@ibb.uni-stuttgart.de
// SPDX-License-Identifier: LGPL-3.0-or-later

#pragma once

/**
 * \file materialhelpers.hh
 * \brief helper functions used by material model implementations.
 * \ingroup  materials
 */

#include <ranges>

#include <dune/common/float_cmp.hh>

#include <Eigen/Dense>

#include <ikarus/utils/concepts.hh>
#include <ikarus/utils/tensorutils.hh>

namespace Ikarus::Materials {
/**
 * \brief Represents a pair of stress or strain matrix indices (row and column).
 */
struct MatrixIndexPair
{
  Eigen::Index row; ///< Row index.
  Eigen::Index col; ///< Column index.
};
} // namespace Ikarus::Materials
namespace Ikarus::Impl {

/**
 * \brief Helper function to create an array of free Voigt indices.
 * \tparam size The size of the fixed pairs array.
 * \param fixed An array of MatrixIndexPair representing fixed indices.
 * \return std::array<size_t, 6 - size> The array of free Voigt indices.
 */
template <size_t size>
consteval auto createfreeVoigtIndices(const std::array<Materials::MatrixIndexPair, size>& fixed) {
  std::array<size_t, 6 - size> res{};
  std::array<size_t, size> voigtFixedIndices;
  std::ranges::transform(fixed, voigtFixedIndices.begin(), [](auto pair) { return toVoigt(pair.row, pair.col); });
  std::ranges::sort(voigtFixedIndices);
  std::ranges::set_difference(std::ranges::iota_view(size_t(0), size_t(6)), voigtFixedIndices, res.begin());
  std::ranges::sort(res);
  return res;
}

/**
 * \brief Helper function to create an array of fixed Voigt indices.
 * \tparam size The size of the fixed pairs array.
 * \param fixed An array of MatrixIndexPair representing fixed indices.
 * \return std::array<size_t, size> The array of fixed Voigt indices.
 */
template <size_t size>
consteval auto createFixedVoigtIndices(const std::array<Materials::MatrixIndexPair, size>& fixed) {
  std::array<size_t, size> fixedIndices;
  std::ranges::transform(fixed, fixedIndices.begin(), [](auto pair) { return toVoigt(pair.row, pair.col); });
  std::ranges::sort(fixedIndices);
  return fixedIndices;
}

/**
 * \brief Helper function to count the number of diagonal indices in the fixed pairs array.
 * \tparam size The size of the fixed pairs array.
 * \param fixed An array of MatrixIndexPair representing fixed indices.
 * \return constexpr size_t The number of diagonal indices.
 */
template <size_t size>
constexpr size_t countDiagonalIndices(const std::array<Materials::MatrixIndexPair, size>& fixed) {
  size_t count = 0;
  for (auto v : fixed) {
    if (v.col == v.row)
      ++count;
  }
  return count;
}

/**
 * \brief Converts the input measure from voigt notation to tensor notation if its not already in tensor notation.
 *
 * \tparam Derived The derived type of the input matrix.
 * \param v The input measure.
 * \param isStrain Flag indicating whether the vector represents a strain (default is true).
 * \return decltype(auto) The converted matrix.
 */
template <typename Derived>
auto maybeFromVoigt(const Eigen::MatrixBase<Derived>& v, bool isStrain = true) {
  if constexpr (Concepts::EigenVector<Derived>) { // receiving vector means Voigt notation
    return fromVoigt(v.derived(), isStrain).eval();
  } else
    return v.derived().eval();
}

/**
 * \brief Converts the input measure from tensor notation to voigt notation if its not already in voigt notation
 *
 * \tparam Derived The derived type of the input matrix.
 * \param v the input measure
 * \param isStrain Flag indicating whether the vector represents a strain (default is true).
 * \return decltype(auto) The converted vector.
 */
template <typename Derived>
auto maybeToVoigt(const Eigen::MatrixBase<Derived>& v, bool isStrain = false) {
  if constexpr (Concepts::EigenVector<Derived>) { // receiving vector means Voigt notation
    return v.derived().eval();
  } else
    return toVoigt(v.derived(), isStrain).eval();
}

/**
 * \brief Checks weather a scalar is greater then zero. Here its used in the context of checking the determinant of the
 * right Cauchy Green tensor C.
 *
 * \param eps the epislon for comparison (defaults to 1e-10)
 * \param det the evaluated determinant
 */
template <typename ScalarType>
void checkPositiveOrAbort(ScalarType det, double eps = 1e-10) {
  if (Dune::FloatCmp::le(static_cast<double>(det), 0.0, eps)) {
    std::cerr << "Determinant of right Cauchy Green tensor C must be greater than zero. detC = " +
                     std::to_string(static_cast<double>(det));
    abort();
  }
}

/**
 * \brief Computes the principal stretches of the input strain matrix C.
 *
 * \tparam ScalarType The ScalarType
 * \tparam Derived The derived type of the input matrix
 * \param C the input strain matrix
 * \param options should be either `Eigen::ComputeEigenvectors` or `Eigen::EigenvaluesOnly`
 * \return auto pair of principalstretches and corresponding eigenvectors (if `Eigen::EigenvaluesOnly` the
 * eigenvectors is a zero matrix )
 */
template <typename ScalarType, typename Derived>
auto principalStretches(const Eigen::MatrixBase<Derived>& C, int options = Eigen::ComputeEigenvectors) {
  Eigen::SelfAdjointEigenSolver<Derived> eigensolver{};

  eigensolver.compute(C, options);

  if (eigensolver.info() != Eigen::Success)
    DUNE_THROW(Dune::MathError, "Failed to compute eigenvalues and eigenvectors of C.");

  auto& eigenvalues  = eigensolver.eigenvalues();
  auto& eigenvectors = options == Eigen::ComputeEigenvectors ? eigensolver.eigenvectors() : Derived::Zero();

  auto principalStretches = eigenvalues.cwiseSqrt().eval();
  return std::make_pair(principalStretches, eigenvectors);
}

/**
 * \brief Computes the determinant of a matrix through its principal values (i.e. eigenvalues).
 *
 * \tparam Vector the type of the vector of principal stretches
 * \param principalValues the principal values.
 * \return auto The determinant.
 */
template <Concepts::EigenVector3 Vector>
inline Vector::Scalar determinantFromPrincipalValues(const Vector& principalValues) {
  return principalValues.prod();
}

/**
 * \brief Computes the deviatoric part of the principal stretches as \f$ \bar{\lambda_i} = \lambda_i^{-\frac{1}{3}}
 * \f$
 *
 * \tparam Vector the type of the vector of principal stretches
 * \param lambda the total principal stretches
 * \return Vector the deviatoric principal stretches
 */
template <Concepts::EigenVector3 Vector>
inline Vector deviatoricStretches(const Vector& lambda) {
  using ScalarType = typename Vector::Scalar;
  ScalarType J     = determinantFromPrincipalValues(lambda);
  ScalarType Jmod  = pow(J, -1.0 / 3.0);
  return Jmod * lambda;
}

/**
 * \brief Computes the invariants from the principal stretches
 *
 * \tparam Vector the type of the vector of principal stretches
 * \param lambda the total principal stretches
 * \return Vector the invariants
 */
template <Concepts::EigenVector3 Vector>
inline Vector invariants(const Vector& lambda) {
  using ScalarType   = typename Vector::Scalar;
  auto lambdaSquared = lambda.array().square();
  auto invariants    = Vector::Zero().eval();

  invariants[0] = lambdaSquared.sum();
  invariants[1] =
      lambdaSquared[0] * lambdaSquared[1] + lambdaSquared[1] * lambdaSquared[2] + lambdaSquared[0] * lambdaSquared[2];
  invariants[2] = determinantFromPrincipalValues(lambdaSquared);

  return invariants;
}
} // namespace Ikarus::Impl