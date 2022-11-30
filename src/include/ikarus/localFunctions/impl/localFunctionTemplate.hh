// SPDX-FileCopyrightText: 2022 The Ikarus Developers mueller@ibb.uni-stuttgart.de
// SPDX-License-Identifier: LGPL-2.1-or-later

#pragma once

#include "localFunctionInterface.h"

#include <concepts>

#include <Eigen/Core>
#include <Eigen/Dense>

#include <ikarus/LocalBasis/localBasis.hh>
#include <ikarus/utils/LinearAlgebraHelper.hh>

namespace Ikarus {

  template <typename DuneBasis, typename CoeffContainer>
  class LocalFunctionTemplate : public LocalFunctionInterface<LocalFunctionTemplate<DuneBasis, CoeffContainer>> {
    using Base = LocalFunctionInterface<LocalFunctionTemplate<DuneBasis, CoeffContainer>>;

  public:
    friend Base;
    LocalFunctionTemplate(const Ikarus::LocalBasis<DuneBasis>& basis_, const CoeffContainer& coeffs_)
        : basis{basis_}, coeffs{coeffs_}, coeffsAsMat{Ikarus:: ::viewAsEigenMatrixFixedDyn(coeffs)} {}

    using Traits = LocalFunctionTraits<LocalFunctionTemplate>;

    /** \brief Type used for coordinates */
    using ctype = typename Traits::ctype;
    /** \brief Dimension of the coeffs */
    static constexpr int valueSize = Traits::valueSize;
    /** \brief Dimension of the grid */
    static constexpr int gridDim = Traits::gridDim;
    /** \brief Type for the return value */
    using FunctionReturnType = typename Traits::FunctionReturnType;
    /** \brief Type for the directional derivatives */
    using AlongType = typename Traits::AlongType;
    ;
    /** \brief Type for the coordinates to store the return value */
    using GlobalE = typename FunctionReturnType::CoordinateType;
    /** \brief Type for the Jacobian matrix */
    using Jacobian = typename Traits::Jacobian;
    /** \brief Type for a column of the Jacobian matrix */
    using JacobianColType = typename Traits::JacobianColType;
    /** \brief Type for the derivatives wrt. the coefficients */
    using CoeffDerivMatrix = typename Traits::CoeffDerivMatrix;
    /** \brief Type for ansatz function values */
    using AnsatzFunctionType = typename Traits::AnsatzFunctionType;
    /** \brief Type for the Jacobian of the ansatz function values */
    using AnsatzFunctionJacobian = typename Traits::AnsatzFunctionJacobian;
    /** \brief Matrix to transform the ansatz function Jacobian to world coordinates*/
    using TransformMatrix = typename Traits::TransformMatrix;

    auto& coefficientsRef() { return coeffs; }

  private:
    FunctionReturnType evaluateFunctionImpl(const Eigen::VectorXd& N) const { return FunctionReturnType{}; }

    Jacobian evaluateDerivativeWRTSpaceAllImpl(const AnsatzFunctionType& N, const AnsatzFunctionJacobian& dN) const {
      return Jacobian{};
    }

    JacobianColType evaluateDerivativeWRTSpaceSingleImpl(const AnsatzFunctionType& N, const AnsatzFunctionJacobian& dN,
                                                         int spaceIndex) const {
      return JacobianColType{};
    }

    CoeffDerivMatrix evaluateDerivativeWRTCoeffsImpl(const AnsatzFunctionType& N,
                                                     [[maybe_unused]] const AnsatzFunctionJacobian&,
                                                     int coeffsIndex) const {
      return CoeffDerivMatrix{};
    }

    CoeffDerivMatrix evaluateSecondDerivativeWRTCoeffs(const AnsatzFunctionType& N,
                                                       [[maybe_unused]] const AnsatzFunctionJacobian&,
                                                       const AlongType& along,
                                                       const std::array<size_t, 2>& coeffsIndex) const {
      return CoeffDerivMatrix{};
    }

    std::array<CoeffDerivMatrix, gridDim> evaluateDerivativeWRTCoeffsANDSpatialImpl(
        const AnsatzFunctionType& N, [[maybe_unused]] const AnsatzFunctionJacobian& dN, int coeffsIndex) const {
      return std::array<CoeffDerivMatrix, gridDim>{};
    }

    CoeffDerivMatrix evaluateDerivativeWRTCoeffsANDSpatialSingleImpl(const AnsatzFunctionType& N,
                                                                     [[maybe_unused]] const AnsatzFunctionJacobian& dN,
                                                                     int coeffsIndex, const int spatialIndex) const {
      return CoeffDerivMatrix{};
    }

    CoeffDerivMatrix evaluateThirdDerivativeWRTCoeffsTwoTimesAndSpatialImpl(
        const AnsatzFunctionType& N, [[maybe_unused]] const AnsatzFunctionJacobian& dN, const AlongType& along,
        const std::array<size_t, 2>& coeffsIndex) const {
      return CoeffDerivMatrix{};
    }

    CoeffDerivMatrix evaluateThirdDerivativeWRTCoeffsTwoTimesAndSpatialSingleImpl(
        const AnsatzFunctionType& N, [[maybe_unused]] const AnsatzFunctionJacobian& dN, const AlongType& along,
        const std::array<size_t, 2>& coeffsIndex, const int spatialIndex) const {
      return CoeffDerivMatrix{};
    }

    const Ikarus::LocalBasis<DuneBasis>& basis;
    CoeffContainer coeffs;
    const decltype(Ikarus:: ::viewAsEigenMatrixFixedDyn(coeffs)) coeffsAsMat;
  };

  template <typename DuneBasis, typename CoeffContainer>
  struct LocalFunctionTraits<LocalFunctionTemplate<DuneBasis, CoeffContainer>> {
    /** \brief Type used for coordinates */
    using ctype = typename CoeffContainer::value_type::ctype;
    /** \brief Dimension of the coeffs */
    static constexpr int valueSize = CoeffContainer::value_type::valueSize;
    /** \brief Dimension of the grid */
    static constexpr int gridDim = Ikarus::LocalBasis<DuneBasis>::gridDim;
    /** \brief Type for the return value */
    using FunctionReturnType = typename CoeffContainer::value_type;
    /** \brief Type for the Jacobian matrix */
    using Jacobian = Eigen::Matrix<ctype, valueSize, gridDim>;
    /** \brief Type for the derivatives wrt. the coefficients */
    using CoeffDerivMatrix = Eigen::Matrix<ctype, valueSize, valueSize>;
    /** \brief Type for the Jacobian of the ansatz function values */
    using AnsatzFunctionJacobian = typename Ikarus::LocalBasis<DuneBasis>::JacobianType;
    /** \brief Type for ansatz function values */
    using AnsatzFunctionType = typename Ikarus::LocalBasis<DuneBasis>::AnsatzFunctionType;
    /** \brief Type for the points for evaluation, usually the integration points */
    using DomainType = typename DuneBasis::Traits::DomainType;
    /** \brief Matrix to transform the ansatz function Jacobian to world coordinates*/
    using TransformMatrix = Eigen::Matrix<ctype, gridDim, gridDim>;
    /** \brief Type for a column of the Jacobian matrix */
    using JacobianColType = typename Eigen::internal::plain_col_type<Jacobian>::type;
    /** \brief Type for the directional derivatives */
    using AlongType = Eigen::Vector<ctype, valueSize>;
  };

}  // namespace Ikarus
