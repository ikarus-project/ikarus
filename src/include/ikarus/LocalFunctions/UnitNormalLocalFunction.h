//
// Created by Alex on 21.04.2021.
//

#pragma once

#include "LocalFunctionInterface.h"

#include <concepts>
#include <iostream>

#include <Eigen/Core>
#include <Eigen/Dense>

#include <ikarus/LocalBasis/localBasis.h>
#include <ikarus/utils/LinearAlgebraHelper.h>

namespace Ikarus {

  template <typename DuneBasis, typename CoeffContainer>
  class UnitNormalLocalFunction : public LocalFunctionInterface<UnitNormalLocalFunction<DuneBasis, CoeffContainer>> {
    using Base = LocalFunctionInterface<UnitNormalLocalFunction<DuneBasis, CoeffContainer>>;

  public:
    friend Base;
    UnitNormalLocalFunction(const Ikarus::LocalBasis<DuneBasis>& basis_, const CoeffContainer& coeffs_)
        : basis{basis_}, coeffs{coeffs_}, coeffsAsMat{Ikarus::LinearAlgebra::viewAsEigenMatrixFixedDyn(coeffs)} {}

    using Traits = LocalFunctionTraits<UnitNormalLocalFunction>;

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
    /** \brief Type for the coordinates to store the return value */
    using GlobalE = typename FunctionReturnType::CoordinateType;
    /** \brief Type for the Jacobian matrix */
    using Jacobian = typename Traits::Jacobian;
    /** \brief Type for a column of the Jacobian matrix */
    using JacobianColType = typename Traits::JacobianColType;
    /** \brief Type for the derivatives wrt. the coeffiecients */
    using CoeffDerivMatrix = typename Traits::CoeffDerivMatrix;
    /** \brief Type for ansatz function values */
    using AnsatzFunctionType = typename Traits::AnsatzFunctionType;
    /** \brief Type for the Jacobian of the ansatz function values */
    using AnsatzFunctionJacobian = typename Traits::AnsatzFunctionJacobian;
    /** \brief Matrix to transform the ansatz function Jacobian to world coordinates*/
    using TransformMatrix = typename Traits::TransformMatrix;

    auto& coefficientsRef() { return coeffs; }

  private:
    Jacobian evaluateDerivativeWRTSpaceAllImpl(const AnsatzFunctionType& N, const AnsatzFunctionJacobian& dN) const {
      return Jacobian();
    }

    JacobianColType evaluateDerivativeWRTSpaceSingleImpl(const AnsatzFunctionType& N, const AnsatzFunctionJacobian& dN,
                                                         int spaceIndex) const {
      return JacobianColType();
    }

    CoeffDerivMatrix evaluateDerivativeWRTCoeffsImpl(const AnsatzFunctionType& N,
                                                     [[maybe_unused]] const AnsatzFunctionJacobian&,
                                                     int coeffsIndex) const {
      return CoeffDerivMatrix();
    }

    CoeffDerivMatrix evaluateSecondDerivativeWRTCoeffs(const AnsatzFunctionType& N,
                                                       [[maybe_unused]] const AnsatzFunctionJacobian&,
                                                       const AlongType& along,
                                                       const std::array<size_t, 2>& coeffsIndex) const {
      return CoeffDerivMatrix();
    }

    std::array<CoeffDerivMatrix, gridDim> evaluateDerivativeWRTCoeffsANDSpatialImpl(
        const AnsatzFunctionType& N, [[maybe_unused]] const AnsatzFunctionJacobian& dN, int coeffsIndex) const {
      return std::array<CoeffDerivMatrix, gridDim>();
    }

    CoeffDerivMatrix evaluateDerivativeWRTCoeffsANDSpatialSingleImpl(const AnsatzFunctionType& N,
                                                         [[maybe_unused]] const AnsatzFunctionJacobian& dN,
                                                         int coeffsIndex, const int spatialIndex) const {
      return CoeffDerivMatrix();
    }

    std::array<CoeffDerivMatrix, gridDim> evaluateThirdDerivativeWRTCoeffsTwoTimesAndSpatialImpl(
        const AnsatzFunctionType& N, [[maybe_unused]] const AnsatzFunctionJacobian& dN, const AlongType& along,
        const std::array<size_t, 2>& coeffsIndex) const {
      return std::array<CoeffDerivMatrix, gridDim>();
    }

    CoeffDerivMatrix evaluateThirdDerivativeWRTCoeffsTwoTimesAndSpatialSingleImpl(
        const AnsatzFunctionType& N, [[maybe_unused]] const AnsatzFunctionJacobian& dN, const AlongType& along,
        const std::array<size_t, 2>& coeffsIndex, const int spatialIndex) const {
      return CoeffDerivMatrix();
    }

    FunctionReturnType evaluateFunctionImpl(const AnsatzFunctionType& N) const {
      return (coeffsAsMat * N).normalize();
    }

    const Ikarus::LocalBasis<DuneBasis>& basis;
    CoeffContainer coeffs;
    const decltype(Ikarus::LinearAlgebra::viewAsEigenMatrixFixedDyn(coeffs)) coeffsAsMat;
  };

  template <typename DuneBasis, typename CoeffContainer>
  struct LocalFunctionTraits<UnitNormalLocalFunction<DuneBasis, CoeffContainer>> {
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
    /** \brief Type for the derivatives wrt. the coeffiecients */
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
