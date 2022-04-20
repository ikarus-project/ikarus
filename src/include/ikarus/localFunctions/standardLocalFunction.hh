//
// Created by Alex on 21.04.2021.
//

#pragma once

#include "localFunctionInterface.hh"
#include "localFunctionExpression.h"

#include <concepts>
#include <iostream>

#include <Eigen/Core>
#include <Eigen/Dense>

#include <ikarus/localBasis/localBasis.hh>
#include <ikarus/utils/linearAlgebraHelper.hh>

namespace Ikarus {

  template <typename DuneBasis, typename CoeffContainer>
  class StandardLocalFunction : public LocalFunctionExpression<StandardLocalFunction<DuneBasis, CoeffContainer>> {
    using Base = LocalFunctionExpression<StandardLocalFunction<DuneBasis, CoeffContainer>>;

  public:
    friend Base;
    StandardLocalFunction(const Ikarus::LocalBasis<DuneBasis>& p_basis, const CoeffContainer& coeffs_)
        : basis_{p_basis}, coeffs{coeffs_}, coeffsAsMat{Ikarus::LinearAlgebra::viewAsEigenMatrixFixedDyn(coeffs)} {}

    using Traits = LocalFunctionTraits<StandardLocalFunction>;
    /** \brief Type used for coordinates */
    using ctype = typename Traits::ctype;
    //    /** \brief Dimension of the coeffs */
    static constexpr int valueSize = Traits::valueSize;
    /** \brief Dimension of the grid */
    static constexpr int gridDim = Traits::gridDim;
    /** \brief Type for coordinate vector in world space */
    using FunctionReturnType = typename Traits::FunctionReturnType;
    /** \brief Type for the directional derivatives */
    using AlongType = typename Traits::AlongType;
    /** \brief Type for the coordinates to store the return value */
    using GlobalE = typename FunctionReturnType::CoordinateType;
    /** \brief Type for the Jacobian matrix */
    using Jacobian = typename Traits::Jacobian;
    /** \brief Type for a column of the Jacobian matrix */
    using JacobianColType = typename Traits::JacobianColType;
    /** \brief Type for the derivatives wrT the coeffiecients */
    using CoeffDerivMatrix = typename Traits::CoeffDerivMatrix;
    /** \brief Type for ansatz function values */
    using AnsatzFunctionType = typename Traits::AnsatzFunctionType;
    /** \brief Type for the Jacobian of the ansatz function values */
    using AnsatzFunctionJacobian = typename Traits::AnsatzFunctionJacobian;

    auto& coefficientsRef() { return coeffs; }

    const Ikarus::LocalBasis<DuneBasis>& basis() const
    {
      return basis_;
    }

  private:
    Jacobian evaluateDerivativeWRTSpaceAllImpl(const AnsatzFunctionType& N, const AnsatzFunctionJacobian& dN) const {
      return coeffsAsMat
             * dN.template cast<ctype>();  // The cast here is only necessary since the autodiff types are not working
                                           // otherwise, see Issue https://github.com/autodiff/autodiff/issues/73
    }

    JacobianColType evaluateDerivativeWRTSpaceSingleImpl(const AnsatzFunctionType&, const AnsatzFunctionJacobian& dN,
                                                         int spaceIndex) const {
      return coeffsAsMat * dN.col(spaceIndex);
    }

    CoeffDerivMatrix evaluateDerivativeWRTCoeffsImpl(const AnsatzFunctionType& N,
                                                     [[maybe_unused]] const AnsatzFunctionJacobian&,
                                                     int coeffsIndex) const {
      CoeffDerivMatrix mat;
      mat.setIdentity(valueSize);
      mat.diagonal() *= N[coeffsIndex];
      return mat;
    }

    std::array<CoeffDerivMatrix, gridDim> evaluateDerivativeWRTCoeffsANDSpatialImpl(
        const AnsatzFunctionType& N, [[maybe_unused]] const AnsatzFunctionJacobian& dN, int coeffsIndex) const {
      std::array<CoeffDerivMatrix, gridDim> Warray;
      for (int dir = 0; dir < gridDim; ++dir) {
        Warray[dir].setIdentity(valueSize);
        Warray[dir].diagonal() *= dN(coeffsIndex, dir);
      }

      return Warray;
    }

    CoeffDerivMatrix evaluateDerivativeWRTCoeffsANDSpatialSingleImpl(const AnsatzFunctionType& N,
                                                                     [[maybe_unused]] const AnsatzFunctionJacobian& dN,
                                                                     int coeffsIndex, const int spatialIndex) const {
      CoeffDerivMatrix W;
      W.setIdentity(valueSize);
      W.diagonal() *= dN(coeffsIndex, spatialIndex);

      return W;
    }

    FunctionReturnType evaluateFunctionImpl(const AnsatzFunctionType& N) const {
      return FunctionReturnType(coeffsAsMat * N);
    }


    const Ikarus::LocalBasis<DuneBasis>& basis_;
    CoeffContainer coeffs;
    const decltype(Ikarus::LinearAlgebra::viewAsEigenMatrixFixedDyn(coeffs)) coeffsAsMat;
  };

  template <typename DuneBasis, typename CoeffContainer>
  struct LocalFunctionTraits<StandardLocalFunction<DuneBasis, CoeffContainer>> {
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
    using CoeffDerivMatrix = Eigen::DiagonalMatrix<ctype, valueSize>;
    /** \brief Type for the Jacobian of the ansatz function values */
    using AnsatzFunctionJacobian = typename Ikarus::LocalBasis<DuneBasis>::JacobianType;
    /** \brief Type for ansatz function values */
    using AnsatzFunctionType = typename Ikarus::LocalBasis<DuneBasis>::AnsatzFunctionType;
    /** \brief Type for the points for evaluation, usually the integration points */
    using DomainType = typename DuneBasis::Traits::DomainType;
    /** \brief Type for a column of the Jacobian matrix */
    using JacobianColType = typename Eigen::internal::plain_col_type<Jacobian>::type;
    /** \brief Type for the directional derivatives */
    using AlongType = Eigen::Vector<ctype, valueSize>;
  };

}  // namespace Ikarus
