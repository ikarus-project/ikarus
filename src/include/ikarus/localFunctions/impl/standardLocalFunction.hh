/*
 * This file is part of the Ikarus distribution (https://github.com/IkarusRepo/Ikarus).
 * Copyright (c) 2022. The Ikarus developers.
 *
 * This library is free software; you can redistribute it and/or
 *  modify it under the terms of the GNU Lesser General Public
 *  License as published by the Free Software Foundation; either
 *  version 2.1 of the License, or (at your option) any later version.
 *
 * This library is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
 * Lesser General Public License for more details.
 *
 * You should have received a copy of the GNU Lesser General Public
 * License along with this library; if not, write to the Free Software
 * Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301 USA
 */

#pragma once

#include "clonableLocalFunction.hh"

#include <concepts>

#include <dune/common/indices.hh>

#include <Eigen/Core>
#include <Eigen/Dense>

#include <ikarus/localBasis/localBasis.hh>
#include <ikarus/localFunctions/localFunctionHelper.hh>
#include <ikarus/localFunctions/localFunctionInterface.hh>
#include <ikarus/utils/linearAlgebraHelper.hh>

namespace Ikarus {

  template <typename DuneBasis, typename CoeffContainer, std::size_t ID = 0>
  class StandardLocalFunction : public LocalFunctionInterface<StandardLocalFunction<DuneBasis, CoeffContainer, ID>>,
                                public ClonableLocalFunction<StandardLocalFunction<DuneBasis, CoeffContainer, ID>> {
    using Interface = LocalFunctionInterface<StandardLocalFunction<DuneBasis, CoeffContainer, ID>>;

  public:
    friend Interface;
    friend ClonableLocalFunction<StandardLocalFunction>;

    constexpr StandardLocalFunction(const Ikarus::LocalBasis<DuneBasis>& p_basis, const CoeffContainer& coeffs_,
                                    Dune::template index_constant<ID> = Dune::template index_constant<std::size_t(0)>{})
        : basis_{p_basis}, coeffs{coeffs_}, coeffsAsMat{Ikarus::viewAsEigenMatrixFixedDyn(coeffs)} {}

    static constexpr bool isLeaf = true;
    using Ids                    = Dune::index_constant<ID>;

    template <size_t ID_ = 0>
    static constexpr int orderID = ID_ == ID ? linear : constant;

    template <typename LocalFunctionEvaluationArgs_, typename LocalFunctionImpl_>
    friend auto evaluateDerivativeImpl(const LocalFunctionInterface<LocalFunctionImpl_>& f,
                                       const LocalFunctionEvaluationArgs_& localFunctionArgs);

    template <typename LocalFunctionEvaluationArgs_, typename LocalFunctionImpl_>
    friend auto evaluateFunctionImpl(const LocalFunctionInterface<LocalFunctionImpl_>& f,
                                     const LocalFunctionEvaluationArgs_& localFunctionArgs);

    using Traits = LocalFunctionTraits<StandardLocalFunction<DuneBasis, CoeffContainer, ID>>;
    /** \brief Type used for coordinates */
    using ctype = typename Traits::ctype;
    //    /** \brief Dimension of the coeffs */
    static constexpr int valueSize = Traits::valueSize;
    /** \brief Dimension of the correction size of coeffs */
    static constexpr int correctionSize = Traits::correctionSize;
    /** \brief Dimension of the grid */
    static constexpr int gridDim = Traits::gridDim;
    /** \brief Type for coordinate vector in world space */
    using FunctionReturnType = typename Traits::FunctionReturnType;
    /** \brief The manifold where the function values lives in */
    using Manifold = typename Traits::Manifold;
    /** \brief Type for the directional derivatives */
    using AlongType = typename Traits::AlongType;
    /** \brief Type for the Jacobian matrix */
    using Jacobian = typename Traits::Jacobian;
    /** \brief Type for a column of the Jacobian matrix */
    using JacobianColType = typename Traits::JacobianColType;
    /** \brief Type for the derivatives wrT the coefficients */
    using CoeffDerivMatrix = typename Traits::CoeffDerivMatrix;
    /** \brief Type for ansatz function values */
    using AnsatzFunctionType = typename Traits::AnsatzFunctionType;
    /** \brief Type for the Jacobian of the ansatz function values */
    using AnsatzFunctionJacobian = typename Traits::AnsatzFunctionJacobian;

    const auto& coefficientsRef() const { return coeffs; }
    auto& coefficientsRef() { return coeffs; }

    template <typename OtherType>
    struct Rebind {
      using other = StandardLocalFunction<
          DuneBasis, typename Std::Rebind<CoeffContainer, typename Manifold::template Rebind<OtherType>::other>::other,
          ID>;
    };

    const Ikarus::LocalBasis<DuneBasis>& basis() const { return basis_; }

  private:
    template <typename DomainTypeOrIntegrationPointIndex, typename... TransformArgs>
    FunctionReturnType evaluateFunctionImpl(const DomainTypeOrIntegrationPointIndex& ipIndexOrPosition,
                                            [[maybe_unused]] const TransformWith<TransformArgs...>&) const {
      const auto& N = evaluateFunctionWithIPorCoord(ipIndexOrPosition, basis_);

      return FunctionReturnType(coeffsAsMat * N);
    }

    template <typename DomainTypeOrIntegrationPointIndex, typename... TransformArgs>
    Jacobian evaluateDerivativeWRTSpaceAllImpl(const DomainTypeOrIntegrationPointIndex& ipIndexOrPosition,
                                               const TransformWith<TransformArgs...>& transArgs) const {
      const auto& dNraw = evaluateDerivativeWithIPorCoord(ipIndexOrPosition, basis_);
      maytransformDerivatives(dNraw, dNTransformed, transArgs);
      return coeffsAsMat
             * dNTransformed
                   .template cast<ctype>();  // The cast here is only necessary since the autodiff types are not working
                                             // otherwise, see Issue https://github.com/autodiff/autodiff/issues/73
    }

    template <typename DomainTypeOrIntegrationPointIndex, typename... TransformArgs>
    JacobianColType evaluateDerivativeWRTSpaceSingleImpl(const DomainTypeOrIntegrationPointIndex& ipIndexOrPosition,
                                                         int spaceIndex,
                                                         const TransformWith<TransformArgs...>& transArgs) const {
      const auto& dNraw = evaluateDerivativeWithIPorCoord(ipIndexOrPosition, basis_);
      maytransformDerivatives(dNraw, dNTransformed, transArgs);

      return coeffsAsMat * dNTransformed.col(spaceIndex);
    }

    template <typename DomainTypeOrIntegrationPointIndex, typename... TransformArgs>
    CoeffDerivMatrix evaluateDerivativeWRTCoeffsImpl(const DomainTypeOrIntegrationPointIndex& ipIndexOrPosition,
                                                     int coeffsIndex, const TransformWith<TransformArgs...>&) const {
      const auto& N = evaluateFunctionWithIPorCoord(ipIndexOrPosition, basis_);
      CoeffDerivMatrix mat;
      mat.setIdentity(valueSize);
      mat.diagonal() *= N[coeffsIndex];
      return mat;
    }

    template <typename DomainTypeOrIntegrationPointIndex, typename... TransformArgs>
    std::array<CoeffDerivMatrix, gridDim> evaluateDerivativeWRTCoeffsANDSpatialImpl(
        const DomainTypeOrIntegrationPointIndex& ipIndexOrPosition, int coeffsIndex,
        const TransformWith<TransformArgs...>& transArgs) const {
      const auto& dNraw = evaluateDerivativeWithIPorCoord(ipIndexOrPosition, basis_);
      maytransformDerivatives(dNraw, dNTransformed, transArgs);
      std::array<CoeffDerivMatrix, gridDim> Warray;
      for (int dir = 0; dir < gridDim; ++dir) {
        Warray[dir].setIdentity(valueSize);
        Warray[dir].diagonal() *= dNTransformed(coeffsIndex, dir);
      }

      return Warray;
    }

    template <typename DomainTypeOrIntegrationPointIndex, typename... TransformArgs>
    CoeffDerivMatrix evaluateDerivativeWRTCoeffsANDSpatialSingleImpl(
        const DomainTypeOrIntegrationPointIndex& ipIndexOrPosition, int coeffsIndex, int spatialIndex,
        const TransformWith<TransformArgs...>& transArgs) const {
      const auto& dNraw = evaluateDerivativeWithIPorCoord(ipIndexOrPosition, basis_);
      maytransformDerivatives(dNraw, dNTransformed, transArgs);
      CoeffDerivMatrix W;
      W.setIdentity(valueSize);
      W.diagonal() *= dNTransformed(coeffsIndex, spatialIndex);

      return W;
    }

    mutable AnsatzFunctionJacobian dNTransformed;
    const Ikarus::LocalBasis<DuneBasis>& basis_;
    CoeffContainer coeffs;
    const decltype(Ikarus::viewAsEigenMatrixFixedDyn(coeffs)) coeffsAsMat;
  };

  template <typename DuneBasis, typename CoeffContainer, std::size_t ID>
  struct LocalFunctionTraits<StandardLocalFunction<DuneBasis, CoeffContainer, ID>> {
    /** \brief Type used for coordinates */
    using ctype = typename CoeffContainer::value_type::ctype;
    /** \brief Dimension of the coeffs */
    static constexpr int valueSize = CoeffContainer::value_type::valueSize;
    /** \brief Dimension of the correction size of coeffs */
    static constexpr int correctionSize = CoeffContainer::value_type::correctionSize;
    /** \brief Dimension of the grid */
    static constexpr int gridDim = Ikarus::LocalBasis<DuneBasis>::gridDim;
    /** \brief The manifold where the function values lives in */
    using Manifold = typename CoeffContainer::value_type;
    /** \brief Type for the return value */
    using FunctionReturnType = typename Manifold::CoordinateType;
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
