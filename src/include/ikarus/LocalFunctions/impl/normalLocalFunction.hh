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

  template <typename Geometry, typename DuneBasis, typename CoeffContainer, std::size_t ID = 0>
  class NormalLocalFunction : public LocalFunctionInterface<NormalLocalFunction<Geometry,DuneBasis, CoeffContainer, ID>>,
                                public ClonableLocalFunction<NormalLocalFunction<Geometry,DuneBasis, CoeffContainer, ID>> {
    using Interface = LocalFunctionInterface<NormalLocalFunction<Geometry,DuneBasis, CoeffContainer, ID>>;

  public:
    friend Interface;
    friend ClonableLocalFunction<NormalLocalFunction>;

    constexpr NormalLocalFunction(const Geometry& geo,const Ikarus::LocalBasis<DuneBasis>& p_basis, const CoeffContainer& coeffs_,
                                    Dune::template index_constant<ID> = Dune::template index_constant<std::size_t(0)>{})
        : geo_{geo},basis_{p_basis}, coeffs{coeffs_}, coeffsAsMat{Ikarus::viewAsEigenMatrixFixedDyn(coeffs)} {}

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

    using Traits = LocalFunctionTraits<NormalLocalFunction<Geometry,DuneBasis, CoeffContainer, ID>>;
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
      using other = NormalLocalFunction<
          Geometry,DuneBasis, typename Std::Rebind<CoeffContainer, typename Manifold::template Rebind<OtherType>::other>::other,
          ID>;
    };

    const Ikarus::LocalBasis<DuneBasis>& basis() const { return basis_; }

  private:
    template <typename DomainTypeOrIntegrationPointIndex, typename... TransformArgs>
    FunctionReturnType evaluateFunctionImpl(const DomainTypeOrIntegrationPointIndex& ipIndexOrPosition,
                                            [[maybe_unused]] const TransformWith<TransformArgs...>& transArgs) const {
      const auto& dNraw = evaluateDerivativeWithIPorCoord(ipIndexOrPosition, basis_);
      maytransformDerivatives(dNraw, dNTransformed, transArgs);
      auto a_1a_2 =  (toEigenMatrix(geo_.jacobianTransposed(ipIndexOrPosition))+ coeffsAsMat * dNTransformed.template cast<ctype>()).eval();
      //     [a_1 a_2]  = [A_1 A_2]                     +                   [u_1 u_2]

      return a_1a_2.col(0).cross(a_1a_2.col(1));
    }

    template <typename DomainTypeOrIntegrationPointIndex, typename... TransformArgs>
    Jacobian evaluateDerivativeWRTSpaceAllImpl(const DomainTypeOrIntegrationPointIndex& ipIndexOrPosition,
                                               const TransformWith<TransformArgs...>& transArgs) const {
      // get Second Derivatives of Ansatz Function
      const auto& dNraw = evaluateDerivativeWithIPorCoord(ipIndexOrPosition, basis_);
      maytransformDerivatives(dNraw, dNTransformed, transArgs);
      auto a_1a_2 =  (toEigenMatrix(geo_.jacobianTransposed(ipIndexOrPosition))+ coeffsAsMat * dNTransformed.template cast<ctype>()).eval();
      auto a1 =  a_1a_2.col(0).eval();
      auto a2 =  a_1a_2.col(1).eval();
     auto A_11A_22_A_12 =  toEigenMatrix(geo_.secondDerivativeOfPosition(ipIndexOrPosition));
     auto A_11 =  A_11A_22_A_12.col(0).eval();
     auto A_22 =  A_11A_22_A_12.col(1).eval();
     auto A_12 =  A_11A_22_A_12.col(2).eval();
     const auto& ddNraw = evaluateSecondDerivativeWithIPorCoord(ipIndexOrPosition, basis_);
     static_assert(transArgs.isEmpty,"Transform second derivatives not implemented.");

     auto u_11u_22_u12 = coeffsAsMat * ddNraw.template cast<ctype>();
     auto u_11 =  u_11u_22_u12.col(0).eval();
     auto u_22 = u_11u_22_u12.col(1).eval();
     auto u_12 =   u_11u_22_u12.col(2).eval();

     auto a_11 =  (A_11+u_11).eval();
     auto a_22 =  (A_22+u_22).eval();
     auto a_12 =  (A_12+u_12).eval();
     Jacobian res;
     res.col(0) = a_11.cross(a2) + a1.cross(a_12);
     res.col(1) = a_12.cross(a2) + a1.cross(a_22);

     return res;
    }

    template <typename DomainTypeOrIntegrationPointIndex, typename... TransformArgs>
    JacobianColType evaluateDerivativeWRTSpaceSingleImpl(const DomainTypeOrIntegrationPointIndex& ipIndexOrPosition,
                                                         int spaceIndex,
                                                         const TransformWith<TransformArgs...>& transArgs) const {
      const auto& dNraw = evaluateDerivativeWithIPorCoord(ipIndexOrPosition, basis_);
      maytransformDerivatives(dNraw, dNTransformed, transArgs);
      auto a_1a_2 =  (toEigenMatrix(geo_.jacobianTransposed(ipIndexOrPosition))+ coeffsAsMat * dNTransformed.template cast<ctype>()).eval();

      auto A_11A_22_A_12 =  toEigenMatrix(geo_.secondDerivativeOfPosition(ipIndexOrPosition));
      const auto& ddNraw = evaluateSecondDerivativeWithIPorCoord(ipIndexOrPosition, basis_);
      static_assert(transArgs.isEmpty,"Transform second derivatives not implemented.");

      auto u_11u_22_u12 = coeffsAsMat * ddNraw.template cast<ctype>();

      auto a_11a_22_a_12 = A_11A_22_A_12+u_11u_22_u12;
      const int otherIndex = spaceIndex==0? 1: 0;
      Jacobian res;
      //  a_11a_22_a_12[0].cross(a2) + a1.cross(a_11a_22_a_12[2]);
      // -a_11a_22_a_12[1].cross(a1) - a2.cross(a_11a_22_a_12[2]) ;
      return (a_11a_22_a_12[spaceIndex].cross(a_1a_2.col(otherIndex)) + a_1a_2.col(spaceIndex).cross(a_11a_22_a_12[2]))* (spaceIndex==0? 1: -1);

    }

    template <typename DomainTypeOrIntegrationPointIndex, typename... TransformArgs>
    CoeffDerivMatrix evaluateDerivativeWRTCoeffsImpl(const DomainTypeOrIntegrationPointIndex& ipIndexOrPosition,
                                                     int coeffsIndex,
                                                     const TransformWith<TransformArgs...>& transArgs) const {
      const auto& dNraw = evaluateDerivativeWithIPorCoord(ipIndexOrPosition, basis_);
      maytransformDerivatives(dNraw, dNTransformed, transArgs);
      auto a_1a_2 =  (toEigenMatrix(geo_.jacobianTransposed(ipIndexOrPosition))+ coeffsAsMat * dNTransformed.template cast<ctype>()).eval();
      //     [a_1 a_2]  = [A_1 A_2]                     +                   [u_1 u_2]
      Eigen::DiagonalMatrix<ctype,3> a1dC = Eigen::Matrix<ctype,3,3>::Identity()* dNTransformed(coeffsIndex,0);
      Eigen::DiagonalMatrix<ctype,3> a2dC = Eigen::Matrix<ctype,3,3>::Identity()* dNTransformed(coeffsIndex,1);
      Eigen::Matrix<ctype,3,3> a1Hat = toSkewMatrix(a_1a_2.col(0));
      Eigen::Matrix<ctype,3,3> a2Hat = toSkewMatrix(a_1a_2.col(1));

      return a1Hat*a2dC-a2Hat*a1dC;
    }

    template <typename DomainTypeOrIntegrationPointIndex, typename... TransformArgs>
    std::array<CoeffDerivMatrix, gridDim> evaluateDerivativeWRTCoeffsANDSpatialImpl(
        const DomainTypeOrIntegrationPointIndex& ipIndexOrPosition, int coeffsIndex,
        const TransformWith<TransformArgs...>& transArgs) const {
      const auto& dNraw = evaluateDerivativeWithIPorCoord(ipIndexOrPosition, basis_);
      const auto& ddNraw = evaluateSecondDerivativeWithIPorCoord(ipIndexOrPosition, basis_);
      static_assert(transArgs.isEmpty,"Transform second derivatives not implemented.");
      maytransformDerivatives(dNraw, dNTransformed, transArgs);
      auto a_1a_2 =  (toEigenMatrix(geo_.jacobianTransposed(ipIndexOrPosition))+ coeffsAsMat * dNTransformed.template cast<ctype>()).eval();
      //     [a_1 a_2]  = [A_1 A_2]                     +                   [u_1 u_2]
      Eigen::DiagonalMatrix<ctype,3> a1_C = Eigen::Matrix<ctype,3,3>::Identity()* dNTransformed(coeffsIndex,0);
      Eigen::DiagonalMatrix<ctype,3> a2_C = Eigen::Matrix<ctype,3,3>::Identity()* dNTransformed(coeffsIndex,1);

      // get Second Derivatives of Ansatz Function
      auto a1 =  a_1a_2.col(0).eval();
      auto a2 =  a_1a_2.col(1).eval();
      auto A_11A_22_A_12 =  toEigenMatrix(geo_.secondDerivativeOfPosition(ipIndexOrPosition));
      auto A_11 =  A_11A_22_A_12.col(0).eval();
      auto A_22 =  A_11A_22_A_12.col(1).eval();
      auto A_12 =  A_11A_22_A_12.col(2).eval();

      auto u_11u_22_u12 = coeffsAsMat * ddNraw.template cast<ctype>();
      auto u_11 =  u_11u_22_u12.col(0).eval();
      auto u_22 = u_11u_22_u12.col(1).eval();
      auto u_12 =   u_11u_22_u12.col(2).eval();

      auto a_11 =  (A_11+u_11).eval();
      auto a_22 =  (A_22+u_22).eval();
      auto a_12 =  (A_12+u_12).eval();
      Eigen::DiagonalMatrix<ctype,3> a_11_C = Eigen::Matrix<ctype,3,3>::Identity()* ddNraw(coeffsIndex,0);
      Eigen::DiagonalMatrix<ctype,3> a_22_C = Eigen::Matrix<ctype,3,3>::Identity()* ddNraw(coeffsIndex,1);
      Eigen::DiagonalMatrix<ctype,3> a_12_C = Eigen::Matrix<ctype,3,3>::Identity()* ddNraw(coeffsIndex,2);
      std::array<CoeffDerivMatrix, 2> a3_1_C_a3_2_C;
      a3_1_C_a3_2_C[0] = -cross(a2,a_11_C) + cross(a_11,a2_C) -cross(a_12,a1_C) + cross(a1,a_12_C);
      a3_1_C_a3_2_C[1] = -cross(a2,a_12_C) + cross(a_12,a2_C) -cross(a_22,a1_C) + cross(a1,a_22_C);

      return a3_1_C_a3_2_C;
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
    const Geometry& geo_;
    CoeffContainer coeffs;
    const decltype(Ikarus::viewAsEigenMatrixFixedDyn(coeffs)) coeffsAsMat;
  };

  template <typename Geometry,typename DuneBasis, typename CoeffContainer, std::size_t ID>
  struct LocalFunctionTraits<NormalLocalFunction<Geometry,DuneBasis, CoeffContainer, ID>> {
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
