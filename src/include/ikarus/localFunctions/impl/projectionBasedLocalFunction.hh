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

#include <Eigen/Core>
#include <Eigen/Dense>

#include <ikarus/localBasis/localBasis.hh>
#include <ikarus/localFunctions/localFunctionHelper.hh>
#include <ikarus/localFunctions/localFunctionInterface.hh>
#include <ikarus/utils/linearAlgebraHelper.hh>

namespace Ikarus {

  template <typename DuneBasis, typename CoeffContainer, std::size_t ID = 0>
  class ProjectionBasedLocalFunction
      : public LocalFunctionInterface<ProjectionBasedLocalFunction<DuneBasis, CoeffContainer, ID>>,
        public ClonableLocalFunction<ProjectionBasedLocalFunction<DuneBasis, CoeffContainer, ID>> {
    using Interface = LocalFunctionInterface<ProjectionBasedLocalFunction<DuneBasis, CoeffContainer, ID>>;

    template <size_t ID_ = 0>
    static constexpr int orderID = ID_ == ID ? nonLinear : constant;

  public:
    friend Interface;
    friend ClonableLocalFunction<ProjectionBasedLocalFunction>;
    constexpr ProjectionBasedLocalFunction(
        const Ikarus::LocalBasis<DuneBasis>& p_basis, const CoeffContainer& coeffs_,
        Dune::template index_constant<ID> = Dune::template index_constant<std::size_t(0)>{})
        : basis_{p_basis}, coeffs{coeffs_}, coeffsAsMat{Ikarus::viewAsEigenMatrixFixedDyn(coeffs)} {}

    using Traits = LocalFunctionTraits<ProjectionBasedLocalFunction>;

    static constexpr bool isLeaf = true;
    using Ids                    = Dune::index_constant<ID>;

    template <typename LocalFunctionEvaluationArgs_, typename LocalFunctionImpl_>
    friend auto evaluateDerivativeImpl(const LocalFunctionInterface<LocalFunctionImpl_>& f,
                                       const LocalFunctionEvaluationArgs_& localFunctionArgs);

    template <typename LocalFunctionEvaluationArgs_, typename LocalFunctionImpl_>
    friend auto evaluateFunctionImpl(const LocalFunctionInterface<LocalFunctionImpl_>& f,
                                     const LocalFunctionEvaluationArgs_& localFunctionArgs);

    /** \brief Type used for coordinates */
    using ctype = typename Traits::ctype;
    //    /** \brief Dimension of the coeffs */
    static constexpr int valueSize = Traits::valueSize;
    //    /** \brief Dimension of the coeffs */
    static constexpr int correctionSize = Traits::correctionSize;
    /** \brief Dimension of the grid */
    static constexpr int gridDim = Traits::gridDim;
    /** \brief Type for coordinate vector in world space */
    using FunctionReturnType = typename Traits::FunctionReturnType;
    /** \brief Type for the directional derivatives */
    using AlongType = typename Traits::AlongType;
    /** \brief The manifold where the function values lives in */
    using Manifold = typename Traits::Manifold;
    /** \brief Type for the Jacobian matrix */
    using Jacobian = typename Traits::Jacobian;
    /** \brief Type for a column of the Jacobian matrix */
    using JacobianColType = typename Traits::JacobianColType;
    /** \brief Type for the derivatives wrT the coefficients in the local tangent space base*/
    using CoeffDerivMatrix = typename Traits::CoeffDerivMatrix;
    /** \brief Type for the derivatives wrT the coefficients in the embedding space*/
    using CoeffDerivEukMatrix = typename Traits::CoeffDerivEukMatrix;
    /** \brief Type for the derivatives wrT the coefficients in the embedding space on the left and in the tangent space
     * basis on the right*/
    using CoeffDerivEukRieMatrix = Eigen::Matrix<ctype, valueSize, correctionSize>;
    /** \brief Type for ansatz function values */
    using AnsatzFunctionType = typename Traits::AnsatzFunctionType;
    /** \brief Type for the Jacobian of the ansatz function values */
    using AnsatzFunctionJacobian = typename Traits::AnsatzFunctionJacobian;

    const auto& coefficientsRef() const { return coeffs; }
    auto& coefficientsRef() { return coeffs; }

    const Ikarus::LocalBasis<DuneBasis>& basis() const { return basis_; }

    template <typename OtherType>
    struct Rebind {
      using other = ProjectionBasedLocalFunction<
          DuneBasis, typename Std::Rebind<CoeffContainer, typename Manifold::template Rebind<OtherType>::other>::other,
          ID>;
    };

  private:
    static auto tryToCallDerivativeOfProjectionWRTposition(const FunctionReturnType& valE) {
      if constexpr (requires { Manifold::derivativeOfProjectionWRTposition(valE); })
        return Manifold::derivativeOfProjectionWRTposition(valE);
      else
        static_assert(
            requires { Manifold::derivativeOfProjectionWRTposition(valE); },
            " Your passed manifold does not implement derivativeOfProjectionWRTposition.");
    }

    static auto tryToCallSecondDerivativeOfProjectionWRTposition(const FunctionReturnType& valE,
                                                                 const AlongType& along) {
      if constexpr (requires { Manifold::secondDerivativeOfProjectionWRTposition(valE, along); })
        return Manifold::secondDerivativeOfProjectionWRTposition(valE, along);
      else
        static_assert(
            requires { Manifold::secondDerivativeOfProjectionWRTposition(valE, along); },
            " Your passed manifold does not implement derivativeOfProjectionWRTposition.");
    }

    static auto tryToCallThirdDerivativeOfProjectionWRTposition(const FunctionReturnType& valE, const AlongType& along,
                                                                const Eigen::Ref<const AlongType>& along2) {
      if constexpr (requires { Manifold::thirdDerivativeOfProjectionWRTposition(valE, along, along2); })
        return Manifold::thirdDerivativeOfProjectionWRTposition(valE, along, along2);
      else
        static_assert(
            requires { Manifold::thirdDerivativeOfProjectionWRTposition(valE, along, along2); },
            " Your passed manifold does not implement derivativeOfProjectionWRTposition.");
    }

    template <typename DomainTypeOrIntegrationPointIndex, typename... TransformArgs>
    Jacobian evaluateDerivativeWRTSpaceAllImpl(const DomainTypeOrIntegrationPointIndex& ipIndexOrPosition,
                                               const TransformWith<TransformArgs...>& transArgs) const {
      const auto& [N, dNraw] = evaluateFunctionAndDerivativeWithIPorCoord(ipIndexOrPosition, basis_);
      maytransformDerivatives(dNraw, dNTransformed, transArgs);
      Jacobian J              = evaluateEmbeddingJacobianImpl(dNTransformed);
      FunctionReturnType valE = evaluateEmbeddingFunctionImpl(N);
      return tryToCallDerivativeOfProjectionWRTposition(valE) * J;
    }

    template <typename DomainTypeOrIntegrationPointIndex, typename... TransformArgs>
    JacobianColType evaluateDerivativeWRTSpaceSingleImpl(const DomainTypeOrIntegrationPointIndex& ipIndexOrPosition,
                                                         int spaceIndex,
                                                         const TransformWith<TransformArgs...>& transArgs) const {
      const auto& [N, dNraw] = evaluateFunctionAndDerivativeWithIPorCoord(ipIndexOrPosition, basis_);
      maytransformDerivatives(dNraw, dNTransformed, transArgs);
      JacobianColType Jcol    = evaluateEmbeddingJacobianColImpl(dNTransformed, spaceIndex);
      FunctionReturnType valE = evaluateEmbeddingFunctionImpl(N);
      return tryToCallDerivativeOfProjectionWRTposition(valE) * Jcol;
    }

    template <typename DomainTypeOrIntegrationPointIndex, typename... TransformArgs>
    CoeffDerivEukRieMatrix evaluateDerivativeWRTCoeffsImpl(const DomainTypeOrIntegrationPointIndex& ipIndexOrPosition,
                                                           int coeffsIndex,
                                                           const TransformWith<TransformArgs...>&) const {
      const auto& N = evaluateFunctionWithIPorCoord(ipIndexOrPosition, basis_);
      return (evaluateDerivativeWRTCoeffsEukImpl(N, coeffsIndex) * coeffs[coeffsIndex].orthonormalFrame()).eval();
    }

    CoeffDerivEukMatrix evaluateDerivativeWRTCoeffsEukImpl(const AnsatzFunctionType& N, int coeffsIndex) const {
      FunctionReturnType valE = evaluateEmbeddingFunctionImpl(N);
      return (tryToCallDerivativeOfProjectionWRTposition(valE) * N[coeffsIndex]).eval();
    }

    template <typename DomainTypeOrIntegrationPointIndex, typename... AlongArgs, typename... TransformArgs>
    CoeffDerivMatrix evaluateSecondDerivativeWRTCoeffsImpl(const DomainTypeOrIntegrationPointIndex& ipIndexOrPosition,
                                                           const std::array<size_t, 2>& coeffsIndex,
                                                           const Along<AlongArgs...>& alongArgs,
                                                           const TransformWith<TransformArgs...>&) const {
      const auto& N                 = evaluateFunctionWithIPorCoord(ipIndexOrPosition, basis_);
      const FunctionReturnType valE = evaluateEmbeddingFunctionImpl(N);

      CoeffDerivEukMatrix ddt = tryToCallSecondDerivativeOfProjectionWRTposition(valE, std::get<0>(alongArgs.args))
                                * N[coeffsIndex[0]] * N[coeffsIndex[1]];

      if (coeffsIndex[0] == coeffsIndex[1]) {  // Riemannian Hessian Weingarten map correction
        const CoeffDerivEukMatrix dt = evaluateDerivativeWRTCoeffsEukImpl(N, coeffsIndex[0]);
        ddt -= (coeffs[coeffsIndex[0]].getValue().dot(dt * std::get<0>(alongArgs.args)))
               * CoeffDerivEukMatrix::Identity();
      }

      return coeffs[coeffsIndex[0]].orthonormalFrame().transpose() * ddt * coeffs[coeffsIndex[1]].orthonormalFrame();
    }

    template <typename DomainTypeOrIntegrationPointIndex, typename... TransformArgs>
    std::array<CoeffDerivEukRieMatrix, gridDim> evaluateDerivativeWRTCoeffsANDSpatialImpl(
        const DomainTypeOrIntegrationPointIndex& ipIndexOrPosition, int coeffsIndex,
        const TransformWith<TransformArgs...>& transArgs) const {
      std::array<CoeffDerivEukMatrix, gridDim> WarrayEuk
          = evaluateDerivativeWRTCoeffsANDSpatialEukImpl(ipIndexOrPosition, coeffsIndex, transArgs);
      std::array<CoeffDerivEukRieMatrix, gridDim> WarrayRie;
      const auto BLA = coeffs[coeffsIndex].orthonormalFrame().eval();
      for (int dir = 0; dir < gridDim; ++dir)
        WarrayRie[dir] = WarrayEuk[dir] * BLA;

      return WarrayRie;
    }

    template <typename DomainTypeOrIntegrationPointIndex, typename... TransformArgs>
    std::array<CoeffDerivEukMatrix, gridDim> evaluateDerivativeWRTCoeffsANDSpatialEukImpl(
        const DomainTypeOrIntegrationPointIndex& ipIndexOrPosition, int coeffsIndex,
        const TransformWith<TransformArgs...>& transArgs) const {
      const auto& [N, dNraw] = evaluateFunctionAndDerivativeWithIPorCoord(ipIndexOrPosition, basis_);
      maytransformDerivatives(dNraw, dNTransformed, transArgs);
      const FunctionReturnType valE = evaluateEmbeddingFunctionImpl(N);
      const Jacobian J              = evaluateEmbeddingJacobianImpl(dNTransformed);
      const CoeffDerivEukMatrix Pm  = tryToCallDerivativeOfProjectionWRTposition(valE);
      std::array<CoeffDerivEukMatrix, gridDim> Warray;
      for (int dir = 0; dir < gridDim; ++dir) {
        const auto Qi = tryToCallSecondDerivativeOfProjectionWRTposition(valE, J.col(dir));
        Warray[dir]   = Qi * N[coeffsIndex] + Pm * dNTransformed(coeffsIndex, dir);
      }

      return Warray;
    }

    template <typename DomainTypeOrIntegrationPointIndex, typename... TransformArgs>
    CoeffDerivEukRieMatrix evaluateDerivativeWRTCoeffsANDSpatialSingleImpl(
        const DomainTypeOrIntegrationPointIndex& ipIndexOrPosition, int coeffsIndex, int spatialIndex,
        const TransformWith<TransformArgs...>& transArgs) const {
      const CoeffDerivEukMatrix WEuk
          = evaluateDerivativeWRTCoeffsANDSpatialSingleEukImpl(ipIndexOrPosition, coeffsIndex, spatialIndex, transArgs);
      return WEuk * coeffs[coeffsIndex].orthonormalFrame();
    }

    template <typename DomainTypeOrIntegrationPointIndex, typename... TransformArgs>
    CoeffDerivEukMatrix evaluateDerivativeWRTCoeffsANDSpatialSingleEukImpl(
        const DomainTypeOrIntegrationPointIndex& ipIndexOrPosition, int coeffsIndex, int spatialIndex,
        const TransformWith<TransformArgs...>& transArgs) const {
      const auto& [N, dNraw] = evaluateFunctionAndDerivativeWithIPorCoord(ipIndexOrPosition, basis_);
      maytransformDerivatives(dNraw, dNTransformed, transArgs);
      const FunctionReturnType valE = evaluateEmbeddingFunctionImpl(N);
      const JacobianColType Jcol    = evaluateEmbeddingJacobianColImpl(dNTransformed, spatialIndex);
      const CoeffDerivEukMatrix Pm  = tryToCallDerivativeOfProjectionWRTposition(valE);
      CoeffDerivEukMatrix W;
      const auto Qi = tryToCallSecondDerivativeOfProjectionWRTposition(valE, Jcol);
      W             = Qi * N[coeffsIndex] + Pm * dNTransformed(coeffsIndex, spatialIndex);
      return W;
    }

    template <typename DomainTypeOrIntegrationPointIndex, typename... AlongArgs, typename... TransformArgs>
    auto evaluateThirdDerivativeWRTCoeffsTwoTimesAndSpatialImpl(
        const DomainTypeOrIntegrationPointIndex& ipIndexOrPosition, const std::array<size_t, 2>& coeffsIndex,
        const Along<AlongArgs...>& alongArgs, const TransformWith<TransformArgs...>& transArgs) const {
      const auto& [N, dNraw] = evaluateFunctionAndDerivativeWithIPorCoord(ipIndexOrPosition, basis_);
      maytransformDerivatives(dNraw, dNTransformed, transArgs);
      const FunctionReturnType valE = evaluateEmbeddingFunctionImpl(N);
      const Jacobian J              = evaluateEmbeddingJacobianImpl(dNTransformed);
      const auto& along             = std::get<0>(alongArgs.args);
      CoeffDerivEukMatrix ChiArrayEuk;
      ChiArrayEuk.setZero();
      for (int i = 0; i < gridDim; ++i) {
        const auto chi              = tryToCallThirdDerivativeOfProjectionWRTposition(valE, along.col(i), J.col(i));
        const CoeffDerivEukMatrix S = tryToCallSecondDerivativeOfProjectionWRTposition(valE, along.col(i));
        const auto& NI              = N[coeffsIndex[0]];
        const auto& NJ              = N[coeffsIndex[1]];
        const auto& dNIdi           = dNTransformed(coeffsIndex[0], i);
        const auto& dNJdi           = dNTransformed(coeffsIndex[1], i);
        ChiArrayEuk += chi * NI * NJ + S * (dNIdi * NJ + dNJdi * NI);
      }
      if (coeffsIndex[0] == coeffsIndex[1]) {  // Riemannian Hessian Weingarten map correction
        const std::array<CoeffDerivEukMatrix, gridDim> Warray
            = evaluateDerivativeWRTCoeffsANDSpatialEukImpl(ipIndexOrPosition, coeffsIndex[0], transArgs);
        for (int i = 0; i < gridDim; ++i) {
          ChiArrayEuk
              -= coeffs[coeffsIndex[0]].getValue().dot(Warray[i] * along.col(i)) * CoeffDerivEukMatrix::Identity();
        }
      }

      CoeffDerivMatrix ChiArrayRie;
      const auto BLA0T = coeffs[coeffsIndex[0]].orthonormalFrame().transpose().eval();
      const auto BLA1  = coeffs[coeffsIndex[1]].orthonormalFrame();
      ChiArrayRie      = BLA0T * ChiArrayEuk * BLA1;
      return ChiArrayRie;
    }

    template <typename DomainTypeOrIntegrationPointIndex, typename... AlongArgs, typename... TransformArgs>
    CoeffDerivMatrix evaluateThirdDerivativeWRTCoeffsTwoTimesAndSpatialSingleImpl(
        const DomainTypeOrIntegrationPointIndex& ipIndexOrPosition, const std::array<size_t, 2>& coeffsIndex,
        const int spatialIndex, const Along<AlongArgs...>& alongArgs,
        const TransformWith<TransformArgs...>& transArgs) const {
      const auto& [N, dNraw] = evaluateFunctionAndDerivativeWithIPorCoord(ipIndexOrPosition, basis_);
      maytransformDerivatives(dNraw, dNTransformed, transArgs);
      const FunctionReturnType valE = evaluateEmbeddingFunctionImpl(N);
      const Jacobian J              = evaluateEmbeddingJacobianImpl(dNTransformed);
      const auto& along             = std::get<0>(alongArgs.args);
      const CoeffDerivEukMatrix S   = tryToCallSecondDerivativeOfProjectionWRTposition(valE, along);
      const auto chi                = tryToCallThirdDerivativeOfProjectionWRTposition(valE, along, J.col(spatialIndex));
      const auto& NI                = N[coeffsIndex[0]];
      const auto& NJ                = N[coeffsIndex[1]];
      const auto& dNIdi             = dNTransformed(coeffsIndex[0], spatialIndex);
      const auto& dNJdi             = dNTransformed(coeffsIndex[1], spatialIndex);
      CoeffDerivEukMatrix Chi       = chi * NI * NJ + S * (dNIdi * NJ + dNJdi * NI);

      if (coeffsIndex[0] == coeffsIndex[1]) {  // Riemannian Hessian Weingarten map correction
        const CoeffDerivEukMatrix W = evaluateDerivativeWRTCoeffsANDSpatialSingleEukImpl(
            ipIndexOrPosition, coeffsIndex[0], spatialIndex, transArgs);
        Chi -= coeffs[coeffsIndex[0]].getValue().dot(W * along) * CoeffDerivEukMatrix::Identity();
      }
      return (coeffs[coeffsIndex[0]].orthonormalFrame().transpose() * Chi * coeffs[coeffsIndex[1]].orthonormalFrame())
          .eval();
    }

    template <typename DomainTypeOrIntegrationPointIndex, typename... TransformArgs>
    FunctionReturnType evaluateFunctionImpl(const DomainTypeOrIntegrationPointIndex& ipIndexOrPosition,
                                            [[maybe_unused]] const TransformWith<TransformArgs...>&) const {
      const auto& N = evaluateFunctionWithIPorCoord(ipIndexOrPosition, basis_);
      return Manifold(evaluateEmbeddingFunctionImpl(N)).getValue();
    }

    JacobianColType evaluateEmbeddingJacobianColImpl(const AnsatzFunctionJacobian& dN, int spaceIndex) const {
      JacobianColType Jcol = coeffsAsMat * dN.col(spaceIndex);
      return Jcol;
    }

    Jacobian evaluateEmbeddingJacobianImpl(const AnsatzFunctionJacobian& dN) const {
      Jacobian J
          = coeffsAsMat
            * dN.template cast<ctype>();  // The cast here is only necessary since the autodiff types are not working
      // otherwise, see Issue https://github.com/autodiff/autodiff/issues/73
      return J;
    }

    FunctionReturnType evaluateEmbeddingFunctionImpl(const AnsatzFunctionType& N) const { return coeffsAsMat * N; }

    mutable AnsatzFunctionJacobian dNTransformed;
    Ikarus::LocalBasis<DuneBasis> basis_;
    CoeffContainer coeffs;
    const decltype(Ikarus::viewAsEigenMatrixFixedDyn(coeffs)) coeffsAsMat;
  };

  template <typename DuneBasis, typename CoeffContainer, std::size_t ID>
  struct LocalFunctionTraits<ProjectionBasedLocalFunction<DuneBasis, CoeffContainer, ID>> {
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
    /** \brief Type for the derivatives wrt. the coefficients */
    using CoeffDerivMatrix = Eigen::Matrix<ctype, correctionSize, correctionSize>;
    /** \brief Type for the derivatives wrt. the coefficients */
    using CoeffDerivEukMatrix = Eigen::Matrix<ctype, valueSize, valueSize>;
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
