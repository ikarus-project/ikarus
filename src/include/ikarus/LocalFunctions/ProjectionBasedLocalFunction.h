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
  class ProjectionBasedLocalFunction
      : public LocalFunctionInterface<ProjectionBasedLocalFunction<DuneBasis, CoeffContainer>> {
    using Base = LocalFunctionInterface<ProjectionBasedLocalFunction<DuneBasis, CoeffContainer>>;

  public:
    friend Base;
    ProjectionBasedLocalFunction(const Ikarus::LocalBasis<DuneBasis>& basis_, const CoeffContainer& coeffs_)
        : basis{basis_}, coeffs{coeffs_}, coeffsAsMat{Ikarus::LinearAlgebra::viewAsEigenMatrixFixedDyn(coeffs)} {}

    using Traits = LocalFunctionTraits<ProjectionBasedLocalFunction>;

    /** \brief Type used for coordinates */
    using ctype = typename Traits::ctype;
    //    /** \brief Dimension of the coeffs */
    static constexpr int valueSize = Traits::valueSize;
    /** \brief Dimension of the grid */
    static constexpr int gridDim = Traits::gridDim;
    /** \brief Type for coordinate vector in world space */
    using FunctionReturnType = typename Traits::FunctionReturnType;
    /** \brief Type for the directional derivatives */
    using AlongType = Eigen::Vector<ctype, valueSize>;
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
    /** \brief Matrix to transform the ansatz function Jacobian to world coordinates*/
    using TransformMatrix = typename Traits::TransformMatrix;

    auto& coefficientsRef() { return coeffs; }

  private:
    static auto tryToCallDerivativeOfProjectionWRTposition(const GlobalE& valE) {
      if constexpr (requires { FunctionReturnType::derivativeOfProjectionWRTposition(valE); })
        return FunctionReturnType::derivativeOfProjectionWRTposition(valE);
      else
        static_assert(
            requires { FunctionReturnType::derivativeOfProjectionWRTposition(valE); },
            " Your passed manifold does not implement derivativeOfProjectionWRTposition.");
    }

    static auto tryToCallSecondDerivativeOfProjectionWRTposition(const GlobalE& valE, const AlongType& along) {
      if constexpr (requires { FunctionReturnType::secondDerivativeOfProjectionWRTposition(valE, along); })
        return FunctionReturnType::secondDerivativeOfProjectionWRTposition(valE, along);
      else
        static_assert(
            requires { FunctionReturnType::secondDerivativeOfProjectionWRTposition(valE, along); },
            " Your passed manifold does not implement derivativeOfProjectionWRTposition.");
    }

    static auto tryToCallThirdDerivativeOfProjectionWRTposition(const GlobalE& valE, const AlongType& along,
                                                                const Eigen::Ref<const AlongType>& along2) {
      if constexpr (requires { FunctionReturnType::thirdDerivativeOfProjectionWRTposition(valE, along, along2); })
        return FunctionReturnType::thirdDerivativeOfProjectionWRTposition(valE, along, along2);
      else
        static_assert(
            requires { FunctionReturnType::thirdDerivativeOfProjectionWRTposition(valE, along, along2); },
            " Your passed manifold does not implement derivativeOfProjectionWRTposition.");
    }

    Jacobian evaluateDerivativeWRTSpaceAllImpl(const AnsatzFunctionType& N, const AnsatzFunctionJacobian& dN) const {
      Jacobian J   = evaluateEmbeddingJacobianImpl(dN);
      GlobalE valE = evaluateEmbeddingFunctionImpl(N);
      return tryToCallDerivativeOfProjectionWRTposition(valE) * J;
    }

    JacobianColType evaluateDerivativeWRTSpaceSingleImpl(const AnsatzFunctionType& N, const AnsatzFunctionJacobian& dN,
                                                         int spaceIndex) const {
      JacobianColType Jcol = evaluateEmbeddingJacobianColImpl(dN, spaceIndex);
      GlobalE valE         = evaluateEmbeddingFunctionImpl(N);
      return tryToCallDerivativeOfProjectionWRTposition(valE) * Jcol;
    }

    auto evaluateDerivativeWRTCoeffsImpl(const AnsatzFunctionType& N, [[maybe_unused]] const AnsatzFunctionJacobian&,
                                         int coeffsIndex) const {
      GlobalE valE = evaluateEmbeddingFunctionImpl(N);
      return (tryToCallDerivativeOfProjectionWRTposition(valE) * N[coeffsIndex]).eval();
    }

    CoeffDerivMatrix evaluateSecondDerivativeWRTCoeffs(const AnsatzFunctionType& N,
                                                       [[maybe_unused]] const AnsatzFunctionJacobian&,
                                                       const AlongType& along,
                                                       const std::array<size_t, gridDim>& coeffsIndex) const {
      const GlobalE valE = evaluateEmbeddingFunctionImpl(N);
      return tryToCallSecondDerivativeOfProjectionWRTposition(valE, along) * N[coeffsIndex[0]] * N[coeffsIndex[1]];
    }

    auto evaluateDerivativeWRTCoeffsANDSpatialImpl(const AnsatzFunctionType& N,
                                                   [[maybe_unused]] const AnsatzFunctionJacobian& dN,
                                                   int coeffsIndex) const {
      const GlobalE valE        = evaluateEmbeddingFunctionImpl(N);
      const Jacobian J          = evaluateEmbeddingJacobianImpl(dN);
      const CoeffDerivMatrix Pm = tryToCallDerivativeOfProjectionWRTposition(valE);
      std::array<CoeffDerivMatrix, gridDim> Warray;
      for (int dir = 0; dir < gridDim; ++dir) {
        const auto Qi = tryToCallSecondDerivativeOfProjectionWRTposition(valE, J.col(dir));
        Warray[dir]   = Qi * N[coeffsIndex] + Pm * dN(coeffsIndex, dir);
      }

      return Warray;
    }

    auto evaluateDerivativeWRTCoeffsANDSpatialSingleImpl(const AnsatzFunctionType& N,
                                                         [[maybe_unused]] const AnsatzFunctionJacobian& dN,
                                                         int coeffsIndex, const int spatialIndex) const {
      const GlobalE valE         = evaluateEmbeddingFunctionImpl(N);
      const JacobianColType Jcol = evaluateEmbeddingJacobianColImpl(dN, spatialIndex);
      const CoeffDerivMatrix Pm  = tryToCallDerivativeOfProjectionWRTposition(valE);
      CoeffDerivMatrix W;
      const auto Qi = tryToCallSecondDerivativeOfProjectionWRTposition(valE, Jcol);
      W             = Qi * N[coeffsIndex] + Pm * dN(coeffsIndex, spatialIndex);

      return W;
    }

    auto evaluateThirdDerivativeWRTCoeffsTwoTimesAndSpatialImpl(const AnsatzFunctionType& N,
                                                                [[maybe_unused]] const AnsatzFunctionJacobian& dN,
                                                                const AlongType& along,
                                                                const std::array<size_t, gridDim>& coeffsIndex) const {
      const GlobalE valE       = evaluateEmbeddingFunctionImpl(N);
      const Jacobian J         = evaluateEmbeddingJacobianImpl(dN);
      const CoeffDerivMatrix S = tryToCallSecondDerivativeOfProjectionWRTposition(valE, along);
      std::array<CoeffDerivMatrix, gridDim> ChiArray;
      for (int i = 0; i < gridDim; ++i) {
        const auto chi    = tryToCallThirdDerivativeOfProjectionWRTposition(valE, along, J.col(i));
        const auto& NI    = N[coeffsIndex[0]];
        const auto& NJ    = N[coeffsIndex[1]];
        const auto& dNIdi = dN(coeffsIndex[0], i);
        const auto& dNJdi = dN(coeffsIndex[1], i);
        ChiArray[i]       = chi * NI * NJ + S * (dNIdi * NJ + dNJdi * NI);
      }

      return ChiArray;
    }

    auto evaluateThirdDerivativeWRTCoeffsTwoTimesAndSpatialSingleImpl(const AnsatzFunctionType& N,
                                                                      [[maybe_unused]] const AnsatzFunctionJacobian& dN,
                                                                      const AlongType& along,
                                                                      const std::array<size_t, gridDim>& coeffsIndex,
                                                                      const int spatialIndex) const {
      const GlobalE valE       = evaluateEmbeddingFunctionImpl(N);
      const Jacobian J         = evaluateEmbeddingJacobianImpl(dN);
      const CoeffDerivMatrix S = tryToCallSecondDerivativeOfProjectionWRTposition(valE, along);
      CoeffDerivMatrix Chi;
      const auto chi    = tryToCallThirdDerivativeOfProjectionWRTposition(valE, along, J.col(spatialIndex));
      const auto& NI    = N[coeffsIndex[0]];
      const auto& NJ    = N[coeffsIndex[1]];
      const auto& dNIdi = dN(coeffsIndex[0], spatialIndex);
      const auto& dNJdi = dN(coeffsIndex[1], spatialIndex);
      Chi               = chi * NI * NJ + S * (dNIdi * NJ + dNJdi * NI);

      return Chi;
    }

    FunctionReturnType evaluateFunctionImpl(const AnsatzFunctionType& N) const {
      return FunctionReturnType(evaluateEmbeddingFunctionImpl(N));
    }

    JacobianColType evaluateEmbeddingJacobianColImpl(const AnsatzFunctionJacobian& dN, int spaceIndex) const {
      JacobianColType Jcol = coeffsAsMat * dN.col(spaceIndex);
      return Jcol;
    }

    Jacobian evaluateEmbeddingJacobianImpl(const AnsatzFunctionJacobian& dN) const {
      Jacobian J = coeffsAsMat * dN;
      return J;
    }

    GlobalE evaluateEmbeddingFunctionImpl(const Eigen::VectorXd& N) const { return coeffsAsMat * N; }

    const Ikarus::LocalBasis<DuneBasis>& basis;
    CoeffContainer coeffs;
    const decltype(Ikarus::LinearAlgebra::viewAsEigenMatrixFixedDyn(coeffs)) coeffsAsMat;
  };

  template <typename DuneBasis, typename CoeffContainer>
  struct LocalFunctionTraits<ProjectionBasedLocalFunction<DuneBasis, CoeffContainer>> {
      /** \brief Type used for coordinates */
      using ctype                               = typename CoeffContainer::value_type::ctype;
      /** \brief Dimension of the coeffs */
      static constexpr int valueSize = CoeffContainer::value_type::valueSize;
      /** \brief Dimension of the grid */
      static constexpr int gridDim = Ikarus::LocalBasis<DuneBasis>::gridDim;
      /** \brief Type for the return value */
      using FunctionReturnType = typename CoeffContainer::value_type;
      /** \brief Type for the Jacobian matrix */
      using Jacobian               = Eigen::Matrix<ctype, valueSize, gridDim>;
      /** \brief Type for the derivatives wrt. the coeffiecients */
      using CoeffDerivMatrix       = Eigen::Matrix<ctype, valueSize,valueSize>;
      /** \brief Type for the Jacobian of the ansatz function values */
      using AnsatzFunctionJacobian = typename Ikarus::LocalBasis<DuneBasis>::JacobianType;
      /** \brief Type for ansatz function values */
      using AnsatzFunctionType     = typename Ikarus::LocalBasis<DuneBasis>::AnsatzFunctionType;
      /** \brief Type for the points for evaluation, usually the integration points */
      using DomainType             = typename DuneBasis::Traits::DomainType;
      /** \brief Matrix to transform the ansatz function Jacobian to world coordinates*/
      using TransformMatrix        = Eigen::Matrix<ctype, gridDim, gridDim>;
      /** \brief Type for a column of the Jacobian matrix */
      using JacobianColType        = typename Eigen::internal::plain_col_type<Jacobian>::type;
  };

}  // namespace Ikarus
