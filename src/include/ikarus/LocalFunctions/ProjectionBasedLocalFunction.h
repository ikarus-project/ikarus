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
    //    /** \brief Dimension of the coeffs */
    static constexpr int correctionSize = Traits::correctionSize;
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
    /** \brief Type for the derivatives wrT the coefficients in the local tangent space base*/
    using CoeffDerivMatrix = typename Traits::CoeffDerivMatrix;
    /** \brief Type for the derivatives wrT the coefficients in the embedding space*/
    using CoeffDerivEukMatrix = typename Traits::CoeffDerivEukMatrix;
    /** \brief Type for the derivatives wrT the coefficients in the embedding space on the left and in the tangent space basis on the right*/
    using CoeffDerivRieEukMatrix = Eigen::Matrix<ctype, correctionSize, valueSize>;
    /** \brief Type for ansatz function values */
    using AnsatzFunctionType = typename Traits::AnsatzFunctionType;
    /** \brief Type for the Jacobian of the ansatz function values */
    using AnsatzFunctionJacobian = typename Traits::AnsatzFunctionJacobian;

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

    CoeffDerivRieEukMatrix evaluateDerivativeWRTCoeffsImpl(const AnsatzFunctionType& N, [[maybe_unused]] const AnsatzFunctionJacobian& dN,
                                         int coeffsIndex) const {
      return (coeffs[coeffsIndex].orthonormalFrame().transpose()*evaluateDerivativeWRTCoeffsEukImpl(N,dN,coeffsIndex)).eval();
    }


    CoeffDerivEukMatrix evaluateDerivativeWRTCoeffsEukImpl(const AnsatzFunctionType& N, [[maybe_unused]] const AnsatzFunctionJacobian&,
                                         int coeffsIndex) const {
      GlobalE valE = evaluateEmbeddingFunctionImpl(N);
      return (tryToCallDerivativeOfProjectionWRTposition(valE) * N[coeffsIndex]).eval();
    }

    CoeffDerivMatrix evaluateSecondDerivativeWRTCoeffs(const AnsatzFunctionType& N,
                                                       [[maybe_unused]] const AnsatzFunctionJacobian& dN,
                                                       const AlongType& along,
                                                       const std::array<size_t, 2>& coeffsIndex) const {
      const GlobalE valE = evaluateEmbeddingFunctionImpl(N);

      CoeffDerivEukMatrix ddt = tryToCallSecondDerivativeOfProjectionWRTposition(valE, along) * N[coeffsIndex[0]] * N[coeffsIndex[1]];

      if (coeffsIndex[0]==coeffsIndex[1]) { // Riemannian Hessian Weingarten map correction
        const CoeffDerivEukMatrix dt =evaluateDerivativeWRTCoeffsEukImpl(N,dN,coeffsIndex[0]);
          ddt -=  (coeffs[coeffsIndex[0]].getValue().dot(dt*along)) * CoeffDerivEukMatrix::Identity();
      }

      return coeffs[coeffsIndex[0]].orthonormalFrame().transpose()* ddt*coeffs[coeffsIndex[1]].orthonormalFrame();
    }

    std::array<CoeffDerivRieEukMatrix, gridDim> evaluateDerivativeWRTCoeffsANDSpatialImpl(const AnsatzFunctionType& N,
                                                   [[maybe_unused]] const AnsatzFunctionJacobian& dN,
                                                   int coeffsIndex) const {
      std::array<CoeffDerivEukMatrix, gridDim> WarrayEuk = evaluateDerivativeWRTCoeffsANDSpatialEukImpl(N,dN,coeffsIndex);
      std::array<CoeffDerivRieEukMatrix, gridDim> WarrayRie;
      const auto BLAT = coeffs[coeffsIndex].orthonormalFrame().transpose();
      for (int dir = 0; dir < gridDim; ++dir) {
        WarrayRie[dir] =  BLAT*WarrayEuk[dir];
      }
      return WarrayRie;
    }

    auto evaluateDerivativeWRTCoeffsANDSpatialEukImpl(const AnsatzFunctionType& N,
                                                   [[maybe_unused]] const AnsatzFunctionJacobian& dN,
                                                   int coeffsIndex) const {
      const GlobalE valE        = evaluateEmbeddingFunctionImpl(N);
      const Jacobian J          = evaluateEmbeddingJacobianImpl(dN);
      const CoeffDerivEukMatrix Pm = tryToCallDerivativeOfProjectionWRTposition(valE);
      std::array<CoeffDerivEukMatrix, gridDim> Warray;
      for (int dir = 0; dir < gridDim; ++dir) {
        const auto Qi = tryToCallSecondDerivativeOfProjectionWRTposition(valE, J.col(dir));
        Warray[dir]   = Qi * N[coeffsIndex] + Pm * dN(coeffsIndex, dir);
      }

      return Warray;
    }

    CoeffDerivRieEukMatrix evaluateDerivativeWRTCoeffsANDSpatialSingleImpl(const AnsatzFunctionType& N,
                                                         [[maybe_unused]] const AnsatzFunctionJacobian& dN,
                                                         int coeffsIndex, const int spatialIndex) const {
      const CoeffDerivEukMatrix WEuk = evaluateDerivativeWRTCoeffsANDSpatialSingleEukImpl(N,dN,coeffsIndex,spatialIndex);

      return coeffs[coeffsIndex].orthonormalFrame().transpose()*WEuk;
    }

    CoeffDerivEukMatrix evaluateDerivativeWRTCoeffsANDSpatialSingleEukImpl(const AnsatzFunctionType& N,
                                                                     [[maybe_unused]] const AnsatzFunctionJacobian& dN,
                                                                     int coeffsIndex, const int spatialIndex) const {
      const GlobalE valE         = evaluateEmbeddingFunctionImpl(N);
      const JacobianColType Jcol = evaluateEmbeddingJacobianColImpl(dN, spatialIndex);
      const CoeffDerivEukMatrix Pm  = tryToCallDerivativeOfProjectionWRTposition(valE);
      CoeffDerivEukMatrix W;
      const auto Qi = tryToCallSecondDerivativeOfProjectionWRTposition(valE, Jcol);
      W             = Qi * N[coeffsIndex] + Pm * dN(coeffsIndex, spatialIndex);

      return W;
    }


    auto evaluateThirdDerivativeWRTCoeffsTwoTimesAndSpatialImpl(const AnsatzFunctionType& N,
                                                                [[maybe_unused]] const AnsatzFunctionJacobian& dN,
                                                                const AlongType& along,
                                                                const std::array<size_t, 2>& coeffsIndex) const {
      const GlobalE valE       = evaluateEmbeddingFunctionImpl(N);
      const Jacobian J         = evaluateEmbeddingJacobianImpl(dN);
      const CoeffDerivEukMatrix S = tryToCallSecondDerivativeOfProjectionWRTposition(valE, along);
      std::array<CoeffDerivEukMatrix, gridDim> ChiArrayEuk;
      for (int i = 0; i < gridDim; ++i) {
        const auto chi    = tryToCallThirdDerivativeOfProjectionWRTposition(valE, along, J.col(i));
        const auto& NI    = N[coeffsIndex[0]];
        const auto& NJ    = N[coeffsIndex[1]];
        const auto& dNIdi = dN(coeffsIndex[0], i);
        const auto& dNJdi = dN(coeffsIndex[1], i);
        ChiArrayEuk[i]       = chi * NI * NJ + S * (dNIdi * NJ + dNJdi * NI);
      }
      if (coeffsIndex[0]==coeffsIndex[1]) { // Riemannian Hessian Weingarten map correction
        std::array<CoeffDerivEukMatrix, gridDim> Warray =evaluateDerivativeWRTCoeffsANDSpatialEukImpl(N,dN,coeffsIndex[0]);
        for (int i = 0; i < gridDim; ++i) {
          ChiArrayEuk[i] -=  coeffs[coeffsIndex[0]].getValue().dot(Warray[i]*along) * CoeffDerivEukMatrix::Identity();
        }
      }

      std::array<CoeffDerivMatrix, gridDim> ChiArrayRie;
      const auto BLA0T = coeffs[coeffsIndex[0]].orthonormalFrame().transpose();
      const auto BLA1 = coeffs[coeffsIndex[1]].orthonormalFrame();
      for (int i = 0; i < gridDim; ++i)
        ChiArrayRie[i] = BLA0T*ChiArrayEuk[i]*BLA1;
      return ChiArrayRie;
    }

    CoeffDerivMatrix evaluateThirdDerivativeWRTCoeffsTwoTimesAndSpatialSingleImpl(const AnsatzFunctionType& N,
                                                                      [[maybe_unused]] const AnsatzFunctionJacobian& dN,
                                                                      const AlongType& along,
                                                                      const std::array<size_t, 2>& coeffsIndex,
                                                                      const int spatialIndex) const {
      const GlobalE valE       = evaluateEmbeddingFunctionImpl(N);
      const Jacobian J         = evaluateEmbeddingJacobianImpl(dN);
      const CoeffDerivEukMatrix S = tryToCallSecondDerivativeOfProjectionWRTposition(valE, along);
      CoeffDerivEukMatrix Chi;
      const auto chi    = tryToCallThirdDerivativeOfProjectionWRTposition(valE, along, J.col(spatialIndex));
      const auto& NI    = N[coeffsIndex[0]];
      const auto& NJ    = N[coeffsIndex[1]];
      const auto& dNIdi = dN(coeffsIndex[0], spatialIndex);
      const auto& dNJdi = dN(coeffsIndex[1], spatialIndex);
      Chi               = chi * NI * NJ + S * (dNIdi * NJ + dNJdi * NI);

      if (coeffsIndex[0]==coeffsIndex[1]) { // Riemannian Hessian Weingarten map correction
        CoeffDerivEukMatrix W =evaluateDerivativeWRTCoeffsANDSpatialSingleEukImpl(N,dN,coeffsIndex[0],spatialIndex);
        Chi -=  coeffs[coeffsIndex[0]].getValue().dot(W*along) * CoeffDerivEukMatrix::Identity();

      }
      return coeffs[coeffsIndex[0]].orthonormalFrame().transpose()*Chi*coeffs[coeffsIndex[1]].orthonormalFrame();;
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
    using ctype = typename CoeffContainer::value_type::ctype;
    /** \brief Dimension of the coeffs */
    static constexpr int valueSize = CoeffContainer::value_type::valueSize;
    /** \brief Dimension of the correction size of coeffs */
    static constexpr int correctionSize = CoeffContainer::value_type::correctionSize;
    /** \brief Dimension of the grid */
    static constexpr int gridDim = Ikarus::LocalBasis<DuneBasis>::gridDim;
    /** \brief Type for the return value */
    using FunctionReturnType = typename CoeffContainer::value_type;
    /** \brief Type for the Jacobian matrix */
    using Jacobian = Eigen::Matrix<ctype, valueSize, gridDim>;
    /** \brief Type for the derivatives wrt. the coeffiecients */
    using CoeffDerivMatrix = Eigen::Matrix<ctype, correctionSize, correctionSize>;
    /** \brief Type for the derivatives wrt. the coeffiecients */
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
