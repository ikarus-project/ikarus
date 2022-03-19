//
// Created by Alex on 21.04.2021.
//

#pragma once

#include "LocalFunctionInterface.h"

#include <concepts>
#include <iostream>

#include <Eigen/Core>
#include <Eigen/Dense>

#include "ikarus/Interpolators/Interpolator.h"
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
    //    using DomainType = typename Traits::DomainType;
    /** \brief Type used for coordinates */
    using ctype = typename Traits::ctype;
    //    /** \brief Dimension of the coeffs */
    static constexpr int manifoldEmbeddingDim = Traits::manifoldEmbeddingDim;

    /** \brief Dimension of the grid */
    static constexpr int gridDim = Traits::gridDim;

    /** \brief Type for coordinate vector in world space */
    using Manifold  = typename Traits::FunctionReturnType;
    using AlongType = Eigen::Vector<ctype, manifoldEmbeddingDim>;
    using GlobalE   = typename Manifold::CoordinateType;
    /** \brief Type for the transposed Jacobian matrix */
    using Jacobian               = typename Traits::Jacobian;
    using JacobianColType        = typename Traits::JacobianColType;
    using FieldMat               = typename Traits::FieldMat;
    using AnsatzFunctionType     = typename Traits::AnsatzFunctionType;
    using AnsatzFunctionJacobian = typename Traits::AnsatzFunctionJacobian;
    using TransformMatrix        = typename Traits::TransformMatrix;

    void setCoefficients(const CoeffContainer& otherCoeffs) { coeffs = otherCoeffs; }

  private:
    Jacobian evaluateDerivativeWRTSpaceAllImpl(const AnsatzFunctionType& N, const AnsatzFunctionJacobian& dN) const {
      Jacobian J   = evaluateEmbeddingJacobianImpl(dN);
      GlobalE valE = evaluateEmbeddingFunctionImpl(N);
      return Manifold::derivativeOfProjectionWRTposition(valE) * J;
    }

    JacobianColType evaluateDerivativeWRTSpaceSingleImpl(const AnsatzFunctionType& N, const AnsatzFunctionJacobian& dN,
                                                         int spaceIndex) const {
      JacobianColType Jcol = evaluateEmbeddingJacobianColImpl(dN, spaceIndex);
      GlobalE valE         = evaluateEmbeddingFunctionImpl(N);
      return Manifold::derivativeOfProjectionWRTposition(valE) * Jcol;
    }

    auto evaluateDerivativeWRTCoeffsImpl(const AnsatzFunctionType& N, [[maybe_unused]] const AnsatzFunctionJacobian&,
                                         int coeffsIndex) const {
      GlobalE valE = evaluateEmbeddingFunctionImpl(N);
      return (Manifold::derivativeOfProjectionWRTposition(valE) * N[coeffsIndex]).eval();
    }

    auto evaluateSecondDerivativeWRTCoeffs(const long unsigned gpIndex, const Manifold& val, const AlongType& along,
                                           const std::array<size_t, gridDim>& coeffsIndex) const {
      const GlobalE valE = evaluateEmbeddingFunctionImpl(basis.getFunction(gpIndex));
      FieldMat Snn       = Manifold::secondDerivativeOfProjectionWRTposition(valE, along)
                     * basis.getFunction(gpIndex)[coeffsIndex[0]] * basis.getFunction(gpIndex)[coeffsIndex[1]];

      return Snn;
    }

    auto evaluateDerivativeWRTCoeffsANDSpatialImpl(const AnsatzFunctionType& N,
                                                   [[maybe_unused]] const AnsatzFunctionJacobian& dN,
                                                int coeffsIndex) const {
      const GlobalE valE = evaluateEmbeddingFunctionImpl(N);
      const Jacobian J   = evaluateEmbeddingJacobianImpl(dN);
      const FieldMat Pm  = Manifold::derivativeOfProjectionWRTposition(valE);
      std::array<FieldMat, gridDim> Warray;
      for (int dir = 0; dir < gridDim; ++dir) {
        const auto Qi = Manifold::secondDerivativeOfProjectionWRTposition(valE, J.col(dir));
        Warray[dir]   = Qi * N[coeffsIndex] + Pm * dN(coeffsIndex, dir);
      }

      return Warray;
    }

    auto evaluateDerivativeWRTCoeffsANDSpatialSingleImpl(const AnsatzFunctionType& N,
                                                         [[maybe_unused]] const AnsatzFunctionJacobian& dN,
                                                         int coeffsIndex,
                                                         const int spatialIndex) const {
      const GlobalE valE  = evaluateEmbeddingFunctionImpl(N);
      const JacobianColType Jcol = evaluateEmbeddingJacobianColImpl(dN, spatialIndex);
      const FieldMat Pm   = Manifold::derivativeOfProjectionWRTposition(valE);
      FieldMat W;
      const auto Qi = Manifold::secondDerivativeOfProjectionWRTposition(valE, Jcol);
      W             = Qi * N[coeffsIndex] + Pm * dN(coeffsIndex, spatialIndex);

      return W;
    }

    auto evaluateThirdDerivativeWRTCoeffsTwoTimesAndSpatialImpl(const long unsigned gpIndex, const Manifold& val,
                                                                const AlongType& along,
                                                                const TransformMatrix& transformJ,
                                                                const std::array<size_t, gridDim>& coeffsIndex) const {
      const AnsatzFunctionJacobian dN = (basis.getJacobian(gpIndex) * transformJ).eval();
      const GlobalE valE              = evaluateEmbeddingFunctionImpl(basis.getFunction(gpIndex));
      const Jacobian J                = evaluateEmbeddingJacobianImpl(dN);
      const FieldMat S                = Manifold::secondDerivativeOfProjectionWRTposition(valE, along);
      std::array<FieldMat, gridDim> ChiArray;
      for (int i = 0; i < gridDim; ++i) {
        const auto chi    = Manifold::thirdDerivativeOfProjectionWRTposition(valE, along, J.col(i));
        const auto& NI    = basis.getFunction(gpIndex)[coeffsIndex[0]];
        const auto& NJ    = basis.getFunction(gpIndex)[coeffsIndex[1]];
        const auto& dNIdi = dN(coeffsIndex[0], i);
        const auto& dNJdi = dN(coeffsIndex[1], i);
        ChiArray[i]       = chi * NI * NJ + S * (dNIdi * NJ + dNJdi * NI);
      }

      return ChiArray;
    }

    Manifold evaluateFunctionImpl(const AnsatzFunctionType& N) const {
      return Manifold(evaluateEmbeddingFunctionImpl(N));
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
    using ctype                               = typename CoeffContainer::value_type::ctype;
    static constexpr int manifoldEmbeddingDim = CoeffContainer::value_type::valueSize;

    static constexpr int gridDim = Ikarus::LocalBasis<DuneBasis>::gridDim;

    using FunctionReturnType = typename CoeffContainer::value_type;

    using Jacobian               = Eigen::Matrix<ctype, manifoldEmbeddingDim, gridDim>;
    using FieldMat               = Eigen::Matrix<ctype, manifoldEmbeddingDim, manifoldEmbeddingDim>;
    using AnsatzFunctionJacobian = typename Ikarus::LocalBasis<DuneBasis>::JacobianType;
    using AnsatzFunctionType     = typename Ikarus::LocalBasis<DuneBasis>::AnsatzFunctionType;
    using DomainType             = typename DuneBasis::Traits::DomainType;
    using TransformMatrix        = Eigen::Matrix<ctype, gridDim, gridDim>;
    using JacobianColType        = typename Eigen::internal::plain_col_type<Jacobian>::type;
  };

}  // namespace Ikarus
