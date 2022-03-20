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
  class SimpleLocalFunction : public LocalFunctionInterface<SimpleLocalFunction<DuneBasis, CoeffContainer>> {
    using Base = LocalFunctionInterface<SimpleLocalFunction<DuneBasis, CoeffContainer>>;

  public:
    friend Base;
    SimpleLocalFunction(const Ikarus::LocalBasis<DuneBasis>& basis_, const CoeffContainer& coeffs_)
        : basis{basis_}, coeffs{coeffs_}, coeffsAsMat{Ikarus::LinearAlgebra::viewAsEigenMatrixFixedDyn(coeffs)} {}

    using Traits = LocalFunctionTraits<SimpleLocalFunction>;
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

    auto& coefficientsRef() { return coeffs; }

  private:
    Jacobian evaluateDerivativeWRTSpaceAllImpl(const AnsatzFunctionType& N, const AnsatzFunctionJacobian& dN) const {
      return coeffsAsMat * dN;
    }

    JacobianColType evaluateDerivativeWRTSpaceSingleImpl(const AnsatzFunctionType&, const AnsatzFunctionJacobian& dN,
                                                         int spaceIndex) const {
      return coeffsAsMat * dN.col(spaceIndex);
    }

    auto evaluateDerivativeWRTCoeffsImpl(const AnsatzFunctionType& N, [[maybe_unused]] const AnsatzFunctionJacobian&,
                                         int coeffsIndex) const {
      Eigen::DiagonalMatrix<ctype, manifoldEmbeddingDim> mat;
      mat.setIdentity(manifoldEmbeddingDim);
      mat *= N[coeffsIndex];
      return mat;
    }

    auto evaluateDerivativeWRTCoeffsANDSpatialImpl(const AnsatzFunctionType& N,
                                                   [[maybe_unused]] const AnsatzFunctionJacobian& dN,
                                                   int coeffsIndex) const {
      std::array<Eigen::DiagonalMatrix<ctype, manifoldEmbeddingDim>, gridDim> Warray;
      for (int dir = 0; dir < gridDim; ++dir) {
        Warray[dir].setIdentity(manifoldEmbeddingDim);
        Warray[dir] *= dN(coeffsIndex, dir);
      }

      return Warray;
    }

    auto evaluateDerivativeWRTCoeffsANDSpatialSingleImpl(const AnsatzFunctionType& N,
                                                         [[maybe_unused]] const AnsatzFunctionJacobian& dN,
                                                         int coeffsIndex, const int spatialIndex) const {
      Eigen::DiagonalMatrix<ctype, manifoldEmbeddingDim> W;
      W.setIdentity(manifoldEmbeddingDim);
      W *= dN(coeffsIndex, spatialIndex);

      return W;
    }


    Manifold evaluateFunctionImpl(const AnsatzFunctionType& N) const {
      return Manifold(coeffsAsMat * N);
    }


    const Ikarus::LocalBasis<DuneBasis>& basis;
    CoeffContainer coeffs;
    const decltype(Ikarus::LinearAlgebra::viewAsEigenMatrixFixedDyn(coeffs)) coeffsAsMat;
  };

  template <typename DuneBasis, typename CoeffContainer>
  struct LocalFunctionTraits<SimpleLocalFunction<DuneBasis, CoeffContainer>> {
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
