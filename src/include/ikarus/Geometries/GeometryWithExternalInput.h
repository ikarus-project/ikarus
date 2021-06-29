//
// Created by Alex on 21.04.2021.
//

#pragma once

#include <concepts>
#include <span>

#include <Eigen/Core>
#include <Eigen/Dense>

#include "ikarus/Interpolators/Interpolator.h"

namespace Ikarus::Geometry {

  template <typename ct, int wdim, int geodim, std::enable_if_t<geodim <= wdim, bool> = true>
  class GeometryWithExternalInput {
  public:
    /** \brief Type used for coordinates */
    using ctype = ct;

    /** \brief Dimension of the world space */
    static constexpr int coorddimension = wdim;

    /** \brief Dimension of the geometry */
    static constexpr int mydimension = geodim;

    /** \brief Type for local coordinate vector */
    using LocalCoordinate = Eigen::Matrix<ctype, mydimension, 1>;

    /** \brief Type for coordinate vector in world space */
    using GlobalCoordinate = Eigen::Matrix<ctype, coorddimension, 1>;

    /** \brief Type for the transposed Jacobian matrix */
    using JacobianTransposed = Eigen::Matrix<ctype, mydimension, coorddimension>;

    /** \brief Type for the transposed inverse Jacobian matrix */
    using JacobianInverseTransposed = Eigen::Matrix<ctype, coorddimension, mydimension>;

    template <typename DerivedAnsatzFunctionType, typename GlobalCoordinateListType>
    static ctype determinantJacobian(const Eigen::MatrixBase<DerivedAnsatzFunctionType>& dN,
                                     const Eigen::MatrixBase<GlobalCoordinateListType>& nodevalueList) {
      const auto JT = jacobianTransposed(dN, nodevalueList);
      return sqrt((JT * JT.transpose()).determinant());
    }

    template <typename DerivedAnsatzFunctionType, typename GlobalCoordinateListType>
    static JacobianTransposed jacobianTransposed(const Eigen::MatrixBase<DerivedAnsatzFunctionType>& dN,
                                                 const Eigen::MatrixBase<GlobalCoordinateListType>& nodevalueList) {
      static_assert(DerivedAnsatzFunctionType::ColsAtCompileTime == mydimension);
      static_assert(GlobalCoordinateListType::RowsAtCompileTime == coorddimension);
      spdlog::info("{} {}", dN.rows(), nodevalueList.cols());
      assert(dN.rows() == nodevalueList.cols());
      JacobianTransposed JT;
      for (int i = 0; i < JT.rows(); ++i)
        JT.row(i) = interpolate(dN.col(i), nodevalueList).transpose();

      return JT;
    }

    template <typename DerivedAnsatzFunctionType, typename GlobalCoordinateListType>
    static JacobianInverseTransposed jacobianInverseTransposed(
        const Eigen::MatrixBase<DerivedAnsatzFunctionType>& dN,
        const Eigen::MatrixBase<GlobalCoordinateListType>& nodevalueList) {
      return jacobianTransposed(dN, nodevalueList).completeOrthogonalDecomposition().pseudoInverse();
    }
  };

  template <typename ScalarType>
  using ExternalBrickGeometry = GeometryWithExternalInput<ScalarType, 3, 3>;
  template <typename ScalarType>
  using ExternalSurfaceGeometry = GeometryWithExternalInput<ScalarType, 3, 2>;
  template <typename ScalarType>
  using ExternalPlaneGeometry = GeometryWithExternalInput<ScalarType, 2, 2>;
  template <typename ScalarType>
  using ExternalCurve3dGeometry = GeometryWithExternalInput<ScalarType, 3, 1>;
  template <typename ScalarType>
  using ExternalCurve2dGeometry = GeometryWithExternalInput<ScalarType, 2, 1>;
  template <typename ScalarType>
  using ExternalCurve1dGeometry = GeometryWithExternalInput<ScalarType, 1, 1>;
}  // namespace Ikarus::Geometry
