//
// Created by Alex on 21.04.2021.
//

#pragma once

#include <concepts>
#include <span>

#include <Eigen/Core>
#include <Eigen/Dense>

#include <ikarus/Geometries/GeometryWithExternalInput.h>
#include <ikarus/Interpolators/Interpolator.h>

namespace Ikarus::Geometry {

  template <typename ct, int wdim, int geodim, typename ShapeFunctionType>
  requires requires { geodim <= wdim; }
  class SimpleGeometry {
  private:
    using ExternalGeometry = GeometryWithExternalInput<ct, wdim, geodim>;

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

    /** \brief Type for container of the vertex coordinates
     *
     * \note E.g. if the size of the shapefunctions are know at compile time we use a fixed size
     * Matrix If we have a Q1 element in 3D space we have FixedMatrixd<3,4> If we have an element
     * with unknown shape function we have FixedMatrixd<3,Dynamic>
     * */
    using VertexContainerType =
        typename std::conditional<ShapeFunctionType::sizeOfShapeFunctions == Eigen::Dynamic,
                                  Eigen::Matrix<double,mydimension,Eigen::Dynamic>,
                                  Eigen::Matrix<double,mydimension, ShapeFunctionType::sizeOfShapeFunctions> >::type;

    SimpleGeometry(const VertexContainerType& verticesInput) : vertices{verticesInput} {}

    /** \brief Type for the transposed Jacobian matrix */
    using JacobianTransposedType = Eigen::Matrix<ctype, mydimension, coorddimension>;

    /** \brief Type for the transposed inverse Jacobian matrix */
    using JacobianInverseTransposed = Eigen::Matrix<ctype, coorddimension, mydimension>;

    ctype determinantJacobian(const LocalCoordinate& x) {
      return ExternalGeometry::determinantJacobian(shapeFunctions.evaluateJacobian(x), vertices);
    };

    JacobianTransposedType jacobianTransposed(const LocalCoordinate& x) {
      return ExternalGeometry::jacobianTransposed(shapeFunctions.evaluateJacobian(x), vertices);
    }

    JacobianInverseTransposed jacobianInverseTransposed(const LocalCoordinate& x) {
      return ExternalGeometry::jacobianInverseTransposed(shapeFunctions.evaluateJacobian(x), vertices);
    }

    typename ShapeFunctionType::JacobianType transformCurvLinearDerivativesToCartesian(const LocalCoordinate& x) {
      return ExternalGeometry::transformCurvLinearDerivativesToCartesian(shapeFunctions.evaluateJacobian(x), vertices);
    }

  private:
    ShapeFunctionType shapeFunctions;
    VertexContainerType vertices;
  };

  // template<typename ScalarType>
  // using LinearBrickGeometry = SimpleGeometry<ScalarType, 3, 3,>;
  // template<typename ScalarType>
  // using SurfaceGeometry = SimpleGeometry<ScalarType, 3, 2>;
  // template<typename ScalarType>
  // using PlaneGeometry = SimpleGeometry<ScalarType, 2, 2>;
  // template<typename ScalarType>
  // using Curve3dGeometry = SimpleGeometry<ScalarType, 3, 1>;
  // template<typename ScalarType>
  // using Curve2dGeometry = SimpleGeometry<ScalarType, 2, 1>;
  // template<typename ScalarType>
  // using Curve1dGeometry = SimpleGeometry<ScalarType, 1, 1>;
}  // namespace Ikarus::Geometry
