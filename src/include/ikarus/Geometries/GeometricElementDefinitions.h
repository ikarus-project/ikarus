//
// Created by Alex on 09.06.2021.
//

#pragma once
#include "SimpleGeometry.h"

#include <ikarus/AnsatzFunctions/Lagrange.h>

namespace Ikarus::Geometry {
  template <typename ScalarType, int K>
  using BrickGeometryK = Geometry::SimpleGeometry<ScalarType, 3, 3, Ikarus::LagrangeCube<ScalarType, 3, K>>;

  template <typename ScalarType>
  using LinearBrickGeometry = Geometry::SimpleGeometry<ScalarType, 3, 3, Ikarus::LagrangeCube<ScalarType, 3, 1>>;

  template <typename ScalarType>
  using QuadraticBrickGeometry = Geometry::SimpleGeometry<ScalarType, 3, 3, Ikarus::LagrangeCube<ScalarType, 3, 2>>;

  template <typename ScalarType, int K>
  using PlaneGeometryK = Geometry::SimpleGeometry<ScalarType, 2, 2, Ikarus::LagrangeCube<ScalarType, 2, K>>;

  template <typename ScalarType>
  using LinearPlaneGeometry = Geometry::SimpleGeometry<ScalarType, 2, 2, Ikarus::LagrangeCube<ScalarType, 2, 1>>;

  template <typename ScalarType>
  using QuadraticPlaneGeometry = Geometry::SimpleGeometry<ScalarType, 2, 2, Ikarus::LagrangeCube<ScalarType, 2, 2>>;

  // template<typename ScalarType>
  // using LinearBrickGeometry = SimpleGeometry<ScalarType, 3, 3,Ikarus::LagrangeCube<ScalarType,>;
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