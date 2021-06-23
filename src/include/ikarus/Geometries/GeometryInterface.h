//
// Created by ac120950 on 17.07.2020.
//

#pragma once

#include <Eigen/Core>
namespace Ikarus::Concepts {
  template <class GeometryEntityType>
  concept GeometryEntity = requires(
      GeometryEntityType&& geoEntity,
      Eigen::Matrix<typename GeometryEntityType::ctype, Eigen::Dynamic, GeometryEntityType::mydimension> dN,
      Eigen::Matrix<typename GeometryEntityType::ctype, GeometryEntityType::coorddimension, Eigen::Dynamic> x) {
    typename GeometryEntityType::ctype;
    GeometryEntityType::coorddimension;
    GeometryEntityType::mydimension;
    typename GeometryEntityType::LocalCoordinate;
    typename GeometryEntityType::GlobalCoordinate;
    typename GeometryEntityType::JacobianTransposed;
    typename GeometryEntityType::JacobianInverseTransposed;
    { geoEntity.determinantJacobian(dN, x) } -> std::same_as<typename GeometryEntityType::ctype>;
    { geoEntity.jacobianTransposed(dN, x) } -> std::same_as<typename GeometryEntityType::JacobianTransposed>;
    {
      geoEntity.jacobianInverseTransposed(dN, x)
      } -> std::same_as<typename GeometryEntityType::JacobianInverseTransposed>;
  };
}  // namespace Ikarus::Concepts
