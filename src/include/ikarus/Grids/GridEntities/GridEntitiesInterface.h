//
// Created by Alex on 26.05.2021.
//

#pragma once
#include <dune/geometry/type.hh>

namespace Ikarus::Concepts {
  template <typename GridEntityType>
  concept GridEntity = requires(GridEntityType gEntity, unsigned int codim) {
    { gEntity.level() } -> std::same_as<int>;
    gEntity.geometry();
    { gEntity.type() } -> std::same_as<Dune::GeometryType>;
    { gEntity.subEntities(codim) } -> std::same_as<unsigned int>;
    GridEntityType::codimension;
    GridEntityType::dimension;
    GridEntityType::mydimension;
    typename GridEntityType::Geometry;
  };
}  // namespace Ikarus::Concepts
