//
// Created by Alex on 25.05.2021.
//

#pragma once
#include <concepts>

#include <dune/geometry/type.hh>

namespace Ikarus::Concepts {
  template <typename GridFactoryType>
  concept GridFactory = requires(GridFactoryType fac, typename GridFactoryType::CoordinateType vertexPos) {
    typename GridFactoryType::ctype;
    GridFactoryType::dimension;
    GridFactoryType::dimensionworld;
    typename GridFactoryType::VertexType;
    { fac.insertVertex(vertexPos) } -> std::same_as<void>;
    //        { fac.insertElement(std::span<vertexPos>)  }  ->  std::same_as< void>;
  };
}  // namespace Ikarus::Concepts