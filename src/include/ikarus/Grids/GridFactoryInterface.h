//
// Created by Alex on 25.05.2021.
//

#pragma once
#include <concepts>

#include <dune/geometry/type.hh>

namespace Ikarus::Concepts {
  template <typename GridFactoryType>
  concept GridFactory = requires(GridFactoryType fac, typename GridFactoryType::CoordinateType vertexPos) {
    { fac.insertVertex(vertexPos) } -> std::same_as<void>;
    { fac.insertElement(std::span<size_t>) } -> std::same_as<void>;
    { fac.createGrid(std::span<size_t>) } -> std::same_as<GridFactoryType::GridType>;
  };
}  // namespace Ikarus::Concepts