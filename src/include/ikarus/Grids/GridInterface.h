//
// Created by Alex on 25.05.2021.
//

#pragma once
#include <concepts>
#include <ranges>

/**
 * \namespace<Ikarus::Grid> The grid namespace contains all grid entities, grid factories and grids
 */

namespace Ikarus::Concepts {
  template <typename GridType>
  concept Grid = requires(GridType grid) {
    {1 + 1 == 2};
    //          typename GridType::ctype;
    //          typename GridType::GridView;
    //          GridType::dimension;
    //          GridType::dimensionworld;
    //          typename GridType::VertexType;
    //          { grid.leafGridView() } ->  std::same_as<typename GridType::GridView>;
    //            { grid.getNextFreeId() } -> std::integral;
    //            { elements(grid) }  ->  std::ranges::range;
    //            { vertices(grid) }  ->  std::ranges::range;
  };
}  // namespace Ikarus::Concepts