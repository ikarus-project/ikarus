//
// Created by Alex on 12.06.2021.
//

#pragma once
#include <ranges>

namespace Ikarus::Concepts {
  template <typename GridType>
  concept GridView = requires(GridType grid&& int coDim) {
    { grid.begin<coDim>() } -> std::r;
    { grid.begin<coDim>() } -> std::ranges::range;
    //            { grid.getNextFreeId() } -> std::integral;
    //            { elements(grid) }  ->  std::ranges::range;
    //            { vertices(grid) }  ->  std::ranges::range;
  };
}  // namespace Ikarus::Concepts