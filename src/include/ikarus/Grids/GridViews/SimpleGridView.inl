//
// Created by Alex on 25.05.2021.
//

#include <ranges>

#include <ikarus/Grids/GridInterface.h>

namespace Ikarus::Grid {
  template <int dim, int dimworld, Ikarus::Concepts::Grid GridType>
  auto elements(SimpleGridView<dim, dimworld, GridType> &gridView) {
    return std::span(gridView.template begin<0>(), gridView.template end<0>());
  }

  template <int dim, int dimworld, Ikarus::Concepts::Grid GridType, std::enable_if_t<(dim > 1), bool> = true>
  auto edges(SimpleGridView<dim, dimworld, GridType> &gridView) {
    return std::span(gridView.template begin<dim - 1>(), gridView.template end<dim - 1>());
  }

  template <int dim, int dimworld, Ikarus::Concepts::Grid GridType>
  auto vertices(SimpleGridView<dim, dimworld, GridType> &gridView) {
    return std::span(gridView.template begin<dim>(), gridView.template end<dim>());
  }

}  // namespace Ikarus::Grid