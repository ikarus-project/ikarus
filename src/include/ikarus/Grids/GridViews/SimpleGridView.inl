//
// Created by Alex on 25.05.2021.
//

#include <ranges>

#include <ikarus/Grids/GridInterface.h>

namespace Ikarus::Grid {
  template <int dim, int dimworld, Ikarus::Concepts::Grid GridType>
  auto volumes(SimpleGridView<dim, dimworld, GridType> &gridView) {
    return std::span(gridView.template begin<0>(), gridView.template end<0>());
  }

  template <int dim, int dimworld, Ikarus::Concepts::Grid GridType>
  requires requires { dim > 1; }
  auto edges(SimpleGridView<dim, dimworld, GridType> &gridView) {
    return std::span(gridView.template begin<dim - 1>(), gridView.template end<dim - 1>());
  }

  template <int dim, int dimworld, Ikarus::Concepts::Grid GridType>
  requires requires { dim > 2; }
  auto surfaces(SimpleGridView<dim, dimworld, GridType> &gridView) {
    return std::span(gridView.template begin<dim - 2>(), gridView.template end<dim - 2>());
  }

  template <int dim, int dimworld, Ikarus::Concepts::Grid GridType>
  auto vertices(SimpleGridView<dim, dimworld, GridType> &gridView) {
    return std::span(gridView.template begin<dim>(), gridView.template end<dim>());
  }

  template <int dim, int dimworld, Ikarus::Concepts::Grid GridType, size_t dimE>
  requires requires { dim >= dimE; }
  auto entities(SimpleGridView<dim, dimworld, GridType> &gridView, Dune::index_constant<dimE> &&) {
    return std::span(gridView.template begin<dim - dimE>(), gridView.template end<dim - dimE>());
  }

  template <int dim, int dimworld, Ikarus::Concepts::Grid GridType>
  auto rootEntities(SimpleGridView<dim, dimworld, GridType> &gridView) {
    return std::span(gridView.template begin<0>(), gridView.template end<0>());
  }

}  // namespace Ikarus::Grid