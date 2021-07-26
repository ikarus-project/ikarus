//
// Created by Alex on 25.05.2021.
//

#include <ranges>

namespace Ikarus::Grid {
  template <int dim, int dimworld>
  requires(dim >= 3) auto volumes(SimpleGridView<dim, dimworld> &gridView) {
    return std::span(gridView.template begin<dim - 3>(), gridView.template end<dim - 3>());
  }

  template <int dim, int dimworld>
  requires(dim >= 2) auto surfaces(SimpleGridView<dim, dimworld> &gridView) {
    return std::span(gridView.template begin<dim - 2>(), gridView.template end<dim - 2>());
  }
  template <int dim, int dimworld>
  requires(dim >= 1) auto edges(SimpleGridView<dim, dimworld> &gridView) {
    return std::span(gridView.template begin<dim - 1>(), gridView.template end<dim - 1>());
  }

  template <int dim, int dimworld>
  auto vertices(SimpleGridView<dim, dimworld> &gridView) {
    return std::span(gridView.template begin<dim>(), gridView.template end<dim>());
  }

  template <int dim, int dimworld, size_t dimE>
  requires(dim >= dimE) auto entities(SimpleGridView<dim, dimworld> &gridView,
                                      Dune::index_constant<dimE> &&) {
    return std::span(gridView.template begin<dim - dimE>(), gridView.template end<dim - dimE>());
  }

  template <int dim, int dimworld>
  auto rootEntities(SimpleGridView<dim, dimworld> &gridView) {
    return std::span(gridView.template begin<0>(), gridView.template end<0>());
  }

}  // namespace Ikarus::Grid