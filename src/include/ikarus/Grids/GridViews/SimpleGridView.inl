//
// Created by Alex on 25.05.2021.
//

#include <ranges>

namespace Ikarus::Grid {
  template <int dim, int dimworld>
  requires(dim >= 3) auto& volumes(SimpleGridView<dim, dimworld> &gridView) {
    return gridView.template getEntities<dim - 3>();
  }

  template <int dim, int dimworld>
  requires(dim >= 2) auto& surfaces(SimpleGridView<dim, dimworld> &gridView) {
    return gridView.template getEntities<dim - 2>();
  }

  template <int dim, int dimworld>
  requires(dim >= 2) const auto& surfaces(const SimpleGridView<dim, dimworld> &gridView) {
    return gridView.template getEntities<dim - 2>();
  }

  template <int dim, int dimworld>
  requires(dim >= 1) auto& edges(SimpleGridView<dim, dimworld> &gridView) {
    return gridView.template getEntities<dim - 1>();
  }

  template <int dim, int dimworld>
  requires(dim >= 1) const auto& edges(const SimpleGridView<dim, dimworld> &gridView) {
    return gridView.template getEntities<dim - 1>();
  }

  template <int dim, int dimworld>
  auto& vertices(SimpleGridView<dim, dimworld> &gridView) {
    return gridView.template getEntities<dim >();
  }

  template <int dim, int dimworld, size_t dimE>
  requires(dim >= dimE) auto& entities(SimpleGridView<dim, dimworld> &gridView, Dune::index_constant<dimE> &&) {
    return gridView.template getEntities<dim - dimE>();
  }

  template <int dim, int dimworld>
  auto& rootEntities(SimpleGridView<dim, dimworld> &gridView) {
    return gridView.template getEntities<0>();
  }

  template <int dim, int dimworld>
  int SimpleGridView<dim, dimworld>::size(int codim) const {
    assert(codim <= dim && codim > 0 && "The requested codimension can not be higher than the grid dimension.");
    if constexpr (dim >= 1) {
      if (codim == 0)
        return entities(*this, Dune::index_constant<dim>{}).size();
      else if (codim == 1)
        return entities(*this, Dune::index_constant<dim - 1>{}).size();
    }
    if constexpr (dim >= 2)
      if (codim == 2) return entities(*this, Dune::index_constant<dim - 2>{}).size();

    if constexpr (dim >= 3)
      if (codim == 3) return entities(*this, Dune::index_constant<dim - 3>{}).size();
    return 0;
  }

}  // namespace Ikarus::Grid