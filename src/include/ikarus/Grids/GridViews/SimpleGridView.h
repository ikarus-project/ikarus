//
// Created by Alex on 12.06.2021.
//

#pragma once

namespace Ikarus::Grid {

  template <int dim, int dimworld>
  class SimpleGrid;

  template <int dim, int dimworld>
  class SimpleGridView {
  public:
    using GridType = SimpleGrid<dim, dimworld>;

    explicit SimpleGridView(GridType& gridInput, int levelInput = 0) : grid{&gridInput}, level{levelInput} {}

    static constexpr int dimension = GridType::dimension;

    template <int coDim>
    auto begin() {
      return grid->gridEntitiesContainer.template getSubEntities<coDim>(0).begin();
    }

    template <int coDim>
    auto end() {
      return grid->gridEntitiesContainer.template getSubEntities<coDim>(0).end();
    }

  private:
    GridType* grid;
    int level;
  };
}  // namespace Ikarus::Grid

#include "SimpleGridView.inl"