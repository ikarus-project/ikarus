//
// Created by Alex on 12.06.2021.
//

#pragma once
#include <ikarus/Grids/GridViews/IndexSets/SimpleGridIndexSet.h>
namespace Ikarus::Grid {

  template <int dim, int dimworld>
  class SimpleGrid;

  template <int dim, int dimworld>
  class SimpleGridView {
  public:
    using GridType     = SimpleGrid<dim, dimworld>;
    using IndexSetType = SimpleGridIndexSet<GridType>;

    explicit SimpleGridView(GridType& gridInput, int levelInput = 0)
        : grid{&gridInput},
          level{levelInput},
          indexSet_{std::make_unique<SimpleGridIndexSet<GridType>>(gridInput, levelInput)} {}

    static constexpr int dimension      = dim;
    static constexpr int dimensionworld = dimworld;

    template <int coDim>
    auto begin() {
      return grid->gridEntitiesContainer.template getSubEntities<coDim>(0).begin();
    }

    template <int coDim>
    auto& getEntities() {
      return grid->gridEntitiesContainer.template getSubEntities<coDim>(0);
    }

    template <int coDim>
    const auto& getEntities() const {
      return grid->gridEntitiesContainer.template getSubEntities<coDim>(0);
    }

    template <int coDim>
    auto end() {
      return grid->gridEntitiesContainer.template getSubEntities<coDim>(0).end();
    }

    const auto& indexSet() { return *indexSet_; }

    int size(int codim) const;

  private:
    GridType* grid;
    int level;
    const std::shared_ptr<IndexSetType> indexSet_;
  };
}  // namespace Ikarus::Grid

#include "SimpleGridView.inl"