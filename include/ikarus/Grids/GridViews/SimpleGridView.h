//
// Created by Alex on 12.06.2021.
//
#
#pragma once
#include <ikarus/Grids/GridInterface.h>
#include <ikarus/Grids/SimpleGrid/SimpleGrid.h>
namespace Ikarus::Grid {

  template <int dim, int dimworld, Ikarus::Concepts::Grid GridType> class SimpleGridView {
  public:
    explicit SimpleGridView(GridType& gridInput, int levelInput = 0)
        : grid{&gridInput}, level{levelInput} {}

    template <int coDim> auto begin() { return grid->template getSubEntities<coDim>(0).begin(); }

    template <int coDim> auto end() { return grid->template getSubEntities<coDim>(0).end(); }

  private:
    GridType* grid;
    int level;
  };
}  // namespace Ikarus::Grid

#include "SimpleGridView.inl"