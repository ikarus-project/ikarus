//
// Created by Alex on 25.05.2021.
//

#pragma once
#include <dune/geometry/multilineargeometry.hh>

#include <ikarus/Grids/GridEntities/DefaultGridEntities.h>
#include <ikarus/Grids/GridViews/SimpleGridView.h>

namespace Ikarus::Grid {
  template <int dim, int dimworld>
  class SimpleGrid;
  namespace Impl {
    template <int griddim, int wdim, int... codim>
    static std::tuple<std::vector<DefaultGridEntity<griddim, codim, wdim>>...> GridEntityTupleGenerator(
        std::integer_sequence<int, codim...>);
  }

  template <int dim, int dimworld>
  class SimpleGrid {
    static_assert(dim <= dimworld, "The dimension of the grid can not be larger then the embedding space.");

  public:
    using ctype                         = double;
    static constexpr int dimension      = dim;
    static constexpr int dimensionworld = dimworld;
    using GridEntityTuple
        = decltype(Impl::GridEntityTupleGenerator<dim, dimworld>(std::make_integer_sequence<int, dimension + 1>()));
    using VertexType = DefaultGridEntity<dim, dim, dimworld>;
    using RootEntity = DefaultGridEntity<dim, 0, dimworld>;

    auto leafGridView() { return SimpleGridView<dimension, dimensionworld, SimpleGrid>(*this, 0); }

  private:
    auto& getElements(int level = 0) { return getSubEntities<0>(level); }

    auto& getVertices(int level = 0) {
      return getSubEntities<dim>(level);
      ;
    }

    auto& getEdges(int level = 0) { return getSubEntities<dim - 1>(level); }

    auto& getSurfaces(int level = 0) requires requires { dim - 2 > 0; }
    { return getSubEntities<dim - 2>(level); }

    template <int subEnt>
    auto& getSubEntities(int level = 0) {
      static_assert(subEnt >= 0 && subEnt <= dim, "You asked for a non-existing subEntity!");
      return std::get<subEnt>(gridEntities[level]);
    }

    unsigned int getNextFreeId() { return freeIdCounter++; }

    template <Concepts::Grid GridType>
    friend class SimpleGridFactory;
    friend class SimpleGridView<dim, dimworld, SimpleGrid>;

    std::vector<GridEntityTuple> gridEntities;

    unsigned int freeIdCounter{};
  };

}  // namespace Ikarus::Grid

#include <ikarus/Grids/SimpleGrid/SimpleGridFactory.h>
