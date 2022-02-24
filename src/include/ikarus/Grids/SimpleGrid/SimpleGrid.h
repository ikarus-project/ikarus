//
// Created by Alex on 25.05.2021.
//

#pragma once

#include <dune/geometry/multilineargeometry.hh>
#include <dune/grid/common/exceptions.hh>

#include <ikarus/Grids/GridEntities/DefaultGridEntities.h>
#include <ikarus/Grids/GridViews/SimpleGridView.h>

namespace Ikarus::Grid {
  template <int dimension, int dimensionworld>
  class SimpleGridFactory;

  namespace Impl {
    template <int griddim, int wdim, int... codim>
    static std::tuple<std::vector<DefaultGridEntity<griddim, codim, wdim>>...> GridEntityTupleGenerator(
        std::integer_sequence<int, codim...>);

    template <int dim, int wdim>
    class GridEntitiesContainer {
    public:
      GridEntitiesContainer(size_t vertexsize, size_t elementsize) {
        gridEntities.resize(1);
        getVertices().reserve(vertexsize);
        getRootEntities().reserve(elementsize);
      }

      /** Rule of five: Since the gridEntities have pointer to their subentities we can not allow this to be copyable */
      GridEntitiesContainer(const GridEntitiesContainer&)                = delete;
      GridEntitiesContainer& operator=(const GridEntitiesContainer&)     = delete;
      GridEntitiesContainer& operator=(GridEntitiesContainer&&) noexcept = default;
      GridEntitiesContainer(GridEntitiesContainer&&) noexcept            = default;
      ~GridEntitiesContainer()                                           = default;

      /** \brief  The GridEntityTupleGenerator generates a tuple of grid entities which can be of different size
       * depending on the dimension of the grid, e.g. a grid of dimension 3 has volumes, surfaces, edges and vertices,
       * whereas a grid with dimension 1 only has edges and vertices
       */
      using GridEntityTuple
          = decltype(Impl::GridEntityTupleGenerator<dim, wdim>(std::make_integer_sequence<int, dim + 1>()));

      auto& getRootEntities(int level = 0) { return getSubEntities<0>(level); }
      auto& getVolumes(int level = 0) requires(dim >= 3) { return getSubEntities<dim - 3>(level); }
      auto& getSurfaces(int level = 0) requires(dim >= 2) { return getSubEntities<dim - 2>(level); }
      auto& getEdges(int level = 0) requires(dim >= 1) { return getSubEntities<dim - 1>(level); }
      auto& getVertices(int level = 0) { return getSubEntities<dim>(level); }

      template <int subEnt>
      auto& getSubEntities(int level = 0) {
        static_assert(subEnt >= 0 && subEnt <= dim, "You asked for a non-existing subEntity!");
        return std::get<subEnt>(gridEntities[level]);
      }

    private:
      std::vector<GridEntityTuple> gridEntities;
    };
  }  // namespace Impl

  template <int dim, int dimworld>
  class SimpleGrid {
    static_assert(dim <= dimworld, "The dimension of the grid can not be larger then the embedding space.");

  public:
    using ctype                         = double;
    static constexpr int dimension      = dim;
    static constexpr int dimensionworld = dimworld;
    using GridEntitiesContainer         = Impl::GridEntitiesContainer<dim, dimworld>;

    SimpleGrid(GridEntitiesContainer&& gridEntCont, size_t id)
        : gridEntitiesContainer{std::forward<GridEntitiesContainer>(gridEntCont)}, freeIdCounter{id} {}

    auto leafGridView() { return SimpleGridView<dimension, dimensionworld>(*this, 0); }

  private:
    size_t getNextFreeId() { return freeIdCounter++; }
    friend class SimpleGridView<dimension, dimensionworld>;
    friend class SimpleGridIndexSet<SimpleGrid<dimension, dimensionworld>>;

    template <class EntityType>
    size_t getEntityID(const EntityType& entity) const {
      return entity.getID();
    }

    template <class EntityType>
    size_t getEntityID(EntityType* const entity) const {
      return getEntityID(*entity);
    }

    GridEntitiesContainer gridEntitiesContainer;

    size_t freeIdCounter{};
  };

}  // namespace Ikarus::Grid

#include <ikarus/Grids/SimpleGrid/SimpleGridFactory.h>
