//
// Created by Alex on 25.05.2021.
//

#pragma once
#include <dune/geometry/multilineargeometry.hh>

#include <ikarus/Grids/GridEntities/DefaultGridEntities.h>
#include <ikarus/Grids/GridViews/SimpleGridView.h>
#include <dune/grid/common/exceptions.hh>

#include "SimpleGridTypedefs.h"

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

    SimpleGrid(const std::vector<SimpleGridTypedefs::VertexIndexPair<dimensionworld>>& verticesPositions,
               const std::vector<std::vector<size_t>>& edgesVertexIndices,
               const std::vector<std::vector<size_t>>& surfaceVertexIndices,
               const std::vector<std::vector<size_t>>& elementsVertices,
               const std::vector<std::vector<size_t>>& elementEdgeIndices,
               const std::vector<std::vector<size_t>>& elementSurfaceIndices){
      if (verticesPositions.empty())
        DUNE_THROW(Dune::GridError, "verticesPositions vector is empty. Unable to create Grid");
      if (elementsVertices.empty()) DUNE_THROW(Dune::GridError, "elements vector is empty. Unable to create Grid");

      gridEntities.resize(1);  // resize to hold coarsest grid level

      getVertices().reserve(verticesPositions.size());
      getElements().reserve(elementsVertices.size());

      // add vertices to the grid
      for (auto &vert : verticesPositions) {
        typename SimpleGrid<dimension,dimensionworld>::VertexType newVertex(0, vert.vertex, getNextFreeId());
        newVertex.levelIndex = vert.index;
        getVertices().push_back(newVertex);
      }

      // add element and set vertex pointer of elements
      for (auto &eleVertices : elementsVertices) {
        typename SimpleGrid<dimension,dimensionworld>::RootEntity newElement(0, getNextFreeId());

        for (auto &vertID : eleVertices)
          newElement.getChildVertices().emplace_back(&getVertices()[vertID]);

        getElements().push_back(newElement);
      }

      // collect all elements pointers of each vertex
      for (auto &element : getElements())
        for (auto &vert : vertices(element))
          vert->getFatherElements().emplace_back(&element);

      // add edges to the grid
      if constexpr (SimpleGrid<dimension,dimensionworld>::dimension > 1) {
        for (auto &edge : edgesVertexIndices) {
          getSubEntities<dimension - 1>().emplace_back(0, getNextFreeId());
          auto &newEdge = getEdges().back();

          for (auto &&verticesIndicesOfEdge : edge) {
            // add vertex pointers to edge
            newEdge.getChildVertices().push_back(&getVertices()[verticesIndicesOfEdge]);
            // add edge pointers to vertices
            getVertices()[verticesIndicesOfEdge].template getFatherEntities<dimension - 1>().push_back(&newEdge);
          }
        }

        auto eIt = getElements().begin();
        // add edge pointers to elements
        for (auto &elementedgeIndexPerElement : elementEdgeIndices) {
          for (auto &elementedgeIndex : elementedgeIndexPerElement)
            eIt->template getChildEntities<1>().push_back(&getEdges()[elementedgeIndex]);
          ++eIt;
        }
      }
      // add surfaces to the grid
      if constexpr (SimpleGrid<dimension,dimensionworld>::dimension > 2) {
        for (auto &surf : surfaceVertexIndices) {
          getSubEntities<dimension - 2>().emplace_back(0, getNextFreeId());
          auto &newSurface = getSurfaces().back();

          for (auto &&verticesIndicesOfSurface : surf) {
            // add vertex pointers to surface
            newSurface.getChildVertices().push_back(&getVertices()[verticesIndicesOfSurface]);
            // add surface pointers to vertices
            getVertices()[verticesIndicesOfSurface].template getFatherEntities<dimension - 2>().push_back(
                &newSurface);
          }
        }

        auto eIt = getElements().begin();
        // add surface pointers to the elements
        for (auto &elementSurfaceIndexPerElement : elementSurfaceIndices) {
          for (auto &elementSurfaceIndex : elementSurfaceIndexPerElement)
            eIt->template getChildEntities<2>().push_back(&getSurfaces()[elementSurfaceIndex]);
          ++eIt;
        }
      }
    };

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

    //template<int dimension, int dimensionworld>
    //friend class SimpleGridFactory<dimension, dimensionworld>;
    friend class SimpleGridView<dim, dimworld, SimpleGrid>;

    std::vector<GridEntityTuple> gridEntities;

    unsigned int freeIdCounter{};
  };

}  // namespace Ikarus::Grid

#include <ikarus/Grids/SimpleGrid/SimpleGridFactory.h>
