//
// Created by Alex on 25.05.2021.
//

#pragma once
#include <unordered_set>

#include <ikarus/Grids/GridInterface.h>
#include <ikarus/utils/LinearAlgebraTypedefs.h>
#include <ikarus/utils/std/algorithms.h>

namespace Ikarus::Grid {
  template <Concepts::Grid GridType>
  class SimpleGridFactory {
  private:
    /** \brief dimension of the grid */
    static const int dimension = GridType::dimension;

    /** \brief The grid world dimension */
    static const int dimensionworld = GridType::dimensionworld;

    /** \brief Type used by the grid for coordinates */
    using ctype = typename GridType::ctype;

    /** \brief Type used by the grid for the vertex coordinates */
    using VertexCoordinateType = FixedVector<double, dimensionworld>;

  public:
    void insertElement(const Dune::GeometryType& type, const DynArrayXi& vertices);
    void insertVertex(const VertexCoordinateType& pos);

    GridType createGrid();

  private:
    struct VertexIndexPair {
      FixedVector<double, dimensionworld> vertex;
      unsigned int index;
    };

    void insertVertexIndicesinEdge(std::array<int, 2>&& indices) {
      std::ranges::sort(indices);
      auto index = appendUnique(edgesVertexIndices, indices);
      elementEdgeIndices.back().push_back(index);
    }

    /** \brief Counter that creates the vertex indices */
    unsigned int vertexIndex{};
    std::vector<std::array<int, 2>> edgesVertexIndices;
    std::vector<VertexIndexPair> verticesPositions;
    std::vector<DynArrayXi> elementsVertices;
    std::vector<std::vector<size_t>> elementEdgeIndices;
  };

}  // namespace Ikarus::Grid
#include "SimpleGridFactory.inl"