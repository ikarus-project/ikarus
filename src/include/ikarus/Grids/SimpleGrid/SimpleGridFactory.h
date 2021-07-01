//
// Created by Alex on 25.05.2021.
//

#pragma once
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
    void insertElement(const Dune::GeometryType& type, const std::span<size_t> vertices);
    void insertVertex(const VertexCoordinateType& pos);

    GridType createGrid();

  private:
    struct VertexIndexPair {
      FixedVector<double, dimensionworld> vertex;
      size_t index;
    };

    void insertVertexIndicesinEdge(std::vector<size_t>&& indices) {
      std::ranges::sort(indices);
      auto index = Ikarus::stl::appendUnique(edgesVertexIndices, indices);
      elementEdgeIndices.back().push_back(index);
    }

    /** \brief Counter that creates the vertex indices */
    size_t vertexIndex{};
    std::vector<std::vector<size_t>> edgesVertexIndices;
    std::vector<VertexIndexPair> verticesPositions;
    std::vector<std::vector<size_t>> elementsVertices;
    std::vector<std::vector<size_t>> elementEdgeIndices;
  };

}  // namespace Ikarus::Grid
#include "SimpleGridFactory.inl"