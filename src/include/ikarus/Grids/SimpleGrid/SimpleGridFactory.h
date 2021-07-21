//
// Created by Alex on 25.05.2021.
//

#pragma once
#include <ikarus/Geometries/GeometryType.h>
#include <ikarus/utils/LinearAlgebraTypedefs.h>
#include <ikarus/utils/std/algorithms.h>

namespace Ikarus::Grid {
  template <int dimension, int dimensionworld>
  class SimpleGridFactory {
  private:
    /** \brief Type used by the grid for the vertex coordinates */
    using VertexCoordinateType = Eigen::Vector<double, dimensionworld>;

  public:
    void insertElement(Ikarus::GeometryType type, std::span<size_t> vertices);
    void insertVertex(const VertexCoordinateType& pos);

    auto createGrid();

  private:
    struct VertexIndexPair {
      Eigen::Vector<double, dimensionworld> vertex;
      size_t index;
    };
    void insertVertexIndicesinEdge(std::vector<size_t>&& indices);
    void insertVertexIndicesinSurface(std::vector<size_t>&& indices);

    void storeVerticesIndicesOfEdges(Ikarus::GeometryType type, std::span<size_t> verticesIn);
    void storeVerticesIndicesOfSurfaces(Ikarus::GeometryType type, std::span<size_t> verticesIn);

    std::vector<VertexIndexPair> verticesPositions;
    std::vector<std::vector<size_t>> edgesVertexIndices;
    std::vector<std::vector<size_t>> surfaceVertexIndices;
    std::vector<std::vector<size_t>> elementsVertices;
    std::vector<std::vector<size_t>> elementEdgeIndices;
    std::vector<std::vector<size_t>> elementSurfaceIndices;
    size_t vertexIndex{};
  };

}  // namespace Ikarus::Grid
#include "SimpleGridFactory.inl"