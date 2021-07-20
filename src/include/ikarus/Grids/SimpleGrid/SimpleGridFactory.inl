//
// Created by Alex on 26.05.2021.
//
#pragma once
#include <dune/grid/common/exceptions.hh>

#include <ikarus/utils/std/algorithms.h>
#include "SimpleGridTypedefs.h"
namespace Ikarus::Grid {
  template <int dimension, int dimensionworld>
  using GridType = SimpleGrid<dimension,dimensionworld>;

  template <int dimension, int dimensionworld>
  void SimpleGridFactory<dimension,dimensionworld>::insertElement(const Dune::GeometryType &type, const std::span<size_t> verticesIn) {
    if (type.dim() != dimension) DUNE_THROW(Dune::GridError, "The inserted element has wrong dimensions!");

    storeVerticesIndicesOfEdges(type, verticesIn);
    if constexpr (dimension > 1) storeVerticesIndicesOfSurfaces(type, verticesIn);

    elementsVertices.emplace_back(verticesIn.begin(), verticesIn.end());
  }

  template <int dimension, int dimensionworld>
  void SimpleGridFactory<dimension,dimensionworld>::insertVertex(const VertexCoordinateType &pos) {
    verticesPositions.template emplace_back(SimpleGridTypedefs::VertexIndexPair{pos, vertexIndex++});
  }

  template <int dimension, int dimensionworld>
  SimpleGrid<dimension,dimensionworld> SimpleGridFactory<dimension,dimensionworld>::createGrid() {
    return SimpleGrid<dimension,dimensionworld>(verticesPositions,edgesVertexIndices,surfaceVertexIndices,
                                               elementsVertices,elementEdgeIndices,elementSurfaceIndices);
  }

  /**
   * \brief This takes care of the storage of the vertx pointers for each edge
   * The numbering of vertices, edges, surfaces correspond to \cite sander2020dune Fig 5.12, 5.13
   *
   **/
  template <int dimension, int dimensionworld>
  void SimpleGridFactory<dimension,dimensionworld>::storeVerticesIndicesOfEdges(const Dune::GeometryType &type,
                                                                const std::span<size_t> verticesIn) {
    elementEdgeIndices.emplace_back();
    if (type.isLine()) {
      if (verticesIn.size() != 2)
        DUNE_THROW(Dune::GridError, "You have requested to enter a line, but you"
                                        << " have provided " << verticesIn.size() << " vertices!");
    } else if (type.isTriangle()) {
      if (verticesIn.size() != 3)
        DUNE_THROW(Dune::GridError, "You have requested to enter a triangle, but you"
                                        << " have provided " << verticesIn.size() << " vertices!");
      insertVertexIndicesinEdge({verticesIn[0], verticesIn[1]});  // 0
      insertVertexIndicesinEdge({verticesIn[0], verticesIn[2]});  // 1
      insertVertexIndicesinEdge({verticesIn[1], verticesIn[2]});  // 2
    } else if (type.isQuadrilateral()) {
      if (verticesIn.size() != 4)
        DUNE_THROW(Dune::GridError, "You have requested to enter a quadrilateral, but you"
                                        << " have provided " << verticesIn.size() << " vertices!");
      insertVertexIndicesinEdge({verticesIn[2], verticesIn[0]});  // 0
      insertVertexIndicesinEdge({verticesIn[0], verticesIn[1]});  // 1
      insertVertexIndicesinEdge({verticesIn[1], verticesIn[3]});  // 2
      insertVertexIndicesinEdge({verticesIn[3], verticesIn[2]});  // 3

    } else if (type.isTetrahedron()) {
      if (verticesIn.size() != 4)
        DUNE_THROW(Dune::GridError, "You have requested to enter a tetrahedron, but you"
                                        << " have provided " << verticesIn.size() << " vertices!");
      insertVertexIndicesinEdge({verticesIn[0], verticesIn[1]});  // 0
      insertVertexIndicesinEdge({verticesIn[2], verticesIn[0]});  // 1
      insertVertexIndicesinEdge({verticesIn[1], verticesIn[2]});  // 2
      insertVertexIndicesinEdge({verticesIn[3], verticesIn[0]});  // 3
      insertVertexIndicesinEdge({verticesIn[1], verticesIn[3]});  // 4
      insertVertexIndicesinEdge({verticesIn[2], verticesIn[3]});  // 5

    } else if (type.isPyramid()) {
      if (verticesIn.size() != 5)
        DUNE_THROW(Dune::GridError, "You have requested to enter a pyramid, but you"
                                        << " have provided " << verticesIn.size() << " vertices!");
      insertVertexIndicesinEdge({verticesIn[2], verticesIn[0]});  // 0
      insertVertexIndicesinEdge({verticesIn[1], verticesIn[3]});  // 1
      insertVertexIndicesinEdge({verticesIn[0], verticesIn[1]});  // 2
      insertVertexIndicesinEdge({verticesIn[3], verticesIn[2]});  // 3
      insertVertexIndicesinEdge({verticesIn[4], verticesIn[0]});  // 4
      insertVertexIndicesinEdge({verticesIn[1], verticesIn[4]});  // 5
      insertVertexIndicesinEdge({verticesIn[2], verticesIn[4]});  // 6
      insertVertexIndicesinEdge({verticesIn[3], verticesIn[4]});  // 7

    } else if (type.isPrism()) {
      if (verticesIn.size() != 6)
        DUNE_THROW(Dune::GridError, "You have requested to enter a prism, but you"
                                        << " have provided " << verticesIn.size() << " vertices!");
      insertVertexIndicesinEdge({verticesIn[3], verticesIn[0]});  // 0
      insertVertexIndicesinEdge({verticesIn[1], verticesIn[4]});  // 1
      insertVertexIndicesinEdge({verticesIn[2], verticesIn[5]});  // 2
      insertVertexIndicesinEdge({verticesIn[0], verticesIn[1]});  // 3
      insertVertexIndicesinEdge({verticesIn[2], verticesIn[0]});  // 4
      insertVertexIndicesinEdge({verticesIn[1], verticesIn[2]});  // 5
      insertVertexIndicesinEdge({verticesIn[3], verticesIn[4]});  // 6
      insertVertexIndicesinEdge({verticesIn[5], verticesIn[3]});  // 7
      insertVertexIndicesinEdge({verticesIn[4], verticesIn[5]});  // 8

    } else if (type.isHexahedron()) {
      if (verticesIn.size() != 8)
        DUNE_THROW(Dune::GridError, "You have requested to enter a hexahedron, but you"
                                        << " have provided " << verticesIn.size() << " vertices!");
      insertVertexIndicesinEdge({verticesIn[4], verticesIn[0]});  // 0
      insertVertexIndicesinEdge({verticesIn[1], verticesIn[5]});  // 1
      insertVertexIndicesinEdge({verticesIn[2], verticesIn[6]});  // 2
      insertVertexIndicesinEdge({verticesIn[3], verticesIn[7]});  // 3
      insertVertexIndicesinEdge({verticesIn[2], verticesIn[0]});  // 4
      insertVertexIndicesinEdge({verticesIn[1], verticesIn[3]});  // 5
      insertVertexIndicesinEdge({verticesIn[0], verticesIn[1]});  // 6
      insertVertexIndicesinEdge({verticesIn[3], verticesIn[2]});  // 7
      insertVertexIndicesinEdge({verticesIn[4], verticesIn[6]});  // 8
      insertVertexIndicesinEdge({verticesIn[5], verticesIn[7]});  // 9
      insertVertexIndicesinEdge({verticesIn[4], verticesIn[5]});  // 10
      insertVertexIndicesinEdge({verticesIn[7], verticesIn[6]});  // 11

    } else {
      DUNE_THROW(Dune::GridError, "You cannot insert a " << type << " into a SimpleGrid<" << dimensionworld << ">!");
    }
  }

  /**
   * \brief This takes care of the storage of the vertex pointers for each surface
   * The numbering of vertices, edges, surfaces correspond to \cite sander2020dune Fig 5.12, 5.13
   *
   **/
  template <int dimension, int dimensionworld>
  void SimpleGridFactory<dimension,dimensionworld>::storeVerticesIndicesOfSurfaces(const Dune::GeometryType &type,
                                                                         std::span<size_t> verticesIn) {
    elementSurfaceIndices.emplace_back();
    if (type.isLine()) {
      DUNE_THROW(Dune::GridError, "A line does not have surfaces.");
    } else if (type.isTriangle()) {
      if (verticesIn.size() != 3) DUNE_THROW(Dune::GridError, "A triangle does not have surfaces.");
    } else if (type.isQuadrilateral()) {
      if (verticesIn.size() != 4) DUNE_THROW(Dune::GridError, "A qudarilateral does not have surfaces.");

    } else if (type.isTetrahedron()) {
      if (verticesIn.size() != 4)
        DUNE_THROW(Dune::GridError, "You have requested to enter a tetrahedron, but you"
                                        << " have provided " << verticesIn.size() << " vertices!");
      insertVertexIndicesinSurface({verticesIn[0], verticesIn[1], verticesIn[2]});  // 0
      insertVertexIndicesinSurface({verticesIn[0], verticesIn[1], verticesIn[3]});  // 1
      insertVertexIndicesinSurface({verticesIn[0], verticesIn[2], verticesIn[3]});  // 2
      insertVertexIndicesinSurface({verticesIn[1], verticesIn[2], verticesIn[3]});  // 3

    } else if (type.isPyramid()) {
      if (verticesIn.size() != 5)
        DUNE_THROW(Dune::GridError, "You have requested to enter a pyramid, but you"
                                        << " have provided " << verticesIn.size() << " vertices!");
      insertVertexIndicesinSurface({verticesIn[0], verticesIn[1], verticesIn[2], verticesIn[3]});  // 0
      insertVertexIndicesinSurface({verticesIn[0], verticesIn[2], verticesIn[4]});                 // 1
      insertVertexIndicesinSurface({verticesIn[1], verticesIn[3], verticesIn[4]});                 // 2
      insertVertexIndicesinSurface({verticesIn[0], verticesIn[1], verticesIn[4]});                 // 3
      insertVertexIndicesinSurface({verticesIn[2], verticesIn[3], verticesIn[4]});                 // 4

    } else if (type.isPrism()) {
      if (verticesIn.size() != 6)
        DUNE_THROW(Dune::GridError, "You have requested to enter a prism, but you"
                                        << " have provided " << verticesIn.size() << " vertices!");
      insertVertexIndicesinSurface({verticesIn[0], verticesIn[1], verticesIn[3], verticesIn[4]});  // 0
      insertVertexIndicesinSurface({verticesIn[0], verticesIn[2], verticesIn[3], verticesIn[5]});  // 1
      insertVertexIndicesinSurface({verticesIn[1], verticesIn[2], verticesIn[4], verticesIn[5]});  // 2
      insertVertexIndicesinSurface({verticesIn[0], verticesIn[1], verticesIn[2]});                 // 3
      insertVertexIndicesinSurface({verticesIn[3], verticesIn[4], verticesIn[5]});                 // 4

    } else if (type.isHexahedron()) {
      if (verticesIn.size() != 8)
        DUNE_THROW(Dune::GridError, "You have requested to enter a hexahedron, but you"
                                        << " have provided " << verticesIn.size() << " vertices!");
      insertVertexIndicesinSurface({verticesIn[0], verticesIn[2], verticesIn[4], verticesIn[6]});  // 0
      insertVertexIndicesinSurface({verticesIn[1], verticesIn[3], verticesIn[5], verticesIn[7]});  // 1
      insertVertexIndicesinSurface({verticesIn[0], verticesIn[1], verticesIn[4], verticesIn[5]});  // 2
      insertVertexIndicesinSurface({verticesIn[2], verticesIn[3], verticesIn[6], verticesIn[7]});  // 3
      insertVertexIndicesinSurface({verticesIn[0], verticesIn[1], verticesIn[2], verticesIn[3]});  // 4
      insertVertexIndicesinSurface({verticesIn[4], verticesIn[5], verticesIn[6], verticesIn[7]});  // 5

    } else {
      DUNE_THROW(Dune::GridError, "You cannot insert a " << type << " into a SimpleGrid<" << dimensionworld << ">!");
    }
  }

  template <int dimension, int dimensionworld>
  void SimpleGridFactory<dimension,dimensionworld>::insertVertexIndicesinEdge(std::vector<size_t>&& indices) {
    std::ranges::sort(indices);
    auto index = Ikarus::stl::appendUnique(edgesVertexIndices, indices);
    elementEdgeIndices.back().push_back(index);
  }

  template <int dimension, int dimensionworld>
  void SimpleGridFactory<dimension,dimensionworld>::insertVertexIndicesinSurface(std::vector<size_t>&& indices) {
    std::ranges::sort(indices);
    auto index = Ikarus::stl::appendUnique(surfaceVertexIndices, indices);
    elementSurfaceIndices.back().push_back(index);
  }
}  // namespace Ikarus::Grid