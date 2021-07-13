//
// Created by Alex on 26.05.2021.
//
#pragma once
#include <dune/grid/common/exceptions.hh>

#include <ikarus/utils/std/algorithms.h>
namespace Ikarus::Grid {

  template <Concepts::Grid GridType>
  void SimpleGridFactory<GridType>::insertElement(const Dune::GeometryType &type, const std::span<size_t> verticesIn) {
    if (type.dim() != dimension) DUNE_THROW(Dune::GridError, "The inserted element has wrong dimensions!");

    storeVerticesIndicesOfEdges(type, verticesIn);
    if constexpr (dimension > 1) storeVerticesIndicesOfSurfaces(type, verticesIn);

    elementsVertices.emplace_back(verticesIn.begin(), verticesIn.end());
  }

  template <Concepts::Grid GridType>
  void SimpleGridFactory<GridType>::insertVertex(const VertexCoordinateType &pos) {
    verticesPositions.template emplace_back(SimpleGridFactory::VertexIndexPair{pos, vertexIndex++});
  }

  template <Concepts::Grid GridType>
  GridType SimpleGridFactory<GridType>::createGrid() {
    GridType grid;
    if (verticesPositions.empty())
      DUNE_THROW(Dune::GridError, "verticesPositions vector is empty. Unable to create Grid");
    if (elementsVertices.empty()) DUNE_THROW(Dune::GridError, "elements vector is empty. Unable to create Grid");

    grid.gridEntities.resize(1);  // resize to hold coarsest grid level

    grid.getVertices().reserve(verticesPositions.size());
    grid.getElements().reserve(elementsVertices.size());

    // add vertices to the grid
    for (auto &vert : verticesPositions) {
      typename GridType::VertexType newVertex(0, vert.vertex, grid.getNextFreeId());
      newVertex.levelIndex = vert.index;
      grid.getVertices().push_back(newVertex);
    }

    // add element and set vertex pointer of elements
    for (auto &eleVertices : elementsVertices) {
      typename GridType::RootEntity newElement(0, grid.getNextFreeId());

      for (auto &vertID : eleVertices)
        newElement.getChildVertices().emplace_back(&grid.getVertices()[vertID]);

      grid.getElements().push_back(newElement);
    }

    // collect all elements pointers of each vertex
    for (auto &element : grid.getElements())
      for (auto &vert : vertices(element))
        vert->getFatherElements().emplace_back(&element);

    // add edges to the grid
    if constexpr (GridType::dimension > 1) {
      for (auto &edge : edgesVertexIndices) {
        grid.template getSubEntities<dimension - 1>().emplace_back(0, grid.getNextFreeId());
        auto &newEdge = grid.getEdges().back();

        for (auto &&verticesIndicesOfEdge : edge) {
          // add vertex pointers to edge
          newEdge.getChildVertices().push_back(&grid.getVertices()[verticesIndicesOfEdge]);
          // add edge pointers to vertices
          grid.getVertices()[verticesIndicesOfEdge].template getFatherEntities<dimension - 1>().push_back(&newEdge);
        }
      }

      auto eIt = grid.getElements().begin();
      // add edge pointers to elements
      for (auto &elementedgeIndexPerElement : elementEdgeIndices) {
        for (auto &elementedgeIndex : elementedgeIndexPerElement)
          eIt->template getChildEntities<1>().push_back(&grid.getEdges()[elementedgeIndex]);
        ++eIt;
      }
    }
    // add surfaces to the grid
    if constexpr (GridType::dimension > 2) {
      for (auto &surf : surfaceVertexIndices) {
        grid.template getSubEntities<dimension - 2>().emplace_back(0, grid.getNextFreeId());
        auto &newSurface = grid.getSurfaces().back();

        for (auto &&verticesIndicesOfSurface : surf) {
          // add vertex pointers to surface
          newSurface.getChildVertices().push_back(&grid.getVertices()[verticesIndicesOfSurface]);
          // add surface pointers to vertices
          grid.getVertices()[verticesIndicesOfSurface].template getFatherEntities<dimension - 2>().push_back(
              &newSurface);
        }
      }

      auto eIt = grid.getElements().begin();
      // add surface pointers to the elements
      for (auto &elementSurfaceIndexPerElement : elementSurfaceIndices) {
        for (auto &elementSurfaceIndex : elementSurfaceIndexPerElement)
          eIt->template getChildEntities<2>().push_back(&grid.getSurfaces()[elementSurfaceIndex]);
        ++eIt;
      }
    }
    return grid;
  }

  /**
   * \brief This takes care of the storage of the vertx pointers for each edge
   * The numbering of vertices, edges, surfaces correspond to \cite sander2020dune Fig 5.12, 5.13
   *
   **/
  template <Concepts::Grid GridType>
  void SimpleGridFactory<GridType>::insertElement(const Dune::GeometryType &type, const std::span<size_t> verticesIn) {
    if (type.dim() != dimension) DUNE_THROW(Dune::GridError, "The inserted element has wrong dimensions!");

    storeVerticesIndicesOfEdges(type, verticesIn);
    if constexpr (dimension > 1) storeVerticesIndicesOfSurfaces(type, verticesIn);

    elementsVertices.emplace_back(verticesIn.begin(), verticesIn.end());
  }

  template <Concepts::Grid GridType>
  void SimpleGridFactory<GridType>::insertVertex(const VertexCoordinateType &pos) {
    verticesPositions.template emplace_back(SimpleGridFactory::VertexIndexPair{pos, vertexIndex++});
  }

  template <Concepts::Grid GridType>
  GridType SimpleGridFactory<GridType>::createGrid() {
    GridType grid;
    if (verticesPositions.empty())
      DUNE_THROW(Dune::GridError, "verticesPositions vector is empty. Unable to create Grid");
    if (elementsVertices.empty()) DUNE_THROW(Dune::GridError, "elements vector is empty. Unable to create Grid");

    grid.gridEntities.resize(1);  // resize to hold coarsest grid level

    grid.getVertices().reserve(verticesPositions.size());
    grid.getElements().reserve(elementsVertices.size());

    // add vertices to the grid
    for (auto &vert : verticesPositions) {
      typename GridType::VertexType newVertex(0, vert.vertex, grid.getNextFreeId());
      newVertex.levelIndex = vert.index;
      grid.getVertices().push_back(newVertex);
    }

    // add element and set vertex pointer of elements
    for (auto &eleVertices : elementsVertices) {
      typename GridType::RootEntity newElement(0, grid.getNextFreeId());

      for (auto &vertID : eleVertices)
        newElement.getChildVertices().emplace_back(&grid.getVertices()[vertID]);

      grid.getElements().push_back(newElement);
    }

    // set elements pointers of vertex
    for (auto &element : grid.getElements())
      for (auto &vert : vertices(element))
        vert->getFatherElements().emplace_back(&element);

    // add edges to the grid
    if constexpr (GridType::dimension > 1) {
      for (auto &edge : edgesVertexIndices) {
        grid.template getSubEntities<dimension - 1>().emplace_back(0, grid.getNextFreeId());
        auto &newEdge = grid.getEdges().back();

        for (auto &&verticesIndicesOfEdge : edge) {
          newEdge.getChildVertices().push_back(&grid.getVertices()[verticesIndicesOfEdge]);
          grid.getVertices()[verticesIndicesOfEdge].template getFatherEntities<dimension - 1>().push_back(&newEdge);
        }
      }

      auto eIt = grid.getElements().begin();
      // add edges pointers to the elements
      for (auto &elementedgeIndexPerElement : elementEdgeIndices) {
        for (auto &elementedgeIndex : elementedgeIndexPerElement)
          eIt->template getChildEntities<1>().push_back(&grid.getEdges()[elementedgeIndex]);
        ++eIt;
      }
    }
    if constexpr (GridType::dimension > 2) {
      for (auto &surf : surfaceVertexIndices) {
        grid.template getSubEntities<dimension - 2>().emplace_back(0, grid.getNextFreeId());
        auto &newSurface = grid.getSurfaces().back();

        for (auto &&verticesIndicesOfSurface : surf) {
          newSurface.getChildVertices().push_back(&grid.getVertices()[verticesIndicesOfSurface]);
          grid.getVertices()[verticesIndicesOfSurface].template getFatherEntities<dimension - 2>().push_back(
              &newSurface);
        }
      }

      auto eIt = grid.getElements().begin();
      // add surface pointers to the elements
      for (auto &elementSurfaceIndexPerElement : elementSurfaceIndices) {
        for (auto &elementSurfaceIndex : elementSurfaceIndexPerElement)
          eIt->template getChildEntities<2>().push_back(&grid.getSurfaces()[elementSurfaceIndex]);
        ++eIt;
      }
    }
    return grid;
  }

  template <Concepts::Grid GridType>
  void SimpleGridFactory<GridType>::storeVerticesIndicesOfEdges(const Dune::GeometryType &type,
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
  template <Concepts::Grid GridType>
  void SimpleGridFactory<GridType>::storeVerticesIndicesOfSurfaces(const Dune::GeometryType &type,
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

  template <Concepts::Grid GridType>
  void SimpleGridFactory<GridType>::insertVertexIndicesinEdge(std::vector<size_t>&& indices) {
    std::ranges::sort(indices);
    auto index = Ikarus::stl::appendUnique(edgesVertexIndices, indices);
    elementEdgeIndices.back().push_back(index);
  }

  template <Concepts::Grid GridType>
  void SimpleGridFactory<GridType>::insertVertexIndicesinSurface(std::vector<size_t>&& indices) {
    std::ranges::sort(indices);
    auto index = Ikarus::stl::appendUnique(surfaceVertexIndices, indices);
    elementSurfaceIndices.back().push_back(index);
  }
}  // namespace Ikarus::Grid