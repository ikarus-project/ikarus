//
// Created by Alex on 26.05.2021.
//
#pragma once
#include <dune/grid/common/exceptions.hh>

#include <ikarus/utils/std/algorithms.h>
namespace Ikarus::Grid {

  /**
   * \brief This function inserts the index of the element and creates if the dimension of the grid fits
   *  the corresponding edges and surface. The numbering of vertices, edges, surfaces correspond to \cite sander2020dune
   *Fig 5.12, 5.13
   *
   * \param[in] type The geometric type of the element
   * \param[in] verticesIn The indices of the vertices of the element
   *
   **/
  template <Concepts::Grid GridType>
  void SimpleGridFactory<GridType>::insertElement(const Dune::GeometryType &type,const std::span<size_t> verticesIn) {
    if (type.dim() != dimension) DUNE_THROW(Dune::GridError, "The inserted element has wrong dimensions!");

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

      elementsVertices.emplace_back(verticesIn.begin(),verticesIn.end());
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
      newVertex.leafIndex  = vert.index;
      newVertex.levelIndex = vert.index;
      grid.getVertices().push_back(newVertex);
    }

    // add element and set vertex pointer of elements
    for (auto &eleVertices : elementsVertices) {
      typename GridType::ElementType newElement(0, grid.getNextFreeId());

      for (auto &vertID : eleVertices)
        newElement.getChildVertices().emplace_back(&grid.getVertices()[vertID]);

      grid.getElements().push_back(newElement);
    }

     //set elements pointers of vertex
        for (auto &element : grid.getElements())
          for (auto &vert : vertices(element))
            vert->getFatherElements().emplace_back(&element);

    // add edges to the grid
    if constexpr (GridType::dimension > 1) {
      for (auto &edge : edgesVertexIndices) {
        grid.template getSubEntities<dimension - 1>().emplace_back(0, grid.getNextFreeId());
        auto &newEdge = grid.getEdges().back();

        newEdge.getChildVertices().push_back(&grid.getVertices()[edge[0]]);
        newEdge.getChildVertices().push_back(&grid.getVertices()[edge[1]]);
      }

      auto eIt = grid.getElements().begin();
      // add edges pointers to the elements
      for (auto &elementedgeIndexPerElement : elementEdgeIndices) {
        for (auto &elementedgeIndex : elementedgeIndexPerElement)
          eIt->template getChildEntities<1>().push_back(&grid.getEdges()[elementedgeIndex]);
        ++eIt;
      }

       //add edge pointers to vertices
      // set elements pointers of vertex

            for (auto &&vert : grid.getVertices()) {
              auto hasVertex = [&vert](auto &edge) {
                return (std::ranges::find_if(edge.getChildVertices(),
                                             [&vert](auto &vertex) { return vert.getID() == vertex->getID(); })
                        != end(edge.getChildVertices()));
              };

              for (auto &edgeWhichHasTheVertex : std::ranges::filter_view(grid.getEdges(), hasVertex)) {
                vert.template getFatherEntities<dimension - 1>().push_back(&edgeWhichHasTheVertex);
              }
            }
    }

    return grid;
  }

}  // namespace Ikarus::Grid