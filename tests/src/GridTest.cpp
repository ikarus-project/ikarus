//
// Created by Alex on 25.05.2021.
//
#include <gmock/gmock.h>
#include <gtest/gtest.h>

#include "testHelpers.h"

#include <dune/geometry/type.hh>

#include <ikarus/Geometries/GeometryType.h>
#include <ikarus/Grids/GridHelper/griddrawer.h>
#include <ikarus/Grids/SimpleGrid/SimpleGrid.h>

/** @addtogroup Tests
 *  This module includes all tests
 *  @{
 */

/**
 * \addtogroup GridTests
 * \test This test checks the insertion of vertices and elements in a SimpleGrid<2,2>
 *
 * It also checks the correct obtain unique identifiers from getID().
 *
 * The tested grid looks as follows
 *        13           16
 *   2------------3------------5                \n
 *   |            |            |  \             \n
 *   |    7       |   8        |     \ 16       \n
 * 10|   quad   12|    quad  15|       \        \n
 *   |            |            |    9    \      \n
 *   |            |            |triangle  \     \n
 *   0------------1------------4------------6   \n
 *      11              14           17         \n
 */
TEST(GridTest, GridViewTest) {
  using namespace Ikarus::Grid;
  using Grid = SimpleGrid<2, 2>;
  SimpleGridFactory<2, 2> gridFactory;
  using vertexType = Eigen::Vector2d;
  std::vector<vertexType> verticesVec;
  verticesVec.emplace_back(0.0, 0.0);  // 0
  verticesVec.emplace_back(2.0, 0.0);  // 1
  verticesVec.emplace_back(0.0, 2.0);  // 2
  verticesVec.emplace_back(2.0, 2.0);  // 3
  verticesVec.emplace_back(4.0, 0.0);  // 4
  verticesVec.emplace_back(4.0, 2.0);  // 5
  verticesVec.emplace_back(6.0, 0.0);  // 6

  for (auto &&vert : verticesVec)
    gridFactory.insertVertex(vert);

  std::vector<size_t> elementIndices;
  elementIndices.resize(4);
  elementIndices = {0, 1, 2, 3};
  gridFactory.insertElement(Ikarus::GeometryType::linearQuadrilateral, elementIndices);
  elementIndices = {1, 4, 3, 5};
  gridFactory.insertElement(Ikarus::GeometryType::linearQuadrilateral, elementIndices);
  elementIndices.resize(3);
  elementIndices = {4, 6, 5};
  gridFactory.insertElement(Ikarus::GeometryType::linearTriangle, elementIndices);

  Grid grid = gridFactory.createGrid();

  auto gridView = grid.leafGridView();
  EXPECT_EQ(edges(gridView).size(), 9);
  EXPECT_EQ(surfaces(gridView).size(), 3);
  EXPECT_EQ(vertices(gridView).size(), 7);

  int expectedEdgeId = 10;
  std::vector<std::array<int, 2>> expectedEdgeVertexId;
  expectedEdgeVertexId.push_back({0, 2});
  expectedEdgeVertexId.push_back({0, 1});
  expectedEdgeVertexId.push_back({1, 3});
  expectedEdgeVertexId.push_back({2, 3});
  expectedEdgeVertexId.push_back({1, 4});
  expectedEdgeVertexId.push_back({4, 5});
  expectedEdgeVertexId.push_back({3, 5});
  expectedEdgeVertexId.push_back({4, 6});
  expectedEdgeVertexId.push_back({5, 6});

  int elementCounter = 0;
  for (auto &&edge : edges(gridView)) {
    { EXPECT_EQ(edge.getID(), expectedEdgeId++); }

    int vertexCounter = 0;
    for (auto &&vertex : vertices(edge)) {
      EXPECT_EQ(vertex->type(), Ikarus::GeometryType::vertex);
      EXPECT_EQ(vertex->getID(), expectedEdgeVertexId[elementCounter][vertexCounter]);
      ++vertexCounter;
    }
    ++elementCounter;
  }

  std::vector<std::vector<int>> expectedElementEdgeIds;
  expectedElementEdgeIds.push_back({10, 11, 12, 13});
  expectedElementEdgeIds.push_back({12, 14, 15, 16});
  expectedElementEdgeIds.push_back({17, 15, 18});

  int eleCounter = 0;
  for (auto &&singleElement : surfaces(gridView)) {
    int edgeCounter = 0;
    for (auto &&edge : edges(singleElement)) {
      EXPECT_EQ(edge->type(), Ikarus::GeometryType::linearLine);
      EXPECT_EQ(edge->getID(), expectedElementEdgeIds[eleCounter][edgeCounter]);
      ++edgeCounter;
    }
    ++eleCounter;
  }
  auto ele1 = surfaces(gridView).begin();
  EXPECT_THROW([[maybe_unused]] auto e = ele1->subEntities(3), std::logic_error);
  EXPECT_EQ(ele1->subEntities(2), 4);
  EXPECT_EQ(ele1->subEntities(1), 4);
  ++ele1;
  EXPECT_THROW([[maybe_unused]] auto e = ele1->subEntities(3), std::logic_error);
  EXPECT_EQ(ele1->subEntities(2), 4);
  EXPECT_EQ(ele1->subEntities(1), 4);
  ++ele1;
  EXPECT_THROW([[maybe_unused]] auto e = ele1->subEntities(3), std::logic_error);
  EXPECT_EQ(ele1->subEntities(2), 3);
  EXPECT_EQ(ele1->subEntities(1), 3);
}

/**
 * \addtogroup GridTests
 * \test This test checks the insertion of vertices and elements in a SimpleGrid<2,3>
 *
 */
TEST(GridTest, GridView3DSurfaceTest) {
  using namespace Ikarus::Grid;
  using Grid = SimpleGrid<2, 3>;
  SimpleGridFactory<2, 3> gridFactory;
  using vertexType = Eigen::Vector3d;
  std::vector<vertexType> verticesVec;
  verticesVec.emplace_back(0.0, 0.0, -3.0);  // 0
  verticesVec.emplace_back(2.0, 0.0, -3.0);  // 1
  verticesVec.emplace_back(0.0, 2.0, +3.0);  // 2
  verticesVec.emplace_back(2.0, 2.0, -3.0);  // 3
  verticesVec.emplace_back(4.0, 0.0, -3.0);  // 4
  verticesVec.emplace_back(4.0, 2.0, +3.0);  // 5
  verticesVec.emplace_back(6.0, 0.0, -3.0);  // 6

  for (auto &&vert : verticesVec)
    gridFactory.insertVertex(vert);

  std::vector<size_t> elementIndices;
  elementIndices.resize(4);
  elementIndices = {0, 1, 2, 3};
  gridFactory.insertElement(Ikarus::GeometryType::linearQuadrilateral, elementIndices);
  elementIndices = {1, 4, 3, 5};
  gridFactory.insertElement(Ikarus::GeometryType::linearQuadrilateral, elementIndices);
  elementIndices.resize(3);
  elementIndices = {4, 6, 5};
  gridFactory.insertElement(Ikarus::GeometryType::linearTriangle, elementIndices);

  Grid actualGrid = gridFactory.createGrid();
  auto gridView   = actualGrid.leafGridView();
  EXPECT_EQ(edges(gridView).size(), 9);
  EXPECT_EQ(surfaces(gridView).size(), 3);
  EXPECT_EQ(vertices(gridView).size(), 7);

  for (int i = 0; auto &&vertex : vertices(gridView)) {
    EXPECT_EQ(vertex.type(), Ikarus::GeometryType::vertex);
    EXPECT_EQ(vertex.getPosition(), verticesVec[i]);
    ++i;
  }

  auto &&eleIterator = surfaces(gridView).begin();
  EXPECT_EQ(eleIterator->type(), Ikarus::GeometryType::linearQuadrilateral);
  ++eleIterator;
  EXPECT_EQ(eleIterator->type(), Ikarus::GeometryType::linearQuadrilateral);
  ++eleIterator;
  EXPECT_EQ(eleIterator->type(), Ikarus::GeometryType::linearTriangle);

  std::vector<std::vector<int>> expectedElementEdgeIds;
  expectedElementEdgeIds.push_back({10, 11, 12, 13});
  expectedElementEdgeIds.push_back({12, 14, 15, 16});
  expectedElementEdgeIds.push_back({17, 15, 18});

  int eleCounter = 0;
  for (auto &singleElement : surfaces(gridView)) {
    int edgeCounter = 0;
    for (auto &&edge : edges(singleElement)) {
      EXPECT_EQ(edge->type(), Ikarus::GeometryType::linearLine);
      EXPECT_EQ(edge->getID(), expectedElementEdgeIds[eleCounter][edgeCounter]);
      ++edgeCounter;
    }
    ++eleCounter;
  }
  auto ele1 = surfaces(gridView).begin();
  EXPECT_THROW([[maybe_unused]] auto e = ele1->subEntities(3), std::logic_error);
  EXPECT_EQ(ele1->subEntities(2), 4);
  EXPECT_EQ(ele1->subEntities(1), 4);
  ++ele1;
  EXPECT_THROW([[maybe_unused]] auto e = ele1->subEntities(3), std::logic_error);
  EXPECT_EQ(ele1->subEntities(2), 4);
  EXPECT_EQ(ele1->subEntities(1), 4);
  ++ele1;
  EXPECT_THROW([[maybe_unused]] auto e = ele1->subEntities(3), std::logic_error);
  EXPECT_EQ(ele1->subEntities(2), 3);
  EXPECT_EQ(ele1->subEntities(1), 3);
}

/**
 * \addtogroup GridTests
 * \test This test checks the insertion of vertices and elements in a SimpleGrid<2,3>
 *
 */
TEST(GridTest, GridView3DSolidTest) {
  using namespace Ikarus::Grid;
  using Grid = SimpleGrid<3, 3>;
  SimpleGridFactory<3, 3> gridFactory;
  using vertexType = Eigen::Vector3d;
  std::vector<vertexType> verticesVec;
  verticesVec.emplace_back(0.0, 0.0, -3.0);  // 0
  verticesVec.emplace_back(2.0, 0.0, -3.0);  // 1
  verticesVec.emplace_back(0.0, 2.0, -3.0);  // 2
  verticesVec.emplace_back(2.0, 2.0, -3.0);  // 3
  verticesVec.emplace_back(0.0, 0.0, +3.0);  // 4
  verticesVec.emplace_back(2.0, 0.0, +3.0);  // 5
  verticesVec.emplace_back(0.0, 2.0, +3.0);  // 6
  verticesVec.emplace_back(2.0, 2.0, +3.0);  // 7
  verticesVec.emplace_back(4.0, 0.0, +3.0);  // 8

  for (auto &&vert : verticesVec)
    gridFactory.insertVertex(vert);

  std::vector<size_t> elementIndices;
  elementIndices.resize(8);
  elementIndices = {0, 1, 2, 3, 4, 5, 6, 7};
  gridFactory.insertElement(Ikarus::GeometryType::linearHexahedron, elementIndices);
  elementIndices.resize(4);
  elementIndices = {1, 8, 3, 5};
  gridFactory.insertElement(Ikarus::GeometryType::linearTetrahedron, elementIndices);

  Grid actualGrid = gridFactory.createGrid();

  auto gridView = actualGrid.leafGridView();
  // Element   Edge        VertexIDs
  std::vector<std::vector<std::array<int, 2>>> expectedElementEdgeVertexId;
  expectedElementEdgeVertexId.emplace_back();
  expectedElementEdgeVertexId[0].push_back({0, 4});  // 0
  expectedElementEdgeVertexId[0].push_back({1, 5});  // 1
  expectedElementEdgeVertexId[0].push_back({2, 6});  // 2
  expectedElementEdgeVertexId[0].push_back({3, 7});  // 3
  expectedElementEdgeVertexId[0].push_back({0, 2});  // 4
  expectedElementEdgeVertexId[0].push_back({1, 3});  // 5
  expectedElementEdgeVertexId[0].push_back({0, 1});  // 6
  expectedElementEdgeVertexId[0].push_back({2, 3});  // 7
  expectedElementEdgeVertexId[0].push_back({4, 6});  // 8
  expectedElementEdgeVertexId[0].push_back({5, 7});  // 9
  expectedElementEdgeVertexId[0].push_back({4, 5});  // 10
  expectedElementEdgeVertexId[0].push_back({6, 7});  // 11
  expectedElementEdgeVertexId.emplace_back();
  expectedElementEdgeVertexId[1].push_back({1, 8});  // 0
  expectedElementEdgeVertexId[1].push_back({1, 3});  // 1
  expectedElementEdgeVertexId[1].push_back({3, 8});  // 2
  expectedElementEdgeVertexId[1].push_back({1, 5});  // 3
  expectedElementEdgeVertexId[1].push_back({5, 8});  // 4
  expectedElementEdgeVertexId[1].push_back({3, 5});  // 4
  EXPECT_TRUE(!edges(gridView).empty());
  EXPECT_TRUE(!volumes(gridView).empty());

  for (int EleIter = 0; auto &&ele : volumes(gridView)) {
    EXPECT_TRUE(!edges(ele).empty());
    for (int edgeIter = 0; auto &&edge : edges(ele)) {
      EXPECT_TRUE(!vertices(edge).empty());
      for (int i = 0; auto &&verticesOfEdge : vertices(edge)) {
        EXPECT_EQ(verticesOfEdge->getID(), expectedElementEdgeVertexId[EleIter][edgeIter][i]);
        EXPECT_THAT(verticesOfEdge->getPosition(),
                    EigenApproxEqual(verticesVec[expectedElementEdgeVertexId[EleIter][edgeIter][i]], 1e-15));
        ++i;
      }
      ++edgeIter;
    }
    ++EleIter;
  }

  auto &&eleIterator = volumes(gridView).begin();
  EXPECT_EQ(eleIterator->type(), Ikarus::GeometryType::linearHexahedron);
  ++eleIterator;
  EXPECT_EQ(eleIterator->type(), Ikarus::GeometryType::linearTetrahedron);

  std::vector<int> expectedEdgesAtVertex{3, 4, 3, 5, 3, 5, 3, 3, 3};
  for (int i = 0; auto &&vertex : vertices(gridView))
    EXPECT_EQ(edges(vertex).size(), expectedEdgesAtVertex[i++]);

  auto ele1 = volumes(gridView).begin();
  EXPECT_EQ(ele1->subEntities(3), 8);
  EXPECT_EQ(ele1->subEntities(2), 12);
  EXPECT_EQ(ele1->subEntities(1), 6);

  // surface tests
  std::vector<std::vector<std::vector<int>>> expectedElementSurfaceVertexId;
  expectedElementSurfaceVertexId.emplace_back();
  expectedElementSurfaceVertexId[0].push_back({0, 2, 4, 6});  // 0
  expectedElementSurfaceVertexId[0].push_back({1, 3, 5, 7});  // 1
  expectedElementSurfaceVertexId[0].push_back({0, 1, 4, 5});  // 2
  expectedElementSurfaceVertexId[0].push_back({2, 3, 6, 7});  // 3
  expectedElementSurfaceVertexId[0].push_back({0, 1, 2, 3});  // 4
  expectedElementSurfaceVertexId[0].push_back({4, 5, 6, 7});  // 5

  expectedElementSurfaceVertexId.emplace_back();
  expectedElementSurfaceVertexId[1].push_back({1, 3, 8});  // 0
  expectedElementSurfaceVertexId[1].push_back({1, 5, 8});  // 1
  expectedElementSurfaceVertexId[1].push_back({1, 3, 5});  // 2
  expectedElementSurfaceVertexId[1].push_back({3, 5, 8});  // 3

  for (int EleIter = 0; auto &&ele : volumes(gridView)) {
    EXPECT_TRUE(!edges(ele).empty());
    for (int surfIter = 0; auto &&surf : surfaces(ele)) {
      EXPECT_TRUE(!vertices(surf).empty());
      for (int i = 0; auto &&verticesOfSurface : vertices(surf)) {
        EXPECT_EQ(verticesOfSurface->getID(), expectedElementSurfaceVertexId[EleIter][surfIter][i]);
        EXPECT_THAT(verticesOfSurface->getPosition(),
                    EigenApproxEqual(verticesVec[expectedElementSurfaceVertexId[EleIter][surfIter][i]], 1e-15));
        ++i;
      }
      ++surfIter;
    }
    ++EleIter;
  }
  std::vector<int> expectedSurfacesAtVertex{3, 6, 3, 6, 3, 6, 3, 3, 3};
  for (int i = 0; auto &&vertex : vertices(gridView))
    EXPECT_EQ(surfaces(vertex).size(), expectedSurfacesAtVertex[i++]);
}

TEST(GridTest, GridInsertionException) {
  using namespace Ikarus::Grid;
  SimpleGridFactory<2, 2> gridFactory;
  using vertexType = Eigen::Vector2d;
  std::vector<vertexType> verticesVec;
  verticesVec.emplace_back(0.0, 0.0);  // 0
  verticesVec.emplace_back(2.0, 0.0);  // 1
  verticesVec.emplace_back(0.0, 2.0);  // 2
  verticesVec.emplace_back(2.0, 2.0);  // 3

  for (auto &&vert : verticesVec)
    gridFactory.insertVertex(vert);

  std::vector<size_t> elementIndices;
  elementIndices.resize(4);
  elementIndices = {0, 1, 2, 3};
  EXPECT_THROW(gridFactory.insertElement(Ikarus::GeometryType::linearTriangle, elementIndices), Dune::GridError);

  EXPECT_THROW(gridFactory.insertElement(Ikarus::GeometryType::linearHexahedron, elementIndices), Dune::GridError);

  EXPECT_THROW(gridFactory.insertElement(Ikarus::GeometryType::linearLine, elementIndices), Dune::GridError);
  elementIndices.resize(3);
  elementIndices = {0, 1, 2};

  EXPECT_THROW(gridFactory.insertElement(Ikarus::GeometryType::linearQuadrilateral, elementIndices), Dune::GridError);
}

TEST(GridTest, GridEmptyGridCreation) {
  using namespace Ikarus::Grid;
  SimpleGridFactory<2, 2> gridFactory;
  EXPECT_THROW(gridFactory.createGrid(), Dune::GridError);
  gridFactory.insertVertex({2.0, 1.0});
  EXPECT_THROW(gridFactory.createGrid(), Dune::GridError);
}
/*\@}*/
