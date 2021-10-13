//
// Created by Alex on 21.04.2021.
//

#include <gmock/gmock.h>
#include <gtest/gtest.h>

#include <fstream>
#include <vector>

#include <dune/geometry/type.hh>

#include <Eigen/Core>

#include <ikarus/FEManager/DefaultFEManager.h>
#include <ikarus/FiniteElements/ElasticityFE.h>
#include <ikarus/Geometries/GeometryType.h>
#include <ikarus/Grids/GridData.h>
#include <ikarus/Grids/SimpleGrid/SimpleGrid.h>
#include <ikarus/Variables/VariableVector.h>

class TestFE {
public:
  static void initialize() {}
  [[nodiscard]] static Ikarus::FiniteElements::IFiniteElement::DofPairVectorType getEntityVariablePairs() {
    return Ikarus::FiniteElements::IFiniteElement::DofPairVectorType{};
  }

  static double calculateScalar(const Ikarus::FiniteElements::ScalarAffordances&,
                                std::vector<Ikarus::Variable::IVariable*>&) {
    return 5;
  }
};

TEST(GridDataInterfaceTest, createDataOnEntities) {
  using namespace Ikarus;
  using namespace Ikarus::Grid;
  using namespace Ikarus::FiniteElements;
  using namespace Ikarus::Variable;

  using Grid = SimpleGrid<2, 2>;
  SimpleGridFactory<2, 2> gridFactory;
  using vertexType = Eigen::Vector2d;
  std::vector<vertexType> verticesVec;
  verticesVec.emplace_back(vertexType{0.0, 0.0});  // 0
  verticesVec.emplace_back(vertexType{2.0, 0.0});  // 1
  verticesVec.emplace_back(vertexType{0.0, 2.0});  // 2
  verticesVec.emplace_back(vertexType{2.0, 2.0});  // 3
  verticesVec.emplace_back(vertexType{4.0, 0.0});  // 4
  verticesVec.emplace_back(vertexType{4.0, 2.0});  // 5

  for (auto&& vert : verticesVec)
    gridFactory.insertVertex(vert);

  std::vector<size_t> elementIndices;
  elementIndices.resize(4);
  elementIndices = {0, 1, 2, 3};
  gridFactory.insertElement(Ikarus::GeometryType::linearQuadrilateral, elementIndices);
  elementIndices = {1, 4, 3, 5};
  gridFactory.insertElement(Ikarus::GeometryType::linearQuadrilateral, elementIndices);

  Grid grid = gridFactory.createGrid();

  auto gridView = grid.leafGridView();
  std::vector<Ikarus::FiniteElements::IFiniteElement> fes;

  for (auto&& element : surfaces(gridView))
    fes.emplace_back(Ikarus::FiniteElements::ElasticityFE(element, gridView.indexSet()));

  GridData gridData(gridView.indexSet());

  gridData.add(data(velocity2d), at(Vertices));
  gridData.add(data(pressure), at(RootEntities));
  gridData.add(data(displacement2d), at(EntitiesWithCoDim<2>()));
  gridData.add(data(edgeLength), at(EntitiesWithCoDim<1>()));

  std::array<VariablesTags, 2> expectedVarTags({velocity2d, displacement2d});
  for (auto& vertex : vertices(gridView)) {
    auto& vertexData = gridData.getData(vertex);
    for (int i = 0; auto& vertexVar : vertexData)
      EXPECT_THAT(getTag(vertexVar), static_cast<int>(expectedVarTags[i++]));
  }

  for (auto& surface : surfaces(gridView)) {
    auto& surfaceData = gridData.getData(surface);
    for (auto& surfaceVar : surfaceData)
      EXPECT_THAT(getTag(surfaceVar), static_cast<int>(pressure));
  }

  for (auto& edge : edges(gridView)) {
    auto& edgeData = gridData.getData(edge);
    for (auto& edgeVar : edgeData)
      EXPECT_THAT(getTag(edgeVar), static_cast<int>(edgeLength));
  }

  gridData.remove(data(velocity2d), at(Vertices));

  for (auto& vertex : vertices(gridView)) {
    auto& vertexData = gridData.getData(vertex);
    for (auto& vertexVar : vertexData)
      EXPECT_THAT(getTag(vertexVar), static_cast<int>(displacement2d));
  }
}

TEST(GridDataInterfaceTest, SimpleIndexSetTest) {
  {
    using namespace Ikarus::Grid;
    using namespace Ikarus::FiniteElements;
    using namespace Ikarus::Variable;

    using Grid = SimpleGrid<2, 2>;
    SimpleGridFactory<2, 2> gridFactory;
    using vertexType = Eigen::Vector2d;
    std::vector<vertexType> verticesVec;
    verticesVec.emplace_back(vertexType{0.0, 0.0});  // 0
    verticesVec.emplace_back(vertexType{2.0, 0.0});  // 1
    verticesVec.emplace_back(vertexType{0.0, 2.0});  // 2
    verticesVec.emplace_back(vertexType{2.0, 2.0});  // 3
    verticesVec.emplace_back(vertexType{4.0, 0.0});  // 4
    verticesVec.emplace_back(vertexType{4.0, 2.0});  // 5

    for (auto&& vert : verticesVec)
      gridFactory.insertVertex(vert);

    std::vector<size_t> elementIndices;
    elementIndices.resize(4);
    elementIndices = {0, 1, 2, 3};
    gridFactory.insertElement(Ikarus::GeometryType::linearQuadrilateral, elementIndices);
    elementIndices = {1, 4, 3, 5};
    gridFactory.insertElement(Ikarus::GeometryType::linearQuadrilateral, elementIndices);

    Grid grid = gridFactory.createGrid();

    auto gridView = grid.leafGridView();

    auto indexSet = gridView.indexSet();

    std::array<std::array<int, 4>, 2> expectedVertexIndices{{{0, 1, 2, 3}, {1, 4, 3, 5}}};
    std::array<std::array<int, 4>, 2> expectedEdgeIndices{{{0, 1, 2, 3}, {2, 4, 5, 6}}};
    for (int surfIndex = 0; auto&& surf : surfaces(gridView)) {
      for (size_t i = 0; i < surf.subEntities(2); ++i)
        EXPECT_THAT(expectedVertexIndices[surfIndex][i], indexSet.subIndex(surf, i, 2));
      for (size_t i = 0; i < surf.subEntities(1); ++i)
        EXPECT_THAT(expectedEdgeIndices[surfIndex][i], indexSet.subIndex(surf, i, 1));
      EXPECT_THROW([[maybe_unused]] auto i = surf.subEntities(0), std::logic_error);
      ++surfIndex;
    }
  }

  {
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

    for (auto&& vert : verticesVec)
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

    auto indexSet = gridView.indexSet();

    std::array<std::vector<int>, 3> expectedVertexIndices{{{0, 1, 2, 3}, {1, 4, 3, 5}, {4, 6, 5}}};
    std::array<std::vector<int>, 3> expectedEdgeIndices{{{0, 1, 2, 3}, {2, 4, 5, 6}, {7, 5, 8}}};
    for (int surfIndex = 0; auto&& surf : surfaces(gridView)) {
      for (size_t i = 0; i < surf.subEntities(2); ++i)
        EXPECT_THAT(expectedVertexIndices[surfIndex][i], indexSet.subIndex(surf, i, 2));
      for (size_t i = 0; i < surf.subEntities(1); ++i)
        EXPECT_THAT(expectedEdgeIndices[surfIndex][i], indexSet.subIndex(surf, i, 1));
      EXPECT_THROW([[maybe_unused]] auto i = surf.subEntities(0), std::logic_error);
      ++surfIndex;
    }
  }

  {
    using namespace Ikarus::Grid;
    using namespace Ikarus;
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

    for (auto&& vert : verticesVec)
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

    auto indexSet = gridView.indexSet();

    std::array<std::vector<int>, 2> expectedVertexIndices{{{0, 1, 2, 3, 4, 5, 6, 7}, {1, 8, 3, 5}}};
    std::array<std::vector<int>, 2> expectedEdgeIndices{
        {{0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11}, {12, 5, 13, 1, 14, 15}}};
    std::array<std::vector<int>, 2> expectedSurfaceIndices{{{0, 1, 2, 3, 4, 5}, {6, 7, 8, 9}}};
    for (int eleIndex = 0; auto&& ele : rootEntities(gridView)) {
      for (size_t i = 0; i < ele.subEntities(3); ++i)
        EXPECT_THAT(expectedVertexIndices[eleIndex][i], indexSet.subIndex(ele, i, 3));
      for (size_t i = 0; i < ele.subEntities(2); ++i)
        EXPECT_THAT(expectedEdgeIndices[eleIndex][i], indexSet.subIndex(ele, i, 2));
      for (size_t i = 0; i < ele.subEntities(1); ++i)
        EXPECT_THAT(expectedSurfaceIndices[eleIndex][i], indexSet.subIndex(ele, i, 1));
      EXPECT_THROW([[maybe_unused]] auto i = ele.subEntities(0), std::logic_error);
      ++eleIndex;
    }

    GridData gridData(indexSet);
    using namespace Ikarus::Variable;
    gridData.add(data(velocity1d), at(Vertices));
    gridData.add(data(pressure, datatuple(edgeLength, displacement1d)), at(Surfaces, Edges));

    std::array<std::vector<VariablesTags>, 2> surfaceExpectedSubVars(
        {{
             velocity1d, velocity1d, velocity1d, velocity1d,  //  4 vertices for each surface of the linearHexahedron
             edgeLength,                                      // 4 edges with 2 vars each
             displacement1d, edgeLength, displacement1d, edgeLength, displacement1d, edgeLength, displacement1d,
             pressure  // the var on the surface itself
         },
         {
             velocity1d,  //  3 vertices for each surface of the lineartetrahedron
             velocity1d, velocity1d,
             edgeLength,  // 3 edges with 2 vars each
             displacement1d, edgeLength, displacement1d, edgeLength, displacement1d,
             pressure  // the var on the surface itself
         }});
    for (int eleIndex = 0; auto&& ele : rootEntities(gridView)) {
      for (auto&& surface : surfaces(ele)) {
        auto surfaceSubData = gridData.getAllSubEntityData(surface);
        EXPECT_THAT(gridData.getData(surface)[0], VariableFactory::createVariable(pressure));
        EXPECT_THAT(1, gridData.getData(surface).size());
        for (int varIndex = 0; auto surfaceSubVars : surfaceSubData)
          EXPECT_THAT(*surfaceSubVars, VariableFactory::createVariable(surfaceExpectedSubVars[eleIndex][varIndex++]));
      }
      ++eleIndex;
    }
    std::vector<Ikarus::FiniteElements::IFiniteElement> fes;

    for (auto&& ge : surfaces(gridView))
      fes.emplace_back(Ikarus::FiniteElements::ElasticityFE(ge, gridView.indexSet()));

    auto dh = Ikarus::FEManager::DefaultFEManager(fes, gridView, gridData);

    for (int eleIndex = 0, surfIndex = 0; auto&& [fe, eledof, eleVars, eleData] : dh.elementDofsVariableDataTuple()) {
      for (int varIndex = 0; auto&& eleSubData : eleData)
        EXPECT_THAT(*eleSubData, VariableFactory::createVariable(surfaceExpectedSubVars[eleIndex][varIndex++]));
      if (surfIndex == 5) ++eleIndex;
      ++surfIndex;
    }
  }
}