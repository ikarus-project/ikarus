//
// Created by Alex on 21.04.2021.
//

#include <gmock/gmock.h>

#include "testHelpers.h"

#include <fstream>
#include <vector>

#include <Eigen/Core>

#include <ikarus/DofManager/DefaultDofManager.h>
#include <ikarus/FiniteElements/ElasticityFE.h>
#include <ikarus/FiniteElements/InterfaceFiniteElement.h>
#include <ikarus/Grids/SimpleGrid/SimpleGrid.h>
#include <ikarus/Assembler/ResidualAssembler.h>


TEST(Assembler, VectorAssemblerTest) {
  using namespace Ikarus::Grid;
  using Grid = SimpleGrid<2, 2>;
  SimpleGridFactory<Grid> gridFactory;
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
  gridFactory.insertElement(Dune::GeometryTypes::quadrilateral, elementIndices);
  elementIndices = {1, 4, 3, 5};
  gridFactory.insertElement(Dune::GeometryTypes::quadrilateral, elementIndices);

  Grid grid = gridFactory.createGrid();

  auto gridView = grid.leafGridView();

  std::vector<Ikarus::FiniteElements::IFiniteElement> fes;

  for (auto&& ge : volumes(gridView))
    fes.emplace_back(Ikarus::FiniteElements::ElasticityFE(ge));

  auto dh = Ikarus::DofManager::DefaultDofManager(fes, gridView);
  dh.createElementDofRelationship();

  //auto vectorAssembler = Ikarus::Assembler:VectorAssembler(dh);

}
