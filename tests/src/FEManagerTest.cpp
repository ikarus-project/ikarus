//
// Created by Alex on 21.04.2021.
//

#include <gmock/gmock.h>

#include "testHelpers.h"

#include <fstream>
#include <vector>

#include <Eigen/Core>

#include <ikarus/FEManager/DefaultFEManager.h>
#include <ikarus/FiniteElements/ElasticityFE.h>
#include <ikarus/FiniteElements/InterfaceFiniteElement.h>
#include <ikarus/Geometries/GeometryType.h>
#include <ikarus/Grids/SimpleGrid/SimpleGrid.h>

TEST(FEManager, FEManagertest) {
  using namespace Ikarus::Grid;
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

  std::vector<Ikarus::FiniteElements::IFiniteElement> feContainer;

  for (auto&& ge : surfaces(gridView))
    feContainer.emplace_back(Ikarus::FiniteElements::ElasticityFE(ge, gridView.indexSet(), 1000, 0.3));

  auto feManager = Ikarus::FEManager::DefaultFEManager(feContainer, gridView);

  auto VariablesOfAllElements            = feManager.elementVariables();
  auto VariablesOfFirstElement           = VariablesOfAllElements[0];
  auto VariablesAtVerticesOfFirstElement = VariablesOfFirstElement.get(Ikarus::EntityType::vertex);

  // throws since there was no data defined should this be the case?
  EXPECT_THROW(feManager.elementIndicesVariableDataTuple(), std::logic_error);
  // TODO richtiges Benutzem im Test auch testen

  GridData gridData(gridView.indexSet());
  gridData.add(data(VariableTags::velocity2d), Ikarus::EntityType::vertex);
  feManager.addData(gridData);

  for (auto&& [fe, dofIndices, vars, data] : feManager.elementIndicesVariableDataTuple())
    EXPECT_TRUE(isType(data.get(Ikarus::EntityType::vertex)[0].get(), VariableTags::velocity2d));

  for (auto&& disp2DAtVertex : VariablesAtVerticesOfFirstElement)
    disp2DAtVertex += Eigen::Vector<double, 2>::UnitX();

  EXPECT_THAT(VariablesAtVerticesOfFirstElement.size(), 4);

  for (auto&& eleDofVecSize : feManager.elementDofVectorSize()) {
    EXPECT_THAT(eleDofVecSize, 8);
  }

  auto VariableList2 = feManager.elementVariables()[1].get(Ikarus::EntityType::vertex);

  EXPECT_THAT(getValue(VariableList2[0]), EigenApproxEqual(Eigen::Vector<double, 2>::UnitX(), 1e-15));
  EXPECT_THAT(getValue(VariableList2[1]), EigenApproxEqual(Eigen::Vector<double, 2>::Zero(), 1e-15));
  EXPECT_THAT(getValue(VariableList2[2]), EigenApproxEqual(Eigen::Vector<double, 2>::UnitX(), 1e-15));
  EXPECT_THAT(getValue(VariableList2[3]), EigenApproxEqual(Eigen::Vector<double, 2>::Zero(), 1e-15));

  for (auto&& disp2DAtVertex : VariablesAtVerticesOfFirstElement)
    disp2DAtVertex -= Eigen::Vector<double, 2>::UnitY();

  EXPECT_THAT(getValue(VariableList2[0]), EigenApproxEqual(Eigen::Vector<double, 2>(1, -1), 1e-15));
  EXPECT_THAT(getValue(VariableList2[1]), EigenApproxEqual(Eigen::Vector<double, 2>::Zero(), 1e-15));
  EXPECT_THAT(getValue(VariableList2[2]), EigenApproxEqual(Eigen::Vector<double, 2>(1, -1), 1e-15));
  EXPECT_THAT(getValue(VariableList2[3]), EigenApproxEqual(Eigen::Vector<double, 2>::Zero(), 1e-15));

  EXPECT_THAT(feManager.numberOfDegreesOfFreedom(), vertices(gridView).size() * 2);
  Eigen::VectorXd D(vertices(gridView).size() * 2);

  auto& x = feManager.getVariables();

  D = Eigen::VectorXd::LinSpaced(D.size(), 0, D.size() - 1);

  std::vector<Eigen::Vector2d> xExpected(6);
  for (int i = 0; i < 6; ++i) {
    if (i < 4)
      xExpected[i] << 1 + D[2 * i], -1 + D[2 * i + 1];
    else
      xExpected[i] << D[2 * i], D[2 * i + 1];
  }

  x += D;

  // check if all variables have the correct values
  for (int i = 0; auto& var : x.getValues())
    EXPECT_THAT(getValue(var), EigenApproxEqual(xExpected[i++], 1e-15));

  auto dofIndicesOfFirstElement  = feManager.elementDofs()[0];
  auto dofIndicesOfSecondElement = feManager.elementDofs()[1];

  std::array<Eigen::ArrayX<size_t>, 2> expectedIndices;
  expectedIndices[0].resize(8);
  expectedIndices[1].resize(8);
  std::iota(expectedIndices[0].begin(), expectedIndices[0].end(), 0);
  expectedIndices[1] << 2, 3, 8, 9, 6, 7, 10, 11;

  EXPECT_THAT(dofIndicesOfFirstElement, EigenExactEqual(expectedIndices[0]));
  EXPECT_THAT(dofIndicesOfSecondElement, EigenExactEqual(expectedIndices[1]));
}
