//
// Created by Alex on 21.04.2021.
//
#define EIGEN_MATRIXBASE_PLUGIN "IBB_Eigen_MatrixBaseAddon.h"

#include <gmock/gmock.h>
#include <gtest/gtest.h>

#include "testHelpers.h"

#include <fstream>
#include <vector>

#include <Eigen/Core>

#include <ikarus/DofHandler/DefaultDofHandler.h>
#include <ikarus/FiniteElements/ElasticityFE.h>
#include <ikarus/FiniteElements/InterfaceFiniteElement.h>
#include <ikarus/Grids/SimpleGrid/SimpleGrid.h>

TEST(DofHandler, DofOwner) {
  using namespace Ikarus::Grid;
  using Grid = SimpleGrid<2, 2>;
  SimpleGridFactory<Grid> gridFactory;
  using vertexType = Ikarus::FixedVector2d;
  std::vector<vertexType> verticesVec;
  verticesVec.emplace_back(vertexType{0.0, 0.0});  // 0
  verticesVec.emplace_back(vertexType{2.0, 0.0});  // 1
  verticesVec.emplace_back(vertexType{0.0, 2.0});  // 2
  verticesVec.emplace_back(vertexType{2.0, 2.0});  // 3
  verticesVec.emplace_back(vertexType{4.0, 0.0});  // 4
  verticesVec.emplace_back(vertexType{4.0, 2.0});  // 5

  for (auto&& vert : verticesVec)
    gridFactory.insertVertex(vert);

  Ikarus::DynArrayXi elementIndices;
  elementIndices.resize(4);
  elementIndices << 0, 1, 2, 3;
  gridFactory.insertElement(Dune::GeometryTypes::quadrilateral, elementIndices);
  elementIndices << 1, 4, 3, 5;
  gridFactory.insertElement(Dune::GeometryTypes::quadrilateral, elementIndices);

  Grid grid = gridFactory.createGrid();

  auto gridView = grid.leafGridView();

  for ([[maybe_unused]] auto&& ge : vertices(gridView))
    std::cout << std::endl;

  Ikarus::Variable::DofOwnerDecorator dof_owner;

  dof_owner.addDof(Ikarus::Variable::DISPLACEMENT3D());
  dof_owner.addDof(Ikarus::Variable::DISPLACEMENT3D());

  update(dof_owner.getDof(Ikarus::Variable::DISPLACEMENT3D()), Eigen::Vector<double, 3>::UnitX());

  EXPECT_THAT(getValue(dof_owner.getDof(Ikarus::Variable::DISPLACEMENT3D())),
              EigenApproxEqual(Eigen::Vector<double, 3>::UnitX(), 1e-15));
  EXPECT_THROW(auto const a = dof_owner.getDof(Ikarus::Variable::DISPLACEMENT1D()), std::logic_error);
}

TEST(DofHandler, DofHandlertest) {
  using namespace Ikarus::Grid;
  using Grid = SimpleGrid<2, 2>;
  SimpleGridFactory<Grid> gridFactory;
  using vertexType = Ikarus::FixedVector2d;
  std::vector<vertexType> verticesVec;
  verticesVec.emplace_back(vertexType{0.0, 0.0});  // 0
  verticesVec.emplace_back(vertexType{2.0, 0.0});  // 1
  verticesVec.emplace_back(vertexType{0.0, 2.0});  // 2
  verticesVec.emplace_back(vertexType{2.0, 2.0});  // 3
  verticesVec.emplace_back(vertexType{4.0, 0.0});  // 4
  verticesVec.emplace_back(vertexType{4.0, 2.0});  // 5

  for (auto&& vert : verticesVec)
    gridFactory.insertVertex(vert);

  Ikarus::DynArrayXi elementIndices;
  elementIndices.resize(4);
  elementIndices << 0, 1, 2, 3;
  gridFactory.insertElement(Dune::GeometryTypes::quadrilateral, elementIndices);
  elementIndices << 1, 4, 3, 5;
  gridFactory.insertElement(Dune::GeometryTypes::quadrilateral, elementIndices);

  Grid grid = gridFactory.createGrid();

  auto gridView = grid.leafGridView();
  std::vector<Ikarus::FiniteElements::IFiniteElement> fes;

  for (auto&& ge : elements(gridView))
    fes.emplace_back(Ikarus::FiniteElements::ElasticityFE(ge));

  auto dh = Ikarus::DofHandler::DefaultDofHandler(fes);
  dh.createDofList();
  //
  //
  ////  u.insert({"red", "#FF0000"});
  ////  u.insert({"green", "#00FF00"});
  ////  u.insert({"blue", "#0000FF"});
  //
  ////  for (auto&& u_ : u) {
  ////    std::cout<<u_.second<<" "<<u_.first<<std::endl;
}
