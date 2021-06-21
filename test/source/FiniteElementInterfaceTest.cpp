//
// Created by Alex on 21.04.2021.
//
#define EIGEN_MATRIXBASE_PLUGIN "IBB_Eigen_MatrixBaseAddon.h"
#include <array>
#include <fstream>
#include <gtest/gtest.h>
#include <vector>

#include <Eigen/Core>

//
//
#include <dune/geometry/type.hh>

#include <ikarus/FiniteElements/ElasticityFE.h>
#include <ikarus/FiniteElements/GenericFiniteElement.h>
#include <ikarus/Geometries/GeometryWithExternalInput.h>
#include <ikarus/Grids/GridEntities/DefaultGridEntities.h>
#include <ikarus/Grids/SimpleGrid/SimpleGrid.h>

class TestFE {
public:
  void initialize() { std::cout << "initcalled" << std::endl; }
};

TEST(FiniteElementInterfaceTest, createGenericFEList) {
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

  for (auto&& vert : verticesVec) gridFactory.insertVertex(vert);

  Ikarus::DynArrayXi elementIndices;
  elementIndices.resize(4);
  elementIndices << 0, 1, 2, 3;
  gridFactory.insertElement(Dune::GeometryTypes::quadrilateral, elementIndices);
  elementIndices << 1, 4, 3, 5;
  gridFactory.insertElement(Dune::GeometryTypes::quadrilateral, elementIndices);

  Grid grid = gridFactory.createGrid();

  auto gridView = grid.leafGridView();
  std::vector<Ikarus::PhysicalElements::GenericFE> fes;

  for (auto&& element : elements(gridView)) fes.emplace_back(Ikarus::PhysicalElements::ElasticityFE(element));

  Ikarus::DynVectord fint{};
  Ikarus::DynMatrixd K{};
  fint.setZero(5);
  K.setZero(5, 5);
  for (auto&& fe : fes) {
    initialize(fe);
    //      std::cout<<"DofSize: "<< dofSize(fe)<<std::endl;
    const auto [fintEle, KEle] = calculateLocalSystem(fe);
    fint += calculateRHS(fe);
    K += calculateLHS(fe);
    std::cout << dofSize(fe) << std::endl;
  }

  Ikarus::PhysicalElements::GenericFE fe((TestFE()));

  initialize(fe);
  //  getDofVector(fe);
  EXPECT_THROW(getDofVector(fe), Dune::InvalidStateException);
}