//
// Created by Alex on 21.04.2021.
//

#include <gmock/gmock.h>

#include "testHelpers.h"

#include <fstream>
#include <vector>

#include <Eigen/Core>

#include <ikarus/Assembler/SimpleAssemblers.h>
#include <ikarus/DofManager/DefaultDofManager.h>
#include <ikarus/FiniteElements/ElasticityFE.h>
#include <ikarus/FiniteElements/FiniteElementPolicies.h>
#include <ikarus/FiniteElements/InterfaceFiniteElement.h>
#include <ikarus/Grids/SimpleGrid/SimpleGrid.h>

TEST(Assembler, SimpleAssemblersTest) {
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

  for (auto&& ge : surfaces(gridView))
    fes.emplace_back(Ikarus::FiniteElements::ElasticityFE(ge));

  auto dh = Ikarus::DofManager::DefaultDofManager(fes, gridView);

  auto vectorAssembler = Ikarus::Assembler::VectorAssembler(dh);
  auto fint            = vectorAssembler.getVector(Ikarus::FiniteElements::forces);
  EXPECT_EQ(fint.size(), 12);
  Eigen::VectorXd fintExpected = (Eigen::VectorXd(12) << 1, 1, 2, 2, 1, 1, 2, 2, 1, 1, 1, 1).finished();
  EXPECT_THAT(fint, EigenApproxEqual(fintExpected, 1e-15));

  auto denseMatrixAssembler = Ikarus::Assembler::DenseMatrixAssembler(dh);
  auto K                    = denseMatrixAssembler.getMatrix(Ikarus::FiniteElements::stiffness);
  EXPECT_EQ(K.rows(), 12);
  EXPECT_EQ(K.cols(), 12);
  Eigen::MatrixXd KExpected
      = (Eigen::MatrixXd(12, 12) <<  // clang-format off
  1, 1, 1, 1, 1, 1, 1, 1, 0, 0, 0, 0,
  1, 1, 1, 1, 1, 1, 1, 1, 0, 0, 0, 0,
  1, 1, 2, 2, 1, 1, 2, 2, 1, 1, 1, 1,
  1, 1, 2, 2, 1, 1, 2, 2, 1, 1, 1, 1,
  1, 1, 1, 1, 1, 1, 1, 1, 0, 0, 0, 0,
  1, 1, 1, 1, 1, 1, 1, 1, 0, 0, 0, 0,
  1, 1, 2, 2, 1, 1, 2, 2, 1, 1, 1, 1,
  1, 1, 2, 2, 1, 1, 2, 2, 1, 1, 1, 1,
  0, 0, 1, 1, 0, 0, 1, 1, 1, 1, 1, 1,
  0, 0, 1, 1, 0, 0, 1, 1, 1, 1, 1, 1,
  0, 0, 1, 1, 0, 0, 1, 1, 1, 1, 1, 1,
  0, 0, 1, 1, 0, 0, 1, 1, 1, 1, 1, 1).finished();  // clang-format on
  EXPECT_THAT(K, EigenApproxEqual(KExpected, 1e-15));

  auto sparseMatrixAssembler = Ikarus::Assembler::SparseMatrixAssembler(dh);
  auto KSparse               = sparseMatrixAssembler.getMatrix(Ikarus::FiniteElements::stiffness);

  EXPECT_THAT(KSparse, EigenApproxEqual(KExpected, 1e-15));
  EXPECT_THAT(KSparse, EigenApproxEqual(K, 1e-15));

  auto scalarAssembler = Ikarus::Assembler::ScalarAssembler(dh);
  auto w               = scalarAssembler.getScalar(Ikarus::FiniteElements::potentialEnergy);
  EXPECT_DOUBLE_EQ(w, 26.0);
}
