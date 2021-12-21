//
// Created by Alex on 21.04.2021.
//

#include <gmock/gmock.h>

#include "testHelpers.h"

#include <fstream>
#include <vector>

#include <Eigen/Core>
#include <ikarus/Assembler/SimpleAssemblers.h>

#include "ikarus/LinearAlgebra/DirichletConditionManager.h"
#include <ikarus/FEManager/DefaultFEManager.h>
#include <ikarus/FiniteElements/ElasticityFE.h>
#include <ikarus/FiniteElements/FiniteElementFunctionConcepts.h>
#include <ikarus/FiniteElements/InterfaceFiniteElement.h>
#include <ikarus/Geometries/GeometryType.h>
#include <ikarus/Grids/SimpleGrid/SimpleGrid.h>

TEST(Assembler, SimpleAssemblersTest) {
  using namespace Ikarus::Grid;
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

  const auto indexSet = gridView.indexSet();

  std::vector<Ikarus::FiniteElements::IFiniteElement> fes;

  for (auto&& ge : surfaces(gridView))
    fes.emplace_back(Ikarus::FiniteElements::ElasticityFE(ge, indexSet, 1000, 0.3));

  auto feManager = Ikarus::FEManager::DefaultFEManager(fes, gridView);

  Ikarus::DirichletConditionManager dirichletConditionManager(feManager);

  dirichletConditionManager.addConstraint(vertices(gridView).front(), 0);
  dirichletConditionManager.addConstraint(vertices(gridView).back(), 1);
  dirichletConditionManager.addConstraint(vertices(gridView).at(3), 1);

  auto vectorAssembler = Ikarus::Assembler::VectorAssembler(feManager, dirichletConditionManager);
  auto fint            = vectorAssembler.getVector(Ikarus::FiniteElements::forces);
  EXPECT_EQ(fint.size(), 12);
  const Eigen::VectorXd fintExpected = (Eigen::VectorXd(12) << 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0).finished();
  EXPECT_THAT(fint, EigenApproxEqual(fintExpected, 1e-15));

  auto denseMatrixAssembler = Ikarus::Assembler::DenseMatrixAssembler(feManager, dirichletConditionManager);
  auto K                    = denseMatrixAssembler.getMatrix(Ikarus::FiniteElements::stiffness);
  EXPECT_EQ(K.rows(), 12);
  EXPECT_EQ(K.cols(), 12);
  const Eigen::MatrixXd KExpected
      = (Eigen::MatrixXd(12, 12) <<  // clang-format off
  494.5054945054946,  178.5714285714286, -302.1978021978022, -13.73626373626373,  54.94505494505495,  13.73626373626373, -247.2527472527473, -178.5714285714286,                  0,                  0,                  0,                  0,
  178.5714285714286,  494.5054945054946,  13.73626373626373,  54.94505494505495, -13.73626373626373, -302.1978021978022, -178.5714285714286, -247.2527472527473,                  0,                  0,                  0,                  0,
 -302.1978021978023,  13.73626373626372,  989.0109890109891,                  0, -247.2527472527473,  178.5714285714286,  109.8901098901099, -7.10542735760e-15, -302.1978021978022, -13.73626373626373, -247.2527472527473, -178.5714285714286,
 -13.73626373626373,  54.94505494505495,                  0,  989.0109890109891,  178.5714285714286, -247.2527472527473, -3.55271367880e-15, -604.3956043956044,  13.73626373626373,  54.94505494505495, -178.5714285714286, -247.2527472527473,
  54.94505494505495, -13.73626373626373, -247.2527472527473,  178.5714285714286,  494.5054945054945, -178.5714285714286, -302.1978021978022,  13.73626373626374,                  0,                  0,                  0,                  0,
  13.73626373626373, -302.1978021978023,  178.5714285714286, -247.2527472527473, -178.5714285714286,  494.5054945054945, -13.73626373626373,  54.94505494505495,                  0,                  0,                  0,                  0,
 -247.2527472527473, -178.5714285714286,  109.8901098901099, -2.48689957516e-14, -302.1978021978022, -13.73626373626373,  989.0109890109891,                  0, -247.2527472527473,  178.5714285714286, -302.1978021978022,  13.73626373626374,
 -178.5714285714286, -247.2527472527473, -1.06581410364e-14, -604.3956043956046,  13.73626373626373,  54.94505494505495,                  0,  989.0109890109891,  178.5714285714286, -247.2527472527473, -13.73626373626372,  54.94505494505496,
                  0,                  0, -302.1978021978023,  13.73626373626373,                  0,                  0, -247.2527472527473,  178.5714285714286,  494.5054945054945, -178.5714285714286,  54.94505494505495, -13.73626373626373,
                  0,                  0, -13.73626373626373,  54.94505494505495,                  0,                  0,  178.5714285714286, -247.2527472527473, -178.5714285714286,  494.5054945054945,  13.73626373626374, -302.1978021978022,
                  0,                  0, -247.2527472527473, -178.5714285714286,                  0,                  0, -302.1978021978022, -13.73626373626372,  54.94505494505495,  13.73626373626373,  494.5054945054945,  178.5714285714286,
                  0,                  0, -178.5714285714286, -247.2527472527473,                  0,                  0,  13.73626373626373,  54.94505494505496, -13.73626373626373, -302.1978021978022,  178.5714285714286,  494.5054945054945).finished();  // clang-format on
  EXPECT_THAT(K, EigenApproxEqual(KExpected, 1e-15));

  auto sparseMatrixAssembler = Ikarus::Assembler::SparseMatrixAssembler(feManager);
  auto KSparse               = sparseMatrixAssembler.getMatrix(Ikarus::FiniteElements::stiffness);

  EXPECT_THAT(KSparse, EigenApproxEqual(KExpected, 1e-15));
  EXPECT_THAT(KSparse, EigenApproxEqual(K, 1e-15));

  auto scalarAssembler = Ikarus::Assembler::ScalarAssembler(feManager);
  auto w               = scalarAssembler.getScalar(Ikarus::FiniteElements::potentialEnergy);
  EXPECT_DOUBLE_EQ(w, 26.0);

  // Reduced tests
  const auto& fintRed = vectorAssembler.getReducedVector(Ikarus::FiniteElements::forces);

  EXPECT_EQ(fintRed.size(), fint.size() - 3);

  Eigen::VectorXd testVector = Eigen::VectorXd::LinSpaced(9, 1, 9);
  const auto fullVector      = vectorAssembler.createFullVector(testVector);

  EXPECT_EQ(fullVector.size(), fint.size());
  const Eigen::VectorXd fullVectorExpected = (Eigen::VectorXd(12) << 0, 1, 2, 3, 4, 5, 6, 0, 7, 8, 9, 0).finished();
  EXPECT_THAT(fullVector, EigenApproxEqual(fullVectorExpected, 1e-15));

  const auto& KRed = denseMatrixAssembler.getReducedMatrix(Ikarus::FiniteElements::stiffness);

  Eigen::MatrixXd KExpectedRed = KExpected;
  std::vector<size_t> keepIndices(dirichletConditionManager.freeIndices().begin(),
                                  dirichletConditionManager.freeIndices().end());

  KExpectedRed = KExpectedRed(keepIndices, keepIndices).eval();
//  KExpectedRed = KExpectedRed(Eigen::all, keepIndices).eval();

  EXPECT_THAT(KRed, EigenApproxEqual(KExpectedRed, 1e-15));
}
