//
// Created by Alex on 21.04.2021.
//
#define EIGEN_MATRIXBASE_PLUGIN "IBB_Eigen_MatrixBaseAddon.h"

#include <gmock/gmock.h>
#include <gtest/gtest.h>

#include <fstream>
#include <vector>

#include <dune/geometry/type.hh>

#include <Eigen/Core>

#include <ikarus/FiniteElements/ElasticityFE.h>
#include <ikarus/FiniteElements/InterfaceFiniteElement.h>
#include <ikarus/Grids/SimpleGrid/SimpleGrid.h>

class TestFE {
public:
  static void initialize() { std::cout << "initcalled" << std::endl; }
  [[nodiscard]] static Ikarus::FiniteElements::IFiniteElement::DofVectorType getEntityVariablePairs() {
    return Ikarus::FiniteElements::IFiniteElement::DofVectorType{};
  }

  static double calculateScalar(const Ikarus::FiniteElements::ElementScalarAffordances&) { return 5; }
};

TEST(FiniteElementInterfaceTest, createGenericFEList) {
  using namespace Ikarus::Grid;
  using namespace Ikarus::FiniteElements;

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

  std::vector<size_t> elementIndices;
  elementIndices.resize(4);
  elementIndices = {0, 1, 2, 3};
  gridFactory.insertElement(Dune::GeometryTypes::quadrilateral, elementIndices);
  elementIndices = {1, 4, 3, 5};
  gridFactory.insertElement(Dune::GeometryTypes::quadrilateral, elementIndices);

  Grid grid = gridFactory.createGrid();

  auto gridView = grid.leafGridView();
  std::vector<Ikarus::FiniteElements::IFiniteElement> fes;

  for (auto&& element : volumes(gridView))
    fes.emplace_back(Ikarus::FiniteElements::ElasticityFE(element));

  Ikarus::DynVectord fint{};
  Ikarus::DynMatrixd K{};
  fint.setZero(8);
  K.setZero(8, 8);
  for (auto&& fe : fes) {
    initialize(fe);

    const auto [KEle, fintEle] = calculateLocalSystem(fe, stiffness, forces);
    EXPECT_EQ(dofSize(fe), 8);
    EXPECT_EQ(calculateVector(fe, forces).size(), 8);
    EXPECT_DOUBLE_EQ(calculateScalar(fe, potentialEnergy), 13.0);
    EXPECT_EQ(calculateMatrix(fe, stiffness).cols(), 8);
    EXPECT_EQ(calculateMatrix(fe, stiffness).rows(), 8);
    EXPECT_THROW(calculateMatrix(fe, mass),std::logic_error);
    EXPECT_EQ(KEle.rows(), 8);
    EXPECT_EQ(KEle.cols(), 8);
    EXPECT_EQ(fintEle.size(), 8);
  }

  Ikarus::FiniteElements::IFiniteElement fe((TestFE()));

  initialize(fe);
  auto entityIDDofPair = getEntityVariablePairs(fe);
  for (auto&& [entityID, var] : entityIDDofPair) {
    std::cout << entityID << " " << var[0] << std::endl;
  }
}