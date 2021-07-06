//
// Created by Alex on 21.04.2021.
//

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
  [[nodiscard]] static Ikarus::FiniteElements::IFiniteElement::DofPairVectorType getEntityVariablePairs() {
    return Ikarus::FiniteElements::IFiniteElement::DofPairVectorType{};
  }

  static double calculateScalar(std::vector<Ikarus::Variable::IVariable*>,
                                const Ikarus::FiniteElements::ScalarAffordances&) {
    return 5;
  }
};

TEST(FiniteElementInterfaceTest, createGenericFEList) {
  using namespace Ikarus::Grid;
  using namespace Ikarus::FiniteElements;

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

  for (auto&& element : volumes(gridView))
    fes.emplace_back(Ikarus::FiniteElements::ElasticityFE(element));

  Eigen::VectorXd fint{};
  Eigen::MatrixXd K{};
  fint.setZero(8);
  K.setZero(8, 8);
  for (auto&& fe : fes) {
    initialize(fe);

    std::vector<Ikarus::Variable::IVariable> vars;
    vars.emplace_back(Ikarus::Variable::VariableFactory::createVariable(Ikarus::Variable::displacement2d));
    vars.emplace_back(Ikarus::Variable::VariableFactory::createVariable(Ikarus::Variable::displacement2d));
    vars.emplace_back(Ikarus::Variable::VariableFactory::createVariable(Ikarus::Variable::displacement2d));
    vars.emplace_back(Ikarus::Variable::VariableFactory::createVariable(Ikarus::Variable::displacement2d));

    std::vector<Ikarus::Variable::IVariable*> varsP;
    varsP.resize(4);
    for (int i = 0; auto& varP : varsP)
      varP = &vars[i++];

    const auto [KEle, fintEle] = calculateLocalSystem(fe, varsP, stiffness, forces);
    EXPECT_EQ(dofSize(fe), 8);
    EXPECT_EQ(calculateVector(fe, varsP, forces).size(), 8);
    EXPECT_DOUBLE_EQ(calculateScalar(fe, varsP, potentialEnergy), 13.0);
    EXPECT_EQ(calculateMatrix(fe, varsP, stiffness).cols(), 8);
    EXPECT_EQ(calculateMatrix(fe, varsP, stiffness).rows(), 8);
    EXPECT_THROW(calculateMatrix(fe, varsP, mass), std::logic_error);
    EXPECT_THROW(calculateLocalSystem(fe, varsP, mass, forces), std::logic_error);
    EXPECT_EQ(KEle.rows(), 8);
    EXPECT_EQ(KEle.cols(), 8);
    EXPECT_EQ(fintEle.size(), 8);
  }

Ikarus::FiniteElements::IFiniteElement fe((TestFE()));

initialize(fe);
const auto entityIDDofPair = getEntityVariablePairs(fes[0]);
std::vector<std::pair<size_t, Ikarus::Variable::VariablesTags>> idtagExpected;
idtagExpected.emplace_back(0, Ikarus::Variable::displacement2d);
idtagExpected.emplace_back(1, Ikarus::Variable::displacement2d);
idtagExpected.emplace_back(2, Ikarus::Variable::displacement2d);
idtagExpected.emplace_back(3, Ikarus::Variable::displacement2d);
for (int i = 0; auto&& [entityID, var] : entityIDDofPair) {
  EXPECT_EQ(entityID, idtagExpected[i].first);
  EXPECT_EQ(var.size(), 1);
  EXPECT_EQ(var[0], idtagExpected[i].second);
  ++i;
}

auto feT{fes[0]};  // test copy assignment
}