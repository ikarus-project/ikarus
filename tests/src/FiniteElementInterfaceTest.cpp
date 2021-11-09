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
#include <ikarus/Geometries/GeometryType.h>
#include <ikarus/Grids/SimpleGrid/SimpleGrid.h>

class TestFE {
public:
  static void initialize() {}
  [[nodiscard]] static Ikarus::FiniteElements::IFiniteElement::DofPairVectorType getEntityVariablePairs() {
    return Ikarus::FiniteElements::IFiniteElement::DofPairVectorType{};
  }

  static double calculateScalar(const Ikarus::FiniteElements::ScalarAffordances&, Ikarus::FEValues&,
                                std::optional<std::reference_wrapper<Ikarus::FEValues>>& data) {
    if (data)
      if (isType(data->get().get(Ikarus::EntityType::vertex)[0], Ikarus::Variable::VariablesTags::displacement2d))
        return getValue(data->get().get(Ikarus::EntityType::vertex)[0])[0]
               * getValue(data->get().get(Ikarus::EntityType::vertex)[0])[1];
    return 5;
  }
};

class TestFE2 {
public:
  static void initialize() {}
  [[nodiscard]] static Ikarus::FiniteElements::IFiniteElement::DofPairVectorType getEntityVariablePairs() {
    return Ikarus::FiniteElements::IFiniteElement::DofPairVectorType{};
  }

  static double calculateScalar(const Ikarus::FiniteElements::ScalarAffordances&, Ikarus::FEValues&) { return 5; }
};

TEST(FiniteElementInterfaceTest, createGenericFEList) {
  using namespace Ikarus::Grid;
  using namespace Ikarus::FiniteElements;

  using Grid = SimpleGrid<2, 2>;
  SimpleGridFactory<2, 2> gridFactory;
  //  using vertexType = Eigen::Vector2d;
  std::vector<Eigen::Vector2d> verticesVec;
  verticesVec.emplace_back(0.0, 0.0);  // 0
  verticesVec.emplace_back(2.0, 0.0);  // 1
  verticesVec.emplace_back(0.0, 2.0);  // 2
  verticesVec.emplace_back(2.0, 2.0);  // 3
  verticesVec.emplace_back(4.0, 0.0);  // 4
  verticesVec.emplace_back(4.0, 2.0);  // 5

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

  Eigen::VectorXd fint{};
  Eigen::MatrixXd K{};
  fint.setZero(8);
  K.setZero(8, 8);
  for (auto&& fe : fes) {
    initialize(fe);
    using namespace Ikarus::Variable;
    std::vector<IVariable> vars;
    vars.emplace_back(VariableFactory::createVariable(displacement2d));
    vars.emplace_back(VariableFactory::createVariable(displacement2d));
    vars.emplace_back(VariableFactory::createVariable(displacement2d));
    vars.emplace_back(VariableFactory::createVariable(displacement2d));

    std::vector<IVariable*> varsP;
    varsP.resize(4);
    for (int i = 0; auto& varP : varsP)
      varP = &vars[i++];
    Ikarus::FEValues feValues;
    feValues.set(Ikarus::EntityType::vertex, vars);

    const auto [KEle, fintEle] = calculateLocalSystem(fe, stiffness, forces, feValues);
    EXPECT_EQ(dofSize(fe), 8);
    EXPECT_EQ(calculateVector(fe, forces, feValues).size(), 8);
    EXPECT_DOUBLE_EQ(calculateScalar(fe, potentialEnergy, feValues), 13.0);
    EXPECT_EQ(calculateMatrix(fe, stiffness, feValues).cols(), 8);
    EXPECT_EQ(calculateMatrix(fe, stiffness, feValues).rows(), 8);
    EXPECT_THROW(calculateMatrix(fe, mass, feValues), std::logic_error);
    EXPECT_THROW(calculateLocalSystem(fe, mass, forces, feValues), std::logic_error);
    EXPECT_EQ(KEle.rows(), 8);
    EXPECT_EQ(KEle.cols(), 8);
    EXPECT_EQ(fintEle.size(), 8);
  }

  Ikarus::FiniteElements::IFiniteElement fe((TestFE()));

  initialize(fe);
  const auto entityIDDofPair = getEntityVariableTuple(fes[0]);
  std::vector<std::pair<size_t, Ikarus::Variable::VariablesTags>> idtagExpected;
  idtagExpected.emplace_back(0, Ikarus::Variable::displacement2d);
  idtagExpected.emplace_back(1, Ikarus::Variable::displacement2d);
  idtagExpected.emplace_back(2, Ikarus::Variable::displacement2d);
  idtagExpected.emplace_back(3, Ikarus::Variable::displacement2d);
  for (int i = 0; auto&& [entityID, entityType, varVec] : entityIDDofPair) {
    EXPECT_EQ(entityID, idtagExpected[i].first);
    EXPECT_EQ(varVec.size(), 1);
    EXPECT_EQ(entityType, Ikarus::EntityType::vertex);
    EXPECT_EQ(varVec[0], idtagExpected[i].second);
    ++i;
  }

  auto feT{fes[0]};  // test copy assignment

  using namespace Ikarus::Variable;
  std::vector<IVariable> vars;
  vars.emplace_back(VariableFactory::createVariable(displacement2d));
  vars.emplace_back(VariableFactory::createVariable(displacement2d));
  vars.emplace_back(VariableFactory::createVariable(displacement2d));
  vars.emplace_back(VariableFactory::createVariable(displacement2d));

  std::vector<IVariable*> varsP;
  varsP.resize(4);
  for (int i = 0; auto& varP : varsP)
    varP = &vars[i++];

  Ikarus::FEValues feValues;
  feValues.set(Ikarus::EntityType::vertex, vars);

  std::vector<IVariable> datas;
  datas.emplace_back(VariableFactory::createVariable(displacement2d));
  datas[0] += Eigen::Vector2d(15, 2);
  std::vector<IVariable*> datasP;
  datasP.resize(1);
  for (int i = 0; auto& dataP : datasP)
    dataP = &datas[i++];

  Ikarus::FEValues feDataValues;
  feDataValues.set(Ikarus::EntityType::vertex, datas);

  EXPECT_THAT(feDataValues.get(Ikarus::EntityType::vertex).size(), 1);
  EXPECT_THAT(feDataValues.get(Ikarus::EntityType::edge).size(), 0);
  EXPECT_THAT(feDataValues.get(Ikarus::EntityType::surface).size(), 0);
  EXPECT_THAT(feDataValues.get(Ikarus::EntityType::volume).size(), 0);
  EXPECT_THAT(feDataValues.get(Ikarus::EntityType::generic).size(), 0);

  EXPECT_DOUBLE_EQ(calculateScalar(fe, potentialEnergy, feValues, feDataValues), 30.0);
  datas[0] = VariableFactory::createVariable(displacement1d);
  EXPECT_DOUBLE_EQ(calculateScalar(fe, potentialEnergy, feValues, feDataValues), 5.0);

  Ikarus::FiniteElements::IFiniteElement fe2((TestFE2()));  // check if element without optional data is accepted

  initialize(fe2);
  EXPECT_DOUBLE_EQ(calculateScalar(fe2, potentialEnergy, feValues), 5.0);
}