//
// Created by Alex on 21.04.2021.
//
#include <gmock/gmock.h>
#include <gtest/gtest.h>

#include "testHelpers.h"

#include <array>
#include <fstream>
#include <unordered_set>
#include <vector>

#include <Eigen/Core>

#include "ikarus/Manifolds/RealTuple.h"
#include "ikarus/Variables/InterfaceVariable.h"
#include "ikarus/Variables/VariableDefinitions.h"
#include "ikarus/utils/std/algorithms.h"

TEST(DefaultVariableTest, RealTupleDisplacement) {
  using namespace Ikarus::Variable;
  auto a = VariableFactory::createVariable(Ikarus::Variable::displacement3d);

  a += Eigen::Vector<double, 3>::UnitX();
  EXPECT_EQ(getValue(a), (Eigen::Vector<double, 3>::UnitX()));

  a += Eigen::Vector<double, 3>::UnitY();
  EXPECT_EQ(getValue(a), (Eigen::Vector<double, 3>(1.0, 1.0, 0.0)));

  auto d = a + Eigen::Vector<double, 3>::UnitY();
  EXPECT_EQ(getValue(d), (Eigen::Vector<double, 3>(1.0, 2.0, 0.0)));

  auto b{a};
  EXPECT_EQ(getValue(b), (Eigen::Vector<double, 3>(1.0, 1.0, 0.0)));

  DISPLACEMENT3D c{DISPLACEMENT3D{Eigen::Vector<double, 3>::UnitZ()}};  // move constructor test
  EXPECT_EQ(c.getValue(), (Eigen::Vector<double, 3>(0.0, 0.0, 1.0)));

  c.setValue(Eigen::Vector<double, 3>(13.0, -5.0, 1.0));
  EXPECT_EQ(c.getValue(), (Eigen::Vector<double, 3>(13.0, -5.0, 1.0)));

  b = a;
  EXPECT_EQ(b, a);
  EXPECT_EQ(getValue(b), getValue(a));
  auto testVec = Eigen::Vector<double, 3>(127.0, -5.0, 1.0);
  setValue(b, testVec);

  EXPECT_EQ(getValue(b), testVec);
}

static constexpr double tol = 1e-15;

TEST(DefaultVariableTest, UnitVectorDirector) {
  using namespace Ikarus::Variable;
  DIRECTOR3D a{DIRECTOR3D::CoordinateType::UnitZ()};
  a.update(Eigen::Vector<double, 2>::UnitX());
  const auto aExpected = Eigen::Vector<double, 3>(1.0 / sqrt(2), 0.0, 1.0 / sqrt(2));
  EXPECT_THAT(a.getValue(), EigenApproxEqual(aExpected, tol));

  a = update(a, Eigen::Vector<double, 2>::UnitY());
  EXPECT_THAT(a.getValue(), EigenApproxEqual(Eigen::Vector<double, 3>(0.5, 1.0 / sqrt(2), 0.5), tol));

  const auto d = update(a, Eigen::Vector<double, 2>::UnitY());

  EXPECT_THAT(d.getValue(), EigenApproxEqual(Eigen::Vector<double, 3>(0.18688672392660707344, 0.97140452079103167815,
                                                                      -0.14644660940672621363),
                                             tol));

  DIRECTOR3D b{a};

  EXPECT_THAT(b.getValue(), EigenApproxEqual(Eigen::Vector<double, 3>(0.5, 1.0 / sqrt(2), 0.5), tol));

  DIRECTOR3D c{DIRECTOR3D{Eigen::Vector<double, 3>::UnitZ() * 2.0}};  // move constructor test
  EXPECT_THAT(c.getValue(), EigenApproxEqual(Eigen::Vector<double, 3>(0.0, 0.0, 1.0), tol));

  c.setValue(Eigen::Vector<double, 3>(13.0, -5.0, 1.0));
  EXPECT_THAT(c.getValue(), EigenApproxEqual(Eigen::Vector<double, 3>(13.0, -5.0, 1.0).normalized(), tol));

  b = a;
  EXPECT_EQ(b, a);
  const auto testVec = Eigen::Vector<double, 3>(127.0, -5.0, 1.0);
  b.setValue(testVec);

  EXPECT_THAT(b.getValue(), EigenApproxEqual(testVec.normalized(), tol));
}

TEST(VariableTest, GenericVariableVectorTest) {
  using namespace Ikarus::Variable;
  DISPLACEMENT3D a;
  DISPLACEMENT2D b;
  DISPLACEMENT1D c;
  DISPLACEMENT3D d;
  DISPLACEMENT2D e;
  DIRECTOR3D f{DIRECTOR3D::CoordinateType::UnitZ()};

  std::vector<IVariable> varVec;
  varVec.emplace_back(a);
  varVec.emplace_back(b);
  varVec.emplace_back(c);
  varVec.emplace_back(d);
  varVec.emplace_back(e);

  EXPECT_EQ(valueSize(varVec), 11);

  EXPECT_EQ(correctionSize(varVec), 11);
  varVec.emplace_back(f);
  EXPECT_EQ(valueSize(varVec), 14);
  EXPECT_EQ(correctionSize(varVec), 13);

  Eigen::Vector2d corr;
  corr << 1, 1;
  varVec[5] += corr;

  Eigen::Vector3d varVec4Expected({1, 1, 1});
  varVec4Expected.normalize();
  EXPECT_THAT(getValue(varVec[5]), EigenApproxEqual(varVec4Expected, tol));

  Ikarus::stl::makeUniqueAndSort(varVec);

  EXPECT_EQ(valueSize(varVec), 9);
  EXPECT_EQ(correctionSize(varVec), 8);

  Eigen::VectorXd correction(correctionSize(varVec));

  correction << 1, 2, 3, 4, 5, 6, 7, 8;
  update(varVec, correction);

  EXPECT_THAT(getValue(varVec[0]), EigenApproxEqual(Eigen::Matrix<double, 1, 1>(1), tol));
  EXPECT_THAT(getValue(varVec[1]), EigenApproxEqual(Eigen::Vector2d(2, 3), tol));
  EXPECT_THAT(getValue(varVec[2]), EigenApproxEqual(Eigen::Vector3d(4, 5, 6), tol));
  EXPECT_THAT(
      getValue(varVec[3]),
      EigenApproxEqual(Eigen::Vector3d(0.41279806929140905325, 0.50645665044957854928, -0.75703329861022516933), tol));
}

#include <ikarus/Grids/GridEntities/DefaultGridEntities.h>
TEST(VariableTest, GenericVariableWithFEandGridEntity) {
  using namespace Ikarus::Variable;
  DISPLACEMENT3D a;

  Ikarus::Grid::DefaultGridEntity<3, 0, 3> gridEntity;
}