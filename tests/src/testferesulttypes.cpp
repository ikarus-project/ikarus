// SPDX-FileCopyrightText: 2021-2024 The Ikarus Developers mueller@ibb.uni-stuttgart.de
// SPDX-License-Identifier: LGPL-3.0-or-later

#include <config.h>

#include "testcommon.hh"

#include <dune/common/test/testsuite.hh>

#include <ikarus/finiteelements/feresulttypes.hh>
#include <ikarus/utils/init.hh>

using Dune::TestSuite;
using namespace Ikarus;

MAKE_ENUM(ElementHasResultType, Symmetric, Unsymmetric, Vector, Scalar, Custom,Crazy,Non_Square)

// Matrices
REGISTER_SIMPLE_SYMMETRIC_RESULTTYPE(Res1, 3, 3,false);
REGISTER_SIMPLE_RESULTTYPE(Res2, 3, 3);

// Vectors
REGISTER_SIMPLE_RESULTTYPE(Res3, 3, 1);

// Scalars
REGISTER_SIMPLE_RESULTTYPE(Res4, 1, 1);

// Dynamic
REGISTER_SIMPLE_RESULTTYPE(Res5, Eigen::Dynamic, Eigen::Dynamic);

// Fixed Size Max
REGISTER_RESERVED_RESULTTYPE(Res6, Eigen::Dynamic, Eigen::Dynamic, 7, 13);

// Non-Square matrix
REGISTER_SIMPLE_RESULTTYPE(Res7, 3, 2);

template <ElementHasResultType rt, bool asVec>
struct DummyElement
{
  static auto calculateAt()
  requires(rt == ElementHasResultType::Symmetric)
  {
    ResultTypeContainer<Res1<double, 1, 1>, asVec> result;
    if constexpr (asVec)
      result=Eigen::Matrix<double, 6, 1>::Zero();
    else
      result=Eigen::Matrix<double, 3, 3>::Zero();
    return result;
  }
  static auto calculateAt()
  requires(rt == ElementHasResultType::Unsymmetric)
  {
    ResultTypeContainer<Res2<double, 1, 1>, asVec> result;
    if constexpr (asVec)
      result=Eigen::Matrix<double, 9, 1>::Zero();
    else
      result=Eigen::Matrix<double, 3, 3>::Zero();
    return result;
  }
  static auto calculateAt()
  requires(rt == ElementHasResultType::Vector)
  {
    ResultTypeContainer<Res3<double, 1, 1>, true> result;
    result=Eigen::Matrix<double, 3, 1>::Zero();
    return result;
  }
  static auto calculateAt()
  requires(rt == ElementHasResultType::Scalar)
  {
    ResultTypeContainer<Res4<double, 1, 1>, true> result;
    result=Eigen::Matrix<double, 1, 1>::Zero();
    return result;
  }
  static auto calculateAt()
  requires(rt == ElementHasResultType::Custom)
  {
    ResultTypeContainer<Res5<double, 1, 1>, asVec> result;
    if constexpr (asVec)
      result=Eigen::MatrixXd::Zero(9, 1);
    else
      result=Eigen::MatrixXd::Zero(3, 3);
    return result;
  }

  static auto calculateAt()
requires(rt == ElementHasResultType::Crazy)
  {
    ResultTypeContainer<Res6<double, 1, 1>, asVec> result;
    if constexpr (asVec)
      result=Eigen::MatrixXd::Zero(66, 1);
    else
      result=Eigen::MatrixXd::Zero(6, 11);
    return result;
  }

  static auto calculateAt()
requires(rt == ElementHasResultType::Non_Square)
  {
    ResultTypeContainer<Res7<double, 1, 1>, asVec> result;
    if constexpr (asVec)
      result=Eigen::Matrix<double,6,1>::Zero();
    else
      result=Eigen::Matrix<double,3,2>::Zero();
    return result;
  }
};

struct Shape
{
  int rows;
  int cols;
};

auto testRTs() {
  TestSuite t("Test FE ResultTypes");

  auto testStuff = [&]<ElementHasResultType rt, bool asVec>(Shape expectedShapeVec, Shape expectedShapeMat) {
    auto resultContainer = DummyElement<rt, asVec>::calculateAt();
    auto resultAsVec     = resultContainer.asVec();
    auto resultAsMat     = resultContainer.asMat(expectedShapeMat.rows,expectedShapeMat.cols);

    TestSuite tSub("Test Case: " + toString(rt));
    tSub.check(expectedShapeVec.rows == resultAsVec.rows())
        << "Expected rows for Vec is " << expectedShapeVec.rows << " but is " << resultAsVec.rows();
    tSub.check(expectedShapeVec.cols == resultAsVec.cols())
        << "Expected cols for Vec is " << expectedShapeVec.cols << " but is " << resultAsVec.cols();

    tSub.check(expectedShapeMat.rows == resultAsMat.rows())
        << "Expected rows for Mat is " << expectedShapeMat.rows << " but is " << resultAsMat.rows();
    tSub.check(expectedShapeMat.cols == resultAsMat.cols())
        << "Expected cols for Mat is " << expectedShapeMat.cols << " but is " << resultAsMat.cols();
    t.subTest(tSub);
  };

  // Symmetric case
  testStuff.operator()<ElementHasResultType::Symmetric, true>({6, 1}, {3, 3});
  testStuff.operator()<ElementHasResultType::Symmetric, false>({6, 1}, {3, 3});

  // Unsymmetric case
  testStuff.operator()<ElementHasResultType::Unsymmetric, true>({9, 1}, {3, 3});
  testStuff.operator()<ElementHasResultType::Unsymmetric, false>({9, 1}, {3, 3});

  // Vector case
  testStuff.operator()<ElementHasResultType::Vector, true>({3, 1}, {3, 1});

  // Scalar case
  testStuff.operator()<ElementHasResultType::Scalar, true>({1, 1}, {1, 1});

  // Dynamic case
  testStuff.operator()<ElementHasResultType::Custom, false>({9, 1}, {3, 3});
  testStuff.operator()<ElementHasResultType::Custom, true>({9, 1}, {3, 3});

  // Crazy case
  testStuff.operator()<ElementHasResultType::Crazy, false>({66, 1}, {6, 11});
  testStuff.operator()<ElementHasResultType::Crazy, true>({66, 1}, {6, 11});

  // Crazy case
  testStuff.operator()<ElementHasResultType::Non_Square, false>({6, 1}, {3, 2});
  testStuff.operator()<ElementHasResultType::Non_Square, true>({6, 1}, {3, 2});

  return t;
}

int main(const int argc, char** argv) {
  Ikarus::init(argc, argv);
  TestSuite t;

  t.subTest(testRTs());

  return t.exit();
}
