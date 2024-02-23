// SPDX-FileCopyrightText: 2021-2024 The Ikarus Developers mueller@ibb.uni-stuttgart.de
// SPDX-License-Identifier: LGPL-3.0-or-later

#include <config.h>

#include "testcommon.hh"

#include <dune/common/test/testsuite.hh>

#include <ikarus/finiteelements/feresulttypes.hh>
#include <ikarus/utils/init.hh>

using Dune::TestSuite;
using namespace Ikarus;

MAKE_ENUM(ElementHasResultType, Symmetric, Unsymmetric, Vector, Scalar, Custom)

// Matrices
REGISTER_SYMMETRIC_RESULTTYPE(Res1, 3, 3);
REGISTER_RESULTTYPE(Res2, 3, 3);

// Vectors
REGISTER_RESULTTYPE(Res3, 3, 1);

// Scalars
REGISTER_RESULTTYPE(Res4, 1, 1);

// Dynamic
REGISTER_RESULTTYPE(Res5, Eigen::Dynamic, Eigen::Dynamic);

template <ElementHasResultType rt, bool asVec>
struct DummyElement
{
  static auto calculateAt()
  requires(rt == ElementHasResultType::Symmetric)
  {
    ResultTypeContainer<Res1<double, 1, 1>, asVec> result;
    if constexpr (asVec)
      result.emplace(Eigen::Matrix<double, 6, 1>::Zero());
    else
      result.emplace(Eigen::Matrix<double, 3, 3>::Zero());
    return result;
  }
  static auto calculateAt()
  requires(rt == ElementHasResultType::Unsymmetric)
  {
    ResultTypeContainer<Res2<double, 1, 1>, asVec> result;
    if constexpr (asVec)
      result.emplace(Eigen::Matrix<double, 9, 1>::Zero());
    else
      result.emplace(Eigen::Matrix<double, 3, 3>::Zero());
    return result;
  }
  static auto calculateAt()
  requires(rt == ElementHasResultType::Vector)
  {
    ResultTypeContainer<Res3<double, 1, 1>, true> result;
    result.emplace(Eigen::Matrix<double, 3, 1>::Zero());
    return result;
  }
  static auto calculateAt()
  requires(rt == ElementHasResultType::Scalar)
  {
    ResultTypeContainer<Res4<double, 1, 1>, true> result;
    result.emplace(Eigen::Matrix<double, 1, 1>::Zero());
    return result;
  }
  static auto calculateAt()
  requires(rt == ElementHasResultType::Custom)
  {
    ResultTypeContainer<Res5<double, 1, 1>, asVec> result;
    if constexpr (asVec)
      result.emplace(Eigen::MatrixXd::Zero(9, 1));
    else
      result.emplace(Eigen::MatrixXd::Zero(3, 3));
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
    auto resultAsMat     = resultContainer.asMat();

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

  return t;
}

int main(const int argc, char** argv) {
  Ikarus::init(argc, argv);
  TestSuite t;

  t.subTest(testRTs());

  return t.exit();
}