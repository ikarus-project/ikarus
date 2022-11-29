// SPDX-FileCopyrightText: 2022 The Ikarus Developers mueller@ibb.uni-stuttgart.de
// SPDX-License-Identifier: LGPL-2.1-or-later

//
//
#include <config.h>

#include <dune/common/parallel/mpihelper.hh>
#include <dune/common/test/testsuite.hh>
using Dune::TestSuite;

#include "testHelpers.hh"

#include <ikarus/localFunctions/meta.hh>
#include <ikarus/utils/traits.hh>
auto localFunctionTraitsTest() {
  TestSuite t("localFunctionTraitsTest");
  using namespace Ikarus::DerivativeDirections;
  auto wrt1 = Ikarus::wrt(spatial(0), coeff(7));

  auto counter = countDerivativesType<decltype(wrt1)>();

  t.check(1 == counter.singleCoeffDerivs);
  t.check(0 == counter.twoCoeffDerivs);
  t.check(1 == counter.spatialDerivs);
  t.check(0 == counter.spatialAllCounter);

  auto wrt2 = Ikarus::wrt(spatialAll, coeff(7, 1));

  counter = countDerivativesType<decltype(wrt2)>();

  t.check(0 == counter.singleCoeffDerivs);
  t.check(1 == counter.twoCoeffDerivs);
  t.check(0 == counter.spatialDerivs);
  t.check(1 == counter.spatialAllCounter);
  return t;
}

auto treePathCoeffs() {
  TestSuite t("TreePathCoeffs");
  using namespace Ikarus::DerivativeDirections;
  using namespace Dune::Indices;
  auto coeffTwo = coeff(_0, 7, _1, 9);

  t.check(0 == coeffTwo.index[_0][0]);
  t.check(7 == coeffTwo.index[_0][1]);
  t.check(1 == coeffTwo.index[_1][0]);
  t.check(9 == coeffTwo.index[_1][1]);

  auto coeffSingle = coeff(_4, 1);

  t.check(4 == coeffSingle.index[_0][0]);
  t.check(1 == coeffSingle.index[_0][1]);
  return t;
}

// Move this to compile only Tests
// auto testTupleFilter() {
//   TestSuite t("testTupleFilter");
//   using namespace Ikarus::DerivativeDirections;
//   using namespace Dune::Indices;
//   auto tup                 = std::make_tuple(_0, _1, _3, Ikarus::arithmetic);
//   auto filteredTupExpected = std::make_tuple(_0, _1, _3);
//   constexpr auto predicate = []<typename Type>(Type) { return Type::value != Ikarus::arithmetic; };
//   auto filteredTup         = Ikarus::Std::filter(tup, predicate);
//
//   static_assert(std::is_same_v<decltype(filteredTup), decltype(filteredTupExpected)>);
//
//   auto tup2         = std::make_tuple(Ikarus::arithmetic, _0, Ikarus::arithmetic, _1, _3, Ikarus::arithmetic);
//   auto filteredTup2 = Ikarus::Std::filter(tup2, predicate);
//   static_assert(std::is_same_v<decltype(filteredTup2), decltype(filteredTupExpected)>);
// }

auto makeNestedTupleFlat() {
  TestSuite t("makeNestedTupleFlat");
  std::vector<int> expectedValues({0, 1, 1, 2, 3, 1, 2, 9});
  std::tuple<std::tuple<int, std::tuple<std::tuple<int, std::tuple<double, float>>, float>>,
             std::tuple<int, std::tuple<double, float>>>
      a({0, {{1, {1.0, 2.0}}, 3.0}}, {1, {2.0, 9.0}});
  const std::tuple<std::tuple<int, std::tuple<std::tuple<int, std::tuple<double, float>>, float>>,
                   std::tuple<int, std::tuple<double, float>>>
      y({0, {{1, {1.0, 2.0}}, 3.0}}, {1, {2.0, 9.0}});
  std::tuple<int, int, double, float, float, int, double, float> aFlat;

  static_assert(std::is_same_v<decltype(aFlat), decltype(Ikarus::Std::makeNestedTupleFlat(a))>);
  static_assert(std::is_same_v<std::tuple<const int&, const int&, const double&, const float&, const float&, const int&,
                                          const double&, const float&>,
                               decltype(Ikarus::Std::makeNestedTupleFlatAndStoreReferences(y))>);
  static_assert(std::is_same_v<std::tuple<int&, int&, double&, float&, float&, int&, double&, float&>,
                               decltype(Ikarus::Std::makeNestedTupleFlatAndStoreReferences(a))>);

  auto reducedTuple            = Ikarus::Std::makeNestedTupleFlatAndStoreReferences(a);
  const auto reducedTupleConst = Ikarus::Std::makeNestedTupleFlatAndStoreReferences(y);

  Dune::Hybrid::forEach(Dune::Hybrid::integralRange(Dune::index_constant<std::tuple_size_v<decltype(a)>>()),
                        [&](const auto i) { t.check(expectedValues[i] == std::get<i>(reducedTuple)); });

  Dune::Hybrid::forEach(Dune::Hybrid::integralRange(Dune::index_constant<std::tuple_size_v<decltype(a)>>()),
                        [&](const auto i) { t.check(expectedValues[i] == std::get<i>(reducedTupleConst)); });
  return t;
}

// move to compile only test
// auto testRebind() {
//   std::vector<int> a;
//   using namespace Ikarus::Std;
//   using reboundVector = Rebind<std::vector<int>, double>::other;
//   reboundVector b;
//   std::vector<double> c;
//   static_assert(!areTypesEqual(a, b));
//   static_assert(areTypesEqual(b, c));
//
//   std::array<int, 4> aA;
//
//   using reboundArray = Rebind<std::array<int, 4>, double>::other;
//   reboundArray bA;
//   std::array<double, 4> cA;
//
//   static_assert(!areTypesEqual(aA, bA));
//   static_assert(areTypesEqual(bA, cA));
// }

int main(int argc, char** argv) {
  Dune::MPIHelper::instance(argc, argv);
  TestSuite t;

  t.subTest(localFunctionTraitsTest());
  t.subTest(treePathCoeffs());
  t.subTest(makeNestedTupleFlat());

  return t.exit();
}
