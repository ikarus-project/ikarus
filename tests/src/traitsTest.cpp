//
// Created by Alex on 21.04.2021.
//

#include <gmock/gmock.h>

#include "testHelpers.hh"

#include <ikarus/localFunctions/meta.hh>
#include <ikarus/utils/traits.hh>
TEST(TraitsTest, LocalFunctionTest) {
  using namespace Ikarus::DerivativeDirections;
  auto wrt1 = Ikarus::wrt(spatial(0), coeff(7));

  auto counter = countDerivativesType<decltype(wrt1)>();

  EXPECT_EQ(counter.singleCoeffDerivs, 1);
  EXPECT_EQ(counter.twoCoeffDerivs, 0);
  EXPECT_EQ(counter.spatialDerivs, 1);
  EXPECT_EQ(counter.spatialAllCounter, 0);

  auto wrt2 = Ikarus::wrt(spatialAll, coeff(7, 1));

  counter = countDerivativesType<decltype(wrt2)>();

  EXPECT_EQ(counter.singleCoeffDerivs, 0);
  EXPECT_EQ(counter.twoCoeffDerivs, 1);
  EXPECT_EQ(counter.spatialDerivs, 0);
  EXPECT_EQ(counter.spatialAllCounter, 1);
}


TEST(TraitsTest, TreePathCoeffs) {
  using namespace Ikarus::DerivativeDirections;
  using namespace Dune::Indices;
  auto coeffTwo = coeff(_0, 7,_1,9);

  EXPECT_EQ(coeffTwo.index[_0][0],0);
  EXPECT_EQ(coeffTwo.index[_0][1],7);
  EXPECT_EQ(coeffTwo.index[_1][0],1);
  EXPECT_EQ(coeffTwo.index[_1][1],9);

    auto coeffSingle = coeff(_4, 1);

  EXPECT_EQ(coeffSingle.index[_0][0],4);
  EXPECT_EQ(coeffSingle.index[_0][1],1);

}



TEST(TraitsTest, testTupleFilter) {
  using namespace Ikarus::DerivativeDirections;
  using namespace Dune::Indices;
  auto tup = std::make_tuple(_0,_1,_3,Ikarus::arithmetic);
  auto filteredTupExpected = std::make_tuple(_0,_1,_3);
  constexpr auto predicate = []<typename Type>(Type ){return Type::value!=Ikarus::arithmetic;};
  auto filteredTup = Ikarus::Std::filter(tup,predicate);

  static_assert(std::is_same_v<decltype(filteredTup),decltype(filteredTupExpected)>);

  auto tup2 = std::make_tuple(Ikarus::arithmetic,_0,Ikarus::arithmetic,_1,_3,Ikarus::arithmetic);
  auto filteredTup2 = Ikarus::Std::filter(tup2,predicate);
  static_assert(std::is_same_v<decltype(filteredTup2),decltype(filteredTupExpected)>);

}


TEST(TraitsTest, makeNestedTupleFlat) {
  std::vector<int> expectedValues({0,1,1,2,3,1,2,9});
  std::tuple<std::tuple<int,std::tuple<std::tuple<int,std::tuple<double,float>>,float>>,std::tuple<int,std::tuple<double,float>>> a({0,{{1,{1.0,2.0}},3.0}},{1,{2.0,9.0}});
  const std::tuple<std::tuple<int,std::tuple<std::tuple<int,std::tuple<double,float>>,float>>,std::tuple<int,std::tuple<double,float>>> y({0,{{1,{1.0,2.0}},3.0}},{1,{2.0,9.0}});
  std::tuple<int, int, double, float, float, int, double, float> aFlat;

  static_assert(std::is_same_v<decltype(aFlat),decltype(Ikarus::Std::makeNestedTupleFlat(a))>);
  static_assert(std::is_same_v<std::tuple<const int&, const int&, const double&, const float&, const float&, const int&, const double&, const float&>,decltype(Ikarus::Std::makeNestedTupleFlatAndStoreReferences(y))>);
  static_assert(std::is_same_v<std::tuple< int&,  int&,  double&,  float&,  float&,  int&,  double&,  float&>,decltype(Ikarus::Std::makeNestedTupleFlatAndStoreReferences(a))>);

  auto reducedTuple = Ikarus::Std::makeNestedTupleFlatAndStoreReferences(a);
  const auto reducedTupleConst = Ikarus::Std::makeNestedTupleFlatAndStoreReferences(y);

  Dune::Hybrid::forEach(Dune::Hybrid::integralRange(Dune::index_constant<std::tuple_size_v<decltype(a)>>()), [&](const auto i) {
    EXPECT_EQ(std::get<i>(reducedTuple),expectedValues[i]);
  });

  Dune::Hybrid::forEach(Dune::Hybrid::integralRange(Dune::index_constant<std::tuple_size_v<decltype(a)>>()), [&](const auto i) {
    EXPECT_EQ(std::get<i>(reducedTupleConst),expectedValues[i]);
  });
}


TEST(TraitsTest, testRebind) {
  std::vector<int> a;
  using namespace Ikarus::Std;
  using reboundVector = Rebind<std::vector<int>,double>::other ;
  reboundVector b;
  std::vector<double> c;
  static_assert(!areTypesEqual(a,b));
  static_assert(areTypesEqual(b,c));

  std::array<int,4> aA;

  using reboundArray = Rebind<std::array<int,4>,double>::other ;
  reboundArray bA;
  std::array<double,4> cA;

  static_assert(!areTypesEqual(aA,bA));
  static_assert(areTypesEqual(bA,cA));
}