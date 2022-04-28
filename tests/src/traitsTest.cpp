//
// Created by Alex on 21.04.2021.
//

#include <gmock/gmock.h>

#include "testHelpers.hh"

#include <ikarus/localFunctions/meta.hh>
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