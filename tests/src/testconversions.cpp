// SPDX-FileCopyrightText: 2021-2025 The Ikarus Developers mueller@ibb.uni-stuttgart.de
// SPDX-License-Identifier: LGPL-3.0-or-later

#include <config.h>

#include <dune/common/test/testsuite.hh>

#include <ikarus/finiteelements/mechanics/materials/strainconversions.hh>
#include <ikarus/finiteelements/mechanics/materials/stressconversions.hh>
#include <ikarus/utils/init.hh>

using Dune::TestSuite;
using namespace Ikarus;

auto testToAndFromStressConversions() {
  TestSuite t("To and From Stress Transformations");

  auto F = Eigen::Matrix3d{
      {1.0, 0.0, 0.5},
      {0.0, 0.2, 0.5},
      {0.2, 0.5, 1.0}
  };
  auto S = Eigen::Matrix3d{
      { 0.600872, -0.179083, 0},
      {-0.179083,  0.859121, 0},
      {        0,         0, 1}
  };

  // constexpr auto tags = Dune::makeTupleVector(StressTags::PK2, StressTags::PK1, StressTags::Kirchhoff, StressTags::Cauchy);
  // Dune::Hybrid::forEach(tags, [&](auto tag1) {
  //   Dune::Hybrid::forEach(tags, [&](auto tag2) {
  //     auto transformedStress = transformStress<tag1, tag2>(S, F);
  //   });
  // });
  auto stressTagRange = Dune::Hybrid::integralRange(std::integral_constant<int, 2>(), std::integral_constant<int, 4>());
  Dune::Hybrid::forEach(stressTagRange, [&](auto i) {
    constexpr StressTags tag = StressTags(int(i));
    transformStress<StressTags::PK2, tag>(S, F);
    // auto roundTripStress = transformStress<tag, StressTags::PK2>(convertedStress, F);
    // t.check(roundTripStress.isApprox(S, 1e-6), "Round-trip stress conversion failed");
  });


  return t;
}

int main(int argc, char** argv) {
  Ikarus::init(argc, argv);
  TestSuite t;

  t.subTest(testToAndFromStressConversions());

  return t.exit();
}
