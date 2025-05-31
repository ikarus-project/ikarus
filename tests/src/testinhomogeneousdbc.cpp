// SPDX-FileCopyrightText: 2021-2025 The Ikarus Developers ikarus@ibb.uni-stuttgart.de
// SPDX-License-Identifier: LGPL-3.0-or-later

#include <config.h>

#include "testelasticstrip.hh"

#include <dune/common/test/testsuite.hh>

#include <Eigen/Core>

#include <ikarus/finiteelements/fefactory.hh>
#include <ikarus/finiteelements/mechanics/enhancedassumedstrains.hh>
#include <ikarus/finiteelements/mechanics/materials.hh>
#include <ikarus/finiteelements/mixin.hh>
using namespace Ikarus;
using Dune::TestSuite;

int main(int argc, char** argv) {
  Ikarus::init(argc, argv);
  Dune::TestSuite t;

  auto matParameterSVK = toLamesFirstParameterAndShearModulus({.emodul = 100.0, .nu = 0.3});
  Materials::StVenantKirchhoff matSVK(matParameterSVK);

  LamesFirstParameterAndShearModulus matParameterNH = {.lambda = 24.0, .mu = 6.0};
  Materials::NeoHooke matNH(matParameterNH);

  auto reducedMats = Dune::makeTupleVector(planeStrain(matSVK), planeStrain(matNH));

  std::array<std::pair<int, double>, 4> expectedResultsSVK = {
      std::make_pair(6, 1.814746879163122), std::make_pair(6, 1.850732157345016), std::make_pair(10, 1.817539174273711),
      std::make_pair(11, 1.8408528524451402)};

  std::array<std::pair<int, double>, 4> expectedResultsNH = {
      std::make_pair(7, 2.207111977583243), std::make_pair(7, 2.19445187105487), std::make_pair(11, 2.2070109913128926),
      std::make_pair(10, 2.2078938377725192)};

  std::array<decltype(expectedResultsSVK), 2> expectedResults = {expectedResultsSVK, expectedResultsNH};

  auto testRange = Dune::Hybrid::integralRange(std::integral_constant<int, 0>(), std::integral_constant<int, 2>());

  auto testFunctor = [&]<bool useArcLength = false>(DBCOption dbcOption) {
    Dune::Hybrid::forEach(testRange, [&](auto i) {
      t.subTest(
          elasticStripTest<useArcLength>(dbcOption, reducedMats[i], skills(), 1, expectedResults[i][0], 1, true, true));
      t.subTest(
          elasticStripTest<useArcLength>(dbcOption, reducedMats[i], skills(), 1, expectedResults[i][1], 2, true, true));
      if (dbcOption == DBCOption::Reduced) {
        t.checkThrow<Dune::NotImplemented>(
            [&]() {
              elasticStripTest<useArcLength>(dbcOption, reducedMats[i], skills(eas<EAS::DisplacementGradient>(4)), 2,
                                             expectedResults[i][2]);
            },
            "elasticStripTest with EAS should throw a Dune::NotImplemented for DBCOption::Reduced.");
      } else {
        t.subTest(elasticStripTest<useArcLength>(dbcOption, reducedMats[i], skills(eas<EAS::DisplacementGradient>(4)),
                                                 2, expectedResults[i][2], 1, true, true));
        t.subTest(elasticStripTest<useArcLength>(dbcOption, reducedMats[i],
                                                 skills(eas<EAS::DisplacementGradientTransposed>(4)), 2,
                                                 expectedResults[i][3], 1, true, true));
      }
    });
  };

  testFunctor(DBCOption::Full);
  testFunctor(DBCOption::Reduced);

  testFunctor.operator()<true>(DBCOption::Full);
  testFunctor.operator()<true>(DBCOption::Reduced);

  return t.exit();
}
