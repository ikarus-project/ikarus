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
  LamesFirstParameterAndShearModulus matParameterNH = {.lambda = 24000.0, .mu = 6000.0};
  auto matParameterSVK = toLamesFirstParameterAndShearModulus({.emodul = 100.0, .nu = 0.3});
  Materials::StVenantKirchhoff matSVK(matParameterSVK);
  Materials::NeoHooke matNH(matParameterNH);
  auto reducedMats = Dune::makeTupleVector(planeStrain(matSVK), planeStrain(matNH));

  std::array<std::pair<int, double>, 5> expectedResultsSVK = {
      std::make_pair(6, 1.814746879163122), std::make_pair(6, 1.850732157345016), std::make_pair(12, 1.818440236431206),
      std::make_pair(18, 1.817539174254031), std::make_pair(6, 1.814746879163122)};

  std::array<std::pair<int, double>, 5> expectedResultsNH = {
      std::make_pair(8, 2.207111977583243), std::make_pair(8, 2.19445187105487), std::make_pair(12, 2.207),
      std::make_pair(12, 2.207), std::make_pair(12, 2.207)};

  std::array<decltype(expectedResultsSVK), 2> expectedResults = {expectedResultsSVK, expectedResultsNH};

  auto testRange = Dune::Hybrid::integralRange(std::integral_constant<int, 0>(), std::integral_constant<int, 2>());

  auto testFunctor = [&](DBCOption dbcOption) {
    Dune::Hybrid::forEach(testRange, [&](auto i) {
      t.subTest(elasticStripTest(dbcOption, reducedMats[i], skills(), 1, expectedResults[i][0], 1, true, true));
      t.subTest(elasticStripTest(dbcOption, reducedMats[i], skills(), 1, expectedResults[i][1], 2, true, true));
      t.subTest(elasticStripTest(dbcOption, reducedMats[i], skills(eas<EAS::GreenLagrangeStrain>(4)), 2,
                                 expectedResults[i][2], 1, true, true));
      t.subTest(elasticStripTest(dbcOption, reducedMats[i], skills(eas<EAS::DisplacementGradient>(4)), 2,
                                 expectedResults[i][3], 1, true, true));
      t.subTest(elasticStripTest(dbcOption, reducedMats[i], skills(eas<EAS::DisplacementGradientTransposed>(4)), 2,
                                 expectedResults[i][4], 1, true, true));
    });
  };

  testFunctor(DBCOption::Full);
  testFunctor(DBCOption::Reduced);

  return t.exit();
}
