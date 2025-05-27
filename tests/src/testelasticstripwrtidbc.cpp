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
  LamesFirstParameterAndShearModulus matParameter = {.lambda = 24000.0, .mu = 6000.0};
  Materials::NeoHooke matNH(matParameter);
  auto reducedMat = planeStrain(matNH);

  auto testFunctor = [&](DBCOption dbcOption) {
    t.subTest(elasticStripTest(dbcOption, reducedMat, skills(), 1, std::make_pair(8, 2.207), 1, true, true));
    t.subTest(elasticStripTest(dbcOption, reducedMat, skills(), 1, std::make_pair(8, 2.194), 2, true, true));
    t.subTest(elasticStripTest(dbcOption, reducedMat, skills(eas<EAS::GreenLagrangeStrain>(4)), 2,
                               std::make_pair(12, 2.207), 1, true, true));
    t.subTest(elasticStripTest(dbcOption, reducedMat, skills(eas<EAS::DisplacementGradient>(4)), 2,
                               std::make_pair(12, 2.207), 1, true, true));
    t.subTest(elasticStripTest(dbcOption, reducedMat, skills(eas<EAS::DisplacementGradientTransposed>(4)), 2,
                               std::make_pair(12, 2.207), 1, true, true));
  };

  testFunctor(DBCOption::Full);
  testFunctor(DBCOption::Reduced);

  return t.exit();
}
