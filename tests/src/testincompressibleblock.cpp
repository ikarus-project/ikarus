// SPDX-FileCopyrightText: 2021-2025 The Ikarus Developers ikarus@ibb.uni-stuttgart.de
// SPDX-License-Identifier: LGPL-3.0-or-later

#include <config.h>

#include "testincompressibleblock.hh"

#include <dune/common/test/testsuite.hh>

#include <Eigen/Core>

#include <ikarus/finiteelements/fefactory.hh>
#include <ikarus/finiteelements/mechanics/assumedstress.hh>
#include <ikarus/finiteelements/mechanics/displacementpressure.hh>
#include <ikarus/finiteelements/mechanics/materials.hh>
#include <ikarus/finiteelements/mixin.hh>

using namespace Ikarus;
using Dune::TestSuite;

auto incompressibelBlockResults(int nele) {
  if (nele == 2)
    return std::make_pair(17, 45.08467250483267);
  else if (nele == 4)
    return std::make_pair(19, 41.32367771133791);
  else if (nele == 8)
    return std::make_pair(20, 40.56242425245286);
  else
    DUNE_THROW(Dune::NotImplemented, "No TestResult for the given nele");
}

int main(int argc, char** argv) {
  Ikarus::init(argc, argv);
  Dune::TestSuite t;

  double mu     = 1.61148;
  double lambda = 499.92568;
  auto kappa    = lambda + 2 * mu / 3;

  auto matDEV = Materials::makeOgden<1, Ikarus::PrincipalStretchTags::deviatoric>({mu}, {2.0});
  auto matVOL = Materials::makeMaterialLawFromPenaltyFunction(Materials::PVF1());

  using namespace Dune::Functions::BasisFactory;

  auto sk = displacementPressure(matDEV, matVOL, kappa);

  for (auto i : {2, 4, 8})
    t.subTest(incompressibelBlockTest(sk, incompressibelBlockResults(i), i));

  return t.exit();
}
