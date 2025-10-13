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

auto incompressibelBlockResults(int pD, bool continuous, int nele) {
  if (pD == 1) {
    if (nele == 2)
      return std::make_pair(17, 45.08467250483267);
    else if (nele == 4)
      return std::make_pair(19, 41.32367771133791);
    else if (nele == 8)
      return std::make_pair(20, 40.56242425245286);
    else
      DUNE_THROW(Dune::NotImplemented, "No TestResult for the given nele");
  } else if (pD == 2) {
    if (continuous) {
      if (nele == 2)
        return std::make_pair(19, 39.97523916033919);
      else if (nele == 4)
        return std::make_pair(19, 39.92020798401635);
      else
        DUNE_THROW(Dune::NotImplemented, "No TestResult for the given nele");
    } else {
      if (nele == 2)
        return std::make_pair(22, 40.786978335135046);
      else if (nele == 4)
        return std::make_pair(20, 40.14084241096578);
      else
        DUNE_THROW(Dune::NotImplemented, "No TestResult for the given nele");
    }
  } else {
    DUNE_THROW(Dune::NotImplemented, "No TestResult for the given polynomial degree");
  }
}

int main(int argc, char** argv) {
  Ikarus::init(argc, argv);
  Dune::TestSuite t;

  double mu     = 1.61148;
  double lambda = 499.92568;
  auto kappa    = lambda + 2 * mu / 3;

  auto mat = Materials::makeOgden<1, Ikarus::PrincipalStretchTags::deviatoric>({mu}, {2.0}, kappa, Materials::VF12());

  using namespace Dune::Functions::BasisFactory;

  auto sk = displacementPressure(mat);

  for (auto i : {2, 4, 8})
    t.subTest(incompressibelBlockTest<1, 0, false>(sk, incompressibelBlockResults(1, false, i), i));

  // The H2P1 element can have a continuous or discontinuous field for the pressure
  for (auto i : {2, 4}) {
    t.subTest(incompressibelBlockTest<2, 1, true>(sk, incompressibelBlockResults(2, true, i), i));
    t.subTest(incompressibelBlockTest<2, 1, false>(sk, incompressibelBlockResults(2, false, i), i));
  }
  return t.exit();
}
