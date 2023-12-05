// SPDX-FileCopyrightText: 2022 The Ikarus Developers mueller@ibb.uni-stuttgart.de
// SPDX-License-Identifier: LGPL-3.0-or-later

#include <config.h>

#include "testeas.hh"
#include "testfeelement.hh"

#include <dune/common/test/testsuite.hh>
#include <dune/functions/functionspacebases/basistags.hh>
#include <dune/functions/functionspacebases/lagrangebasis.hh>
#include <dune/functions/functionspacebases/powerbasis.hh>

#include <ikarus/finiteelements/mechanics/enhancedassumedstrains.hh>
#include <ikarus/finiteelements/mechanics/linearelastic.hh>
#include <ikarus/utils/init.hh>

using Dune::TestSuite;

template <typename Basis>
using EASElement = Ikarus::EnhancedAssumedStrains<Ikarus::LinearElastic<Basis>>;

int main(int argc, char** argv) {
  Ikarus::init(argc, argv);
  TestSuite t("EAS");

  using namespace Dune::Functions::BasisFactory;
  auto firstOrderLagrangePrePower2Basis = power<2>(lagrange<1>(), FlatInterleaved());
  auto firstOrderLagrangePrePower3Basis = power<3>(lagrange<1>(), FlatInterleaved());
  constexpr auto randomlyDistorted      = CornerDistortionFlag::randomlyDistorted;
  constexpr auto unDistorted            = CornerDistortionFlag::unDistorted;

  t.subTest(testFEElement<EASElement, 2>(firstOrderLagrangePrePower2Basis, "EAS", randomlyDistorted,
                                         Dune::GeometryTypes::cube(2), checkJacobianFunctor));
  t.subTest(testFEElement<EASElement, 3>(firstOrderLagrangePrePower3Basis, "EAS", randomlyDistorted,
                                         Dune::GeometryTypes::cube(3), checkJacobianFunctor));
  t.subTest(testFEElement<EASElement, 2>(firstOrderLagrangePrePower2Basis, "EAS", unDistorted,
                                         Dune::GeometryTypes::cube(2), checkLinearStressFunctor,
                                         checkResultFunctionFunctor));

  return t.exit();
}
