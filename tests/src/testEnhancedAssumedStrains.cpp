// SPDX-FileCopyrightText: 2022 The Ikarus Developers mueller@ibb.uni-stuttgart.de
// SPDX-License-Identifier: LGPL-3.0-or-later

#include <config.h>

#include "testEAS.hh"
#include "testFEElement.hh"

#include <dune/common/test/testsuite.hh>
#include <dune/functions/functionspacebases/basistags.hh>
#include <dune/functions/functionspacebases/lagrangebasis.hh>
#include <dune/functions/functionspacebases/powerbasis.hh>

#include <ikarus/finiteElements/mechanics/enhancedAssumedStrains.hh>
#include <ikarus/finiteElements/mechanics/linearElastic.hh>
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
  const auto randomlyDistorted          = CornerDistortionFlag::randomlyDistorted;
  const auto unDistorted                = CornerDistortionFlag::unDistorted;

  t.subTest(
      testFEElement<EASElement, 2>(firstOrderLagrangePrePower2Basis, "EAS", randomlyDistorted, checkJacobianFunctor));
  t.subTest(
      testFEElement<EASElement, 3>(firstOrderLagrangePrePower3Basis, "EAS", randomlyDistorted, checkJacobianFunctor));
  t.subTest(
      testFEElement<EASElement, 2>(firstOrderLagrangePrePower2Basis, "EAS", unDistorted, checkCauchyStressFunctor));

  return t.exit();
}
