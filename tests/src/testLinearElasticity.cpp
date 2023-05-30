// SPDX-FileCopyrightText: 2022 The Ikarus Developers mueller@ibb.uni-stuttgart.de
// SPDX-License-Identifier: LGPL-3.0-or-later

#include <config.h>

#include "testFEElement.hh"

#include <dune/common/test/testsuite.hh>
#include <dune/functions/functionspacebases/basistags.hh>
#include <dune/functions/functionspacebases/lagrangebasis.hh>
#include <dune/functions/functionspacebases/powerbasis.hh>

#include <ikarus/finiteElements/mechanics/linearElastic.hh>
#include <ikarus/utils/init.hh>

using Dune::TestSuite;

template <typename Basis>
using LinearElasticElement = Ikarus::LinearElastic<Basis>;

int main(int argc, char** argv) {
  Ikarus::init(argc, argv);
  Dune::TestSuite t("LinearElasticity");

  using namespace Dune::Functions::BasisFactory;
  auto firstOrderLagrangePrePower2Basis         = power<2>(lagrange<1>(), FlatInterleaved());
  auto secondOrderLagrangePrePower2Basis        = power<2>(lagrange<2>(), FlatInterleaved());
  auto firstOrderLagrangePrePower3Basis         = power<3>(lagrange<1>(), FlatInterleaved());
  auto secondOrderLagrangePrePower3Basis        = power<3>(lagrange<2>(), FlatInterleaved());
  auto secondOrderLagrangePrePower3BasisBlocked = power<3>(lagrange<2>());
  const auto randomlyDistorted                  = CornerDistortionFlag::randomlyDistorted;
  const auto unDistorted                        = CornerDistortionFlag::unDistorted;

  t.subTest(testFEElement<LinearElasticElement, 2>(firstOrderLagrangePrePower2Basis, "LinearElastic", randomlyDistorted,
                                                   checkGradientFunctor, checkHessianFunctor, checkJacobianFunctor,
                                                   checkFEByAutoDiffFunctor));
  t.subTest(testFEElement<LinearElasticElement, 2>(secondOrderLagrangePrePower2Basis, "LinearElastic",
                                                   randomlyDistorted, checkGradientFunctor, checkHessianFunctor,
                                                   checkJacobianFunctor, checkFEByAutoDiffFunctor));
  t.subTest(testFEElement<LinearElasticElement, 3>(firstOrderLagrangePrePower3Basis, "LinearElastic", randomlyDistorted,
                                                   checkGradientFunctor, checkHessianFunctor, checkJacobianFunctor,
                                                   checkFEByAutoDiffFunctor));
  t.subTest(testFEElement<LinearElasticElement, 3>(secondOrderLagrangePrePower3Basis, "LinearElastic",
                                                   randomlyDistorted, checkGradientFunctor, checkHessianFunctor,
                                                   checkJacobianFunctor, checkFEByAutoDiffFunctor));
  t.subTest(testFEElement<LinearElasticElement, 3>(secondOrderLagrangePrePower3BasisBlocked, "LinearElastic",
                                                   randomlyDistorted, checkGradientFunctor, checkHessianFunctor,
                                                   checkJacobianFunctor, checkFEByAutoDiffFunctor));
  t.subTest(testFEElement<LinearElasticElement, 2>(firstOrderLagrangePrePower2Basis, "LinearElastic", unDistorted,
                                                   checkCauchyStressFunctor));
  return t.exit();
}
