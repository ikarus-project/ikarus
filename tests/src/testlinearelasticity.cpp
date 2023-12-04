// SPDX-FileCopyrightText: 2022 The Ikarus Developers mueller@ibb.uni-stuttgart.de
// SPDX-License-Identifier: LGPL-3.0-or-later

#include <config.h>

#include "testfeelement.hh"

#include <dune/common/test/testsuite.hh>
#include <dune/functions/functionspacebases/basistags.hh>
#include <dune/functions/functionspacebases/lagrangebasis.hh>
#include <dune/functions/functionspacebases/powerbasis.hh>

#include <ikarus/finiteelements/mechanics/linearelastic.hh>
#include <ikarus/utils/init.hh>

using Dune::TestSuite;

template <typename Basis>
using LinearElasticElement = Ikarus::LinearElastic<Basis>;

int main(int argc, char** argv) {
  Ikarus::init(argc, argv);
  TestSuite t("LinearElasticity");

  using namespace Dune::Functions::BasisFactory;
  auto firstOrderLagrangePrePower2Basis         = power<2>(lagrange<1>(), FlatInterleaved());
  auto secondOrderLagrangePrePower2Basis        = power<2>(lagrange<2>(), FlatInterleaved());
  auto firstOrderLagrangePrePower3Basis         = power<3>(lagrange<1>(), FlatInterleaved());
  auto secondOrderLagrangePrePower3Basis        = power<3>(lagrange<2>(), FlatInterleaved());
  auto secondOrderLagrangePrePower3BasisBlocked = power<3>(lagrange<2>());
  constexpr auto randomlyDistorted              = CornerDistortionFlag::randomlyDistorted;
  constexpr auto unDistorted                    = CornerDistortionFlag::unDistorted;
  constexpr auto cubeGeometry                   = ElementGeometryFlag::cube;
  constexpr auto simplexGeometry                = ElementGeometryFlag::simplex;

  // Test cube 2D
  t.subTest(testFEElement<LinearElasticElement, 2>(firstOrderLagrangePrePower2Basis, "LinearElastic", randomlyDistorted,
                                                   cubeGeometry, checkGradientFunctor, checkHessianFunctor,
                                                   checkJacobianFunctor, checkFEByAutoDiffFunctor));
  t.subTest(testFEElement<LinearElasticElement, 2>(
      secondOrderLagrangePrePower2Basis, "LinearElastic", randomlyDistorted, cubeGeometry, checkGradientFunctor,
      checkHessianFunctor, checkJacobianFunctor, checkFEByAutoDiffFunctor));

  t.subTest(testFEElement<LinearElasticElement, 2>(firstOrderLagrangePrePower2Basis, "LinearElastic", unDistorted,
                                                   cubeGeometry, checkLinearStressFunctor, checkResultFunctionFunctor));

  // Test simplex 2D
  t.subTest(testFEElement<LinearElasticElement, 2>(firstOrderLagrangePrePower2Basis, "LinearElastic", randomlyDistorted,
                                                   simplexGeometry, checkGradientFunctor, checkHessianFunctor,
                                                   checkJacobianFunctor, checkFEByAutoDiffFunctor));
  t.subTest(testFEElement<LinearElasticElement, 2>(
      secondOrderLagrangePrePower2Basis, "LinearElastic", randomlyDistorted, simplexGeometry, checkGradientFunctor,
      checkHessianFunctor, checkJacobianFunctor, checkFEByAutoDiffFunctor));

  t.subTest(testFEElement<LinearElasticElement, 2>(firstOrderLagrangePrePower2Basis, "LinearElastic", unDistorted,
                                                   simplexGeometry, checkLinearStressFunctor));

  // Test cube 3D
  t.subTest(testFEElement<LinearElasticElement, 3>(firstOrderLagrangePrePower3Basis, "LinearElastic", randomlyDistorted,
                                                   cubeGeometry, checkGradientFunctor, checkHessianFunctor,
                                                   checkJacobianFunctor, checkFEByAutoDiffFunctor));
  t.subTest(testFEElement<LinearElasticElement, 3>(
      secondOrderLagrangePrePower3Basis, "LinearElastic", randomlyDistorted, cubeGeometry, checkGradientFunctor,
      checkHessianFunctor, checkJacobianFunctor, checkFEByAutoDiffFunctor));
  t.subTest(testFEElement<LinearElasticElement, 3>(
      secondOrderLagrangePrePower3BasisBlocked, "LinearElastic", randomlyDistorted, cubeGeometry, checkGradientFunctor,
      checkHessianFunctor, checkJacobianFunctor, checkFEByAutoDiffFunctor));

  // Test simplex 3D
  t.subTest(testFEElement<LinearElasticElement, 3>(firstOrderLagrangePrePower3Basis, "LinearElastic", randomlyDistorted,
                                                   simplexGeometry, checkGradientFunctor, checkHessianFunctor,
                                                   checkJacobianFunctor, checkFEByAutoDiffFunctor));

  t.subTest(testFEElement<LinearElasticElement, 3>(firstOrderLagrangePrePower3Basis, "LinearElastic", unDistorted,
                                                   simplexGeometry, checkLinearStressFunctor,
                                                   checkResultFunctionFunctor));
  return t.exit();
}
