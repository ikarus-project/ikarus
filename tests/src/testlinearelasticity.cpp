// SPDX-FileCopyrightText: 2021-2024 The Ikarus Developers mueller@ibb.uni-stuttgart.de
// SPDX-License-Identifier: LGPL-3.0-or-later

#include <config.h>

#include "resultcollection.hh"
#include "testfeelement.hh"

#include <dune/common/test/testsuite.hh>
#include <dune/functions/functionspacebases/basistags.hh>
#include <dune/functions/functionspacebases/lagrangebasis.hh>
#include <dune/functions/functionspacebases/powerbasis.hh>

#include <ikarus/finiteelements/feresulttypes.hh>
#include <ikarus/finiteelements/mechanics/linearelastic.hh>
#include <ikarus/io/resultevaluators.hh>
#include <ikarus/utils/init.hh>

using Dune::TestSuite;

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

  // Test cube 2D
  auto linearElasticFunc = [](const Ikarus::YoungsModulusAndPoissonsRatio& mat) { return Ikarus::linearElastic(mat); };
  t.subTest(testFEElement(firstOrderLagrangePrePower2Basis, "LinearElastic", randomlyDistorted,
                          Dune::ReferenceElements<double, 2>::cube(), linearElasticFunc, Ikarus::skills(),
                          checkGradientFunctor, checkHessianFunctor, checkJacobianFunctor, checkFEByAutoDiffFunctor));
  t.subTest(testFEElement(secondOrderLagrangePrePower2Basis, "LinearElastic", randomlyDistorted,
                          Dune::ReferenceElements<double, 2>::cube(), linearElasticFunc, Ikarus::skills(),
                          checkGradientFunctor, checkHessianFunctor, checkJacobianFunctor, checkFEByAutoDiffFunctor));

  t.subTest(testFEElement(
      firstOrderLagrangePrePower2Basis, "LinearElastic", unDistorted, Dune::ReferenceElements<double, 2>::cube(),
      linearElasticFunc, Ikarus::skills(),
      checkCalculateAtFunctorFactory<Ikarus::ResultTypes::linearStress>(linearStressResultsOfSquare),
      checkCalculateAtFunctorFactory<Ikarus::ResultTypes::linearStress, false>(linearStressResultsOfSquare),
      checkResultFunctionFunctorFactory<Ikarus::ResultTypes::linearStress>(linearStressResultsOfSquare),
      checkResultFunctionFunctorFactory<Ikarus::ResultTypes::linearStress, Ikarus::ResultEvaluators::VonMises>(
          linearVonMisesResultsOfSquare),
      checkResultFunctionFunctorFactory<Ikarus::ResultTypes::linearStress,
                                        Ikarus::ResultEvaluators::PrincipalStress<2>>(
          linearPrincipalStressResultsOfSquare)));

  // Test simplex 2D
  t.subTest(testFEElement(firstOrderLagrangePrePower2Basis, "LinearElastic", randomlyDistorted,
                          Dune::ReferenceElements<double, 2>::simplex(), linearElasticFunc, Ikarus::skills(),
                          checkGradientFunctor, checkHessianFunctor, checkJacobianFunctor, checkFEByAutoDiffFunctor));
  t.subTest(testFEElement(secondOrderLagrangePrePower2Basis, "LinearElastic", randomlyDistorted,
                          Dune::ReferenceElements<double, 2>::simplex(), linearElasticFunc, Ikarus::skills(),
                          checkGradientFunctor, checkHessianFunctor, checkJacobianFunctor, checkFEByAutoDiffFunctor));

  // Test simplex 2D
  t.subTest(testFEElement(firstOrderLagrangePrePower2Basis, "LinearElastic", randomlyDistorted,
                          Dune::ReferenceElements<double, 2>::simplex(), linearElasticFunc, Ikarus::skills(),
                          checkGradientFunctor, checkHessianFunctor, checkJacobianFunctor, checkFEByAutoDiffFunctor));
  t.subTest(testFEElement(secondOrderLagrangePrePower2Basis, "LinearElastic", randomlyDistorted,
                          Dune::ReferenceElements<double, 2>::simplex(), linearElasticFunc, Ikarus::skills(),
                          checkGradientFunctor, checkHessianFunctor, checkJacobianFunctor, checkFEByAutoDiffFunctor));

  t.subTest(
      testFEElement(firstOrderLagrangePrePower2Basis, "LinearElastic", unDistorted,
                    Dune::ReferenceElements<double, 2>::simplex(), linearElasticFunc, Ikarus::skills(),
                    checkCalculateAtFunctorFactory<Ikarus::ResultTypes::linearStress>(linearStressResultsOfTriangle)));

  // Test cube 3D
  t.subTest(testFEElement(firstOrderLagrangePrePower3Basis, "LinearElastic", randomlyDistorted,
                          Dune::ReferenceElements<double, 3>::cube(), linearElasticFunc, Ikarus::skills(),
                          checkGradientFunctor, checkHessianFunctor, checkJacobianFunctor, checkFEByAutoDiffFunctor));
  t.subTest(testFEElement(secondOrderLagrangePrePower3Basis, "LinearElastic", randomlyDistorted,
                          Dune::ReferenceElements<double, 3>::cube(), linearElasticFunc, Ikarus::skills(),
                          checkGradientFunctor, checkHessianFunctor, checkJacobianFunctor, checkFEByAutoDiffFunctor));
  t.subTest(testFEElement(secondOrderLagrangePrePower3BasisBlocked, "LinearElastic", randomlyDistorted,
                          Dune::ReferenceElements<double, 3>::cube(), linearElasticFunc, Ikarus::skills(),
                          checkGradientFunctor, checkHessianFunctor, checkJacobianFunctor, checkFEByAutoDiffFunctor));

  t.subTest(testFEElement(
      firstOrderLagrangePrePower3Basis, "LinearElastic", unDistorted, Dune::ReferenceElements<double, 3>::cube(),
      linearElasticFunc, Ikarus::skills(),
      checkCalculateAtFunctorFactory<Ikarus::ResultTypes::linearStress>(linearStressResultsOfCube),
      checkResultFunctionFunctorFactory<Ikarus::ResultTypes::linearStress>(linearStressResultsOfCube),
      checkResultFunctionFunctorFactory<Ikarus::ResultTypes::linearStress, Ikarus::ResultEvaluators::VonMises>(
          linearVonMisesResultsOfCube),
      checkResultFunctionFunctorFactory<Ikarus::ResultTypes::linearStress,
                                        Ikarus::ResultEvaluators::PrincipalStress<3>>(
          linearPrincipalStressResultsOfCube)));

  // Test simplex 3D
  t.subTest(testFEElement(firstOrderLagrangePrePower3Basis, "LinearElastic", randomlyDistorted,
                          Dune::ReferenceElements<double, 3>::simplex(), linearElasticFunc, Ikarus::skills(),
                          checkGradientFunctor, checkHessianFunctor, checkJacobianFunctor, checkFEByAutoDiffFunctor));

  t.subTest(testFEElement(
      firstOrderLagrangePrePower3Basis, "LinearElastic", unDistorted, Dune::ReferenceElements<double, 3>::simplex(),
      linearElasticFunc, Ikarus::skills(),
      checkCalculateAtFunctorFactory<Ikarus::ResultTypes::linearStress>(linearStressResultsOfTetrahedron),
      checkResultFunctionFunctorFactory<Ikarus::ResultTypes::linearStress>(linearStressResultsOfTetrahedron)));
  return t.exit();
}
