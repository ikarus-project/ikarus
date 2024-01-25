// SPDX-FileCopyrightText: 2021-2024 The Ikarus Developers mueller@ibb.uni-stuttgart.de
// SPDX-License-Identifier: LGPL-3.0-or-later

#include <config.h>

#include "resultcollection.hh"
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

  // Test cube 2D
  t.subTest(testFEElement<LinearElasticElement>(firstOrderLagrangePrePower2Basis, "LinearElastic", randomlyDistorted,
                                                Dune::ReferenceElements<double, 2>::cube(), checkGradientFunctor,
                                                checkHessianFunctor, checkJacobianFunctor, checkFEByAutoDiffFunctor));
  t.subTest(testFEElement<LinearElasticElement>(secondOrderLagrangePrePower2Basis, "LinearElastic", randomlyDistorted,
                                                Dune::ReferenceElements<double, 2>::cube(), checkGradientFunctor,
                                                checkHessianFunctor, checkJacobianFunctor, checkFEByAutoDiffFunctor));

  t.subTest(testFEElement<LinearElasticElement>(
      firstOrderLagrangePrePower2Basis, "LinearElastic", unDistorted, Dune::ReferenceElements<double, 2>::cube(),
      checkCalculateAtFunctorFactory<Ikarus::ResultType::linearStress>(linearStressResultsOfSquare),
      checkResultFunctionFunctorFactory<Ikarus::ResultType::linearStress>()));

  // Test simplex 2D
  t.subTest(testFEElement<LinearElasticElement>(firstOrderLagrangePrePower2Basis, "LinearElastic", randomlyDistorted,
                                                Dune::ReferenceElements<double, 2>::simplex(), checkGradientFunctor,
                                                checkHessianFunctor, checkJacobianFunctor, checkFEByAutoDiffFunctor));
  t.subTest(testFEElement<LinearElasticElement>(secondOrderLagrangePrePower2Basis, "LinearElastic", randomlyDistorted,
                                                Dune::ReferenceElements<double, 2>::simplex(), checkGradientFunctor,
                                                checkHessianFunctor, checkJacobianFunctor, checkFEByAutoDiffFunctor));

  // Test simplex 2D
  t.subTest(testFEElement<LinearElasticElement>(firstOrderLagrangePrePower2Basis, "LinearElastic", randomlyDistorted,
                                                Dune::ReferenceElements<double, 2>::simplex(), checkGradientFunctor,
                                                checkHessianFunctor, checkJacobianFunctor, checkFEByAutoDiffFunctor));
  t.subTest(testFEElement<LinearElasticElement>(secondOrderLagrangePrePower2Basis, "LinearElastic", randomlyDistorted,
                                                Dune::ReferenceElements<double, 2>::simplex(), checkGradientFunctor,
                                                checkHessianFunctor, checkJacobianFunctor, checkFEByAutoDiffFunctor));

  t.subTest(testFEElement<LinearElasticElement>(
      firstOrderLagrangePrePower2Basis, "LinearElastic", unDistorted, Dune::ReferenceElements<double, 2>::simplex(),
      checkCalculateAtFunctorFactory<Ikarus::ResultType::linearStress>(linearStressResultsOfTriangle)));

  // Test cube 3D
  t.subTest(testFEElement<LinearElasticElement>(firstOrderLagrangePrePower3Basis, "LinearElastic", randomlyDistorted,
                                                Dune::ReferenceElements<double, 3>::cube(), checkGradientFunctor,
                                                checkHessianFunctor, checkJacobianFunctor, checkFEByAutoDiffFunctor));
  t.subTest(testFEElement<LinearElasticElement>(secondOrderLagrangePrePower3Basis, "LinearElastic", randomlyDistorted,
                                                Dune::ReferenceElements<double, 3>::cube(), checkGradientFunctor,
                                                checkHessianFunctor, checkJacobianFunctor, checkFEByAutoDiffFunctor));
  t.subTest(testFEElement<LinearElasticElement>(secondOrderLagrangePrePower3BasisBlocked, "LinearElastic",
                                                randomlyDistorted, Dune::ReferenceElements<double, 3>::cube(),
                                                checkGradientFunctor, checkHessianFunctor, checkJacobianFunctor,
                                                checkFEByAutoDiffFunctor));

  t.subTest(testFEElement<LinearElasticElement>(
      firstOrderLagrangePrePower3Basis, "LinearElastic", unDistorted, Dune::ReferenceElements<double, 3>::cube(),
      checkCalculateAtFunctorFactory<Ikarus::ResultType::linearStress>(linearStressResultsOfCube),
      checkResultFunctionFunctorFactory<Ikarus::ResultType::linearStress>()));

  // Test simplex 3D
  t.subTest(testFEElement<LinearElasticElement>(firstOrderLagrangePrePower3Basis, "LinearElastic", randomlyDistorted,
                                                Dune::ReferenceElements<double, 3>::simplex(), checkGradientFunctor,
                                                checkHessianFunctor, checkJacobianFunctor, checkFEByAutoDiffFunctor));

  t.subTest(testFEElement<LinearElasticElement>(
      firstOrderLagrangePrePower3Basis, "LinearElastic", unDistorted, Dune::ReferenceElements<double, 3>::simplex(),
      checkCalculateAtFunctorFactory<Ikarus::ResultType::linearStress>(linearStressResultsOfTetrahedron),
      checkResultFunctionFunctorFactory<Ikarus::ResultType::linearStress>()));
  return t.exit();
}
