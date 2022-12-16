// SPDX-FileCopyrightText: 2022 The Ikarus Developers mueller@ibb.uni-stuttgart.de
// SPDX-License-Identifier: LGPL-3.0-or-later

#include <config.h>

#include <dune/common/parallel/mpihelper.hh>
#include <dune/common/test/testsuite.hh>
using Dune::TestSuite;

#include "testFEElement.hh"

#include <dune/functions/functionspacebases/basistags.hh>
#include <dune/functions/functionspacebases/lagrangebasis.hh>
#include <dune/functions/functionspacebases/powerbasis.hh>

#include "ikarus/finiteElements/mechanics/linearElastic.hh"

template <typename Basis>
using LinearElasticElement = Ikarus::LinearElastic<Basis>;

int main(int argc, char** argv) {
  Dune::MPIHelper::instance(argc, argv);
  TestSuite t("LinearElasticity");

  using namespace Dune::Functions::BasisFactory;
  auto firstOrderLagrangePrePower2Basis  = power<2>(lagrange<1>(), FlatInterleaved());
  auto secondOrderLagrangePrePower2Basis = power<2>(lagrange<2>(), FlatInterleaved());
  auto firstOrderLagrangePrePower3Basis  = power<3>(lagrange<1>(), FlatInterleaved());
  auto secondOrderLagrangePrePower3Basis = power<3>(lagrange<2>(), FlatInterleaved());

  t.subTest(testFEElement<LinearElasticElement, 2>(firstOrderLagrangePrePower2Basis, "LinearElastic",
                                                   checkGradientFunctor, checkHessianFunctor, checkJacobianFunctor));
  t.subTest(testFEElement<LinearElasticElement, 2>(secondOrderLagrangePrePower2Basis, "LinearElastic",
                                                   checkGradientFunctor, checkHessianFunctor, checkJacobianFunctor));
  t.subTest(testFEElement<LinearElasticElement, 3>(firstOrderLagrangePrePower3Basis, "LinearElastic",
                                                   checkGradientFunctor, checkHessianFunctor, checkJacobianFunctor));
  t.subTest(testFEElement<LinearElasticElement, 3>(secondOrderLagrangePrePower3Basis, "LinearElastic",
                                                   checkGradientFunctor, checkHessianFunctor, checkJacobianFunctor));

  return t.exit();
}
