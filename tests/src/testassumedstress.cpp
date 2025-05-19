// SPDX-FileCopyrightText: 2021-2025 The Ikarus Developers ikarus@ibb.uni-stuttgart.de
// SPDX-License-Identifier: LGPL-3.0-or-later

#include <config.h>

#include "testassumedstress.hh"

#include "testfeelement.hh"

#include <dune/common/test/testsuite.hh>
#include <dune/functions/functionspacebases/basistags.hh>
#include <dune/functions/functionspacebases/lagrangebasis.hh>
#include <dune/functions/functionspacebases/powerbasis.hh>

#include <ikarus/finiteelements/mechanics/assumedstress.hh>
#include <ikarus/finiteelements/mechanics/linearelastic.hh>
#include <ikarus/utils/init.hh>

using Dune::TestSuite;

int main(int argc, char** argv) {
  Ikarus::init(argc, argv);
  TestSuite t("Assumed Stress");

  using namespace Dune::Functions::BasisFactory;
  using namespace Ikarus::Materials;

  auto firstOrderLagrangePrePower2Basis = power<2>(lagrange<1>(), FlatInterleaved());
  auto firstOrderLagrangePrePower3Basis = power<3>(lagrange<1>(), FlatInterleaved());
  constexpr auto randomlyDistorted      = CornerDistortionFlag::randomlyDistorted;
  constexpr auto unDistorted            = CornerDistortionFlag::unDistorted;

  auto linearElasticFunc3D = [](const Ikarus::YoungsModulusAndPoissonsRatio& parameter) {
    LinearElasticity lin(Ikarus::toLamesFirstParameterAndShearModulus(parameter));
    return Ikarus::linearElastic(lin);
  };
  auto linearElasticFuncPlaneStress = [](const Ikarus::YoungsModulusAndPoissonsRatio& parameter) {
    LinearElasticity lin(Ikarus::toLamesFirstParameterAndShearModulus(parameter));
    auto linPS = planeStress(lin, 1e-12);
    return Ikarus::linearElastic(linPS);
  };
  auto linearElasticFuncPlaneStrain = [](const Ikarus::YoungsModulusAndPoissonsRatio& parameter) {
    LinearElasticity lin(Ikarus::toLamesFirstParameterAndShearModulus(parameter));
    auto linPS = planeStrain(lin);
    return Ikarus::linearElastic(linPS);
  };

  t.subTest(testFEElement(firstOrderLagrangePrePower2Basis, "PS", randomlyDistorted,
                          Dune::ReferenceElements<double, 2>::cube(), linearElasticFuncPlaneStress,
                          Ikarus::skills(Ikarus::assumedStress(5)), Ikarus::AffordanceCollections::elastoStatics,
                          checkJacobianFunctor));

  t.subTest(testFEElement(firstOrderLagrangePrePower2Basis, "PS", randomlyDistorted,
                          Dune::ReferenceElements<double, 2>::cube(), linearElasticFuncPlaneStrain,
                          Ikarus::skills(Ikarus::assumedStress(5)), Ikarus::AffordanceCollections::elastoStatics,
                          checkJacobianFunctor));

  t.subTest(testFEElement(firstOrderLagrangePrePower3Basis, "PS", randomlyDistorted,
                          Dune::ReferenceElements<double, 3>::cube(), linearElasticFunc3D,
                          Ikarus::skills(Ikarus::assumedStress(18)), Ikarus::AffordanceCollections::elastoStatics,
                          checkJacobianFunctor));

  return t.exit();
}
