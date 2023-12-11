// SPDX-FileCopyrightText: 2021-2023 The Ikarus Developers mueller@ibb.uni-stuttgart.de
// SPDX-License-Identifier: LGPL-3.0-or-later

#include <config.h>

#include "testnonlinearelasticity.hh"

using Dune::TestSuite;

int main(int argc, char** argv) {
  Ikarus::init(argc, argv);
  TestSuite t;

  auto matParameter1 = Ikarus::toLamesFirstParameterAndShearModulus({.emodul = 1000, .nu = 0.3});
  auto matParameter2 = Ikarus::toLamesFirstParameterAndShearModulus({.emodul = 1000, .nu = 0.0});

  Ikarus::StVenantKirchhoff matSVK1(matParameter1);
  Ikarus::StVenantKirchhoff matSVK2(matParameter2);
  auto reducedMat = planeStress(matSVK2, 1e-8);

  t.subTest(NonLinearElasticityLoadControlNRandTR<Grids::Alu>(matSVK1));
  t.subTest(NonLinearElasticityLoadControlNRandTR<Grids::Yasp>(matSVK1));
  t.subTest(NonLinearElasticityLoadControlNRandTR<Grids::Iga>(matSVK1));
  t.subTest(GreenLagrangeStrainTest<2>(reducedMat));
  t.subTest(GreenLagrangeStrainTest<3>(matSVK2));
  t.subTest(SingleElementTest(reducedMat));
  t.subTest(checkFEByAutoDiff<2>(reducedMat));
  t.subTest(checkFEByAutoDiff<3>(matSVK1));
  return t.exit();
}
