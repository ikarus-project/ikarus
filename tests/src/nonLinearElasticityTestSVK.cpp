// SPDX-FileCopyrightText: 2022 The Ikarus Developers mueller@ibb.uni-stuttgart.de
// SPDX-License-Identifier: LGPL-3.0-or-later

//
#include <config.h>

#include "common.hh"
#include "nonLinearElasticityTest.hh"
#include "testHelpers.hh"

using Dune::TestSuite;

int main(int argc, char** argv) {
  Ikarus::init(argc, argv);
  TestSuite t;

  auto matParameter = Ikarus::toLamesFirstParameterAndShearModulus({.emodul = 1000, .nu = 0.3});

  Ikarus::StVenantKirchhoff matSVK(matParameter);

  t.subTest(NonLinearElasticityLoadControlNRandTR<Grids::Alu>(matSVK));
  t.subTest(NonLinearElasticityLoadControlNRandTR<Grids::Yasp>(matSVK));
  t.subTest(NonLinearElasticityLoadControlNRandTR<Grids::Iga>(matSVK));
  return t.exit();
}
