// SPDX-FileCopyrightText: 2021-2024 The Ikarus Developers mueller@ibb.uni-stuttgart.de
// SPDX-License-Identifier: LGPL-3.0-or-later

//
#include <config.h>

#include "testnonlinearelasticity.hh"

#include <ikarus/utils/init.hh>

using Dune::TestSuite;

int main(int argc, char** argv) {
  Ikarus::init(argc, argv);
  TestSuite t;

  auto matParameter = Ikarus::toLamesFirstParameterAndShearModulus({.emodul = 1000, .nu = 0.3});

  Ikarus::NeoHooke matNH(matParameter);

  t.subTest(NonLinearElasticityLoadControlNRandTR<Grids::Alu>(matNH));
  t.subTest(NonLinearElasticityLoadControlNRandTR<Grids::Yasp>(matNH));
  t.subTest(NonLinearElasticityLoadControlNRandTR<Grids::IgaSurfaceIn2D>(matNH));
  return t.exit();
}
