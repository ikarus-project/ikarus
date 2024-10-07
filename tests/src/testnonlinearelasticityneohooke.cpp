// SPDX-FileCopyrightText: 2021-2024 The Ikarus Developers mueller@ibb.uni-stuttgart.de
// SPDX-License-Identifier: LGPL-3.0-or-later

//
#include <config.h>

#include "testnonlinearelasticity.hh"

#include <ikarus/utils/init.hh>

using Dune::TestSuite;

int main(int argc, char** argv) {
  using namespace Ikarus;
  Ikarus::init(argc, argv);
  TestSuite t("NonLinearElastic + NeoHooke Test");

  auto matParameter1 = toLamesFirstParameterAndShearModulus({.emodul = 1000, .nu = 0.3});
  auto matParameter2 = toLamesFirstParameterAndShearModulus({.emodul = 1000, .nu = 0.0});

  NeoHooke matNH1(matParameter1);
  NeoHooke matNH2(matParameter2);

  auto planeStressMat1 = planeStress(matNH1, 1e-8);
  auto planeStressMat2 = planeStress(matNH2, 1e-8);

  t.subTest(NonLinearElasticityLoadControlNRandTR<Grids::Alu>(matNH1));
  t.subTest(NonLinearElasticityLoadControlNRandTR<Grids::Yasp>(matNH1));
  t.subTest(NonLinearElasticityLoadControlNRandTR<Grids::IgaSurfaceIn2D>(matNH1));

  autoDiffTest<2>(t, planeStressMat1, " nu != 0");
  autoDiffTest<2>(t, planeStressMat2, " nu = 0");
  autoDiffTest<3>(t, matNH1, " nu != 0");
  autoDiffTest<3>(t, matNH2, " nu = 0");
  return t.exit();
}
