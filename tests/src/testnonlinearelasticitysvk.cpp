// SPDX-FileCopyrightText: 2021-2024 The Ikarus Developers mueller@ibb.uni-stuttgart.de
// SPDX-License-Identifier: LGPL-3.0-or-later

#include <config.h>

#include "testcommon.hh"
#include "testnonlinearelasticity.hh"

#include <ikarus/utils/init.hh>

using Dune::TestSuite;

int main(int argc, char** argv) {
  using namespace Ikarus;
  Ikarus::init(argc, argv);
  TestSuite t("NonLinearElastic + StVenantKirchhoff Test");

  auto matParameter1 = toLamesFirstParameterAndShearModulus({.emodul = 1000, .nu = 0.3});
  auto matParameter2 = toLamesFirstParameterAndShearModulus({.emodul = 1000, .nu = 0.0});

  StVenantKirchhoff matSVK1(matParameter1);
  StVenantKirchhoff matSVK2(matParameter2);
  auto planeStressMat1 = planeStress(matSVK1, 1e-8);
  auto planeStressMat2 = planeStress(matSVK2, 1e-8);
  auto planeStrainMat = planeStrain(matSVK1);

  t.subTest(NonLinearElasticityLoadControlNRandTR<Grids::Alu>(matSVK1));
  t.subTest(NonLinearElasticityLoadControlNRandTR<Grids::Yasp>(matSVK1));
  t.subTest(NonLinearElasticityLoadControlNRandTR<Grids::IgaSurfaceIn2D>(matSVK1));
  t.subTest(GreenLagrangeStrainTest<2>(planeStressMat1));
  t.subTest(GreenLagrangeStrainTest<2>(planeStrainMat));
  t.subTest(GreenLagrangeStrainTest<3>(matSVK1));
  t.subTest(SingleElementTest(planeStressMat2));

  autoDiffTest<2>(t, planeStressMat1, " nu != 0");
  autoDiffTest<2>(t, planeStressMat2, " nu = 0");
  autoDiffTest<3>(t, matSVK1, " nu != 0");
  autoDiffTest<3>(t, matSVK2, " nu = 0");

  {
    auto grid     = createUGGridFromCorners<2>(CornerDistortionFlag::randomlyDistorted);
    auto gridView = grid->leafGridView();
    /// We artificially apply a Neumann load on the complete boundary
    Dune::BitSetVector<1> neumannVertices(gridView.size(2), true);
    BoundaryPatch neumannBoundary(gridView, neumannVertices);
    t.subTest(checkFESByAutoDiff(
        gridView, power<2>(lagrange<1>()),
        Ikarus::AffordanceCollections::elastoStatics));
        skills(Ikarus::nonLinearElastic(planeStrainMat), volumeLoad<2>(vL), neumannBoundaryLoad(&neumannBoundary, nBL)),
  }
  return t.exit();
}
