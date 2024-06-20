// SPDX-FileCopyrightText: 2021-2024 The Ikarus Developers mueller@ibb.uni-stuttgart.de
// SPDX-License-Identifier: LGPL-3.0-or-later

#include <config.h>

#include "checkfebyautodiff.hh"
#include "testcommon.hh"
#include "testnonlinearelasticity.hh"

#include <ikarus/utils/init.hh>

using Dune::TestSuite;

int main(int argc, char** argv) {
  using namespace Ikarus;
  using namespace Dune::Functions::BasisFactory;
  Ikarus::init(argc, argv);
  TestSuite t;

  auto matParameter1 = toLamesFirstParameterAndShearModulus({.emodul = 1000, .nu = 0.3});
  auto matParameter2 = toLamesFirstParameterAndShearModulus({.emodul = 1000, .nu = 0.0});

  StVenantKirchhoff matSVK1(matParameter1);
  StVenantKirchhoff matSVK2(matParameter2);
  auto reducedMat = planeStress(matSVK2, 1e-8);

  t.subTest(NonLinearElasticityLoadControlNRandTR<Grids::Alu>(matSVK1));
  t.subTest(NonLinearElasticityLoadControlNRandTR<Grids::Yasp>(matSVK1));
  t.subTest(NonLinearElasticityLoadControlNRandTR<Grids::IgaSurfaceIn2D>(matSVK1));
  t.subTest(GreenLagrangeStrainTest<2>(reducedMat));
  t.subTest(GreenLagrangeStrainTest<3>(matSVK2));
  t.subTest(SingleElementTest(reducedMat));

  auto vL = []<typename VectorType>([[maybe_unused]] const VectorType& globalCoord, auto& lamb) {
    Eigen::Vector<typename VectorType::field_type, VectorType::dimension> fExt;
    fExt.setZero();
    fExt[1] = 2 * lamb;
    return fExt;
  };

  auto nBL = []<typename VectorType>([[maybe_unused]] const VectorType& globalCoord, auto& lamb) {
    Eigen::Vector<typename VectorType::field_type, VectorType::dimension> fExt;
    fExt.setZero();
    fExt[0] = lamb / 40;
    return fExt;
  };

  {
    auto grid     = createUGGridFromCorners<2>(CornerDistortionFlag::randomlyDistorted);
    auto gridView = grid->leafGridView();
    /// We artificially apply a Neumann load on the complete boundary
    Dune::BitSetVector<1> neumannVertices(gridView.size(2), true);
    BoundaryPatch neumannBoundary(gridView, neumannVertices);
    t.subTest(checkFESByAutoDiff(
        gridView, power<2>(lagrange<1>()),
        skills(Ikarus::nonLinearElastic(reducedMat), volumeLoad<2>(vL), neumannBoundaryLoad(&neumannBoundary, nBL)),
        Ikarus::AffordanceCollections::elastoStatics));
  }

  {
    auto grid     = createUGGridFromCorners<3>(CornerDistortionFlag::randomlyDistorted);
    auto gridView = grid->leafGridView();
    /// We artificially apply a Neumann load on the complete boundary
    Dune::BitSetVector<1> neumannVertices(gridView.size(3), true);
    BoundaryPatch neumannBoundary(gridView, neumannVertices);
    t.subTest(checkFESByAutoDiff(
        gridView, power<3>(lagrange<1>()),
        skills(Ikarus::nonLinearElastic(matSVK1), volumeLoad<3>(vL), neumannBoundaryLoad(&neumannBoundary, nBL)),
        Ikarus::AffordanceCollections::elastoStatics));
  }

  return t.exit();
}
