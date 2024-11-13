// SPDX-FileCopyrightText: 2021-2024 The Ikarus Developers mueller@ibb.uni-stuttgart.de
// SPDX-License-Identifier: LGPL-3.0-or-later

#include <config.h>

#include "testhelpers.hh"

#include <muesli/muesli.h>

#include <dune/common/test/testsuite.hh>

#include <Eigen/Eigenvalues>

#include <ikarus/finiteelements/mechanics/materials.hh>
#include <ikarus/finiteelements/mechanics/materials/muesli/mueslimaterials.hh>
#include <ikarus/finiteelements/physicshelper.hh>
#include <ikarus/utils/init.hh>

using namespace Ikarus;
using namespace Ikarus::Materials;
using Dune::TestSuite;

template <StrainTags strainTag, typename MuesliMAT, typename IkarusMAT>
requires(std::is_base_of_v<muesli::smallStrainMaterial, typename MuesliMAT::MaterialModel> or
         std::is_base_of_v<muesli::finiteStrainMaterial, typename MuesliMAT::MaterialModel>)
auto testMaterials(const MuesliMAT& muesliMat, const IkarusMAT& ikarusMat) {
  TestSuite t(MuesliMAT::name() + " vs " + IkarusMAT::name() + " InputStrainMeasure: " + toString(strainTag));

  Eigen::Matrix3d cc{
      { 0.600872, -0.179083, 0},
      {-0.179083,  0.859121, 0},
      {        0,         0, 1}
  };
  auto c = [&]() {
    if constexpr (strainTag == Ikarus::StrainTags::linear)
      return transformStrain<Ikarus::StrainTags::rightCauchyGreenTensor, Ikarus::StrainTags::greenLagrangian>(cc)
          .eval();
    else
      return transformStrain<Ikarus::StrainTags::rightCauchyGreenTensor, strainTag>(cc).eval();
  }();

  auto energy_muesli  = muesliMat.template storedEnergy<strainTag>(c);
  auto stress_muesli  = muesliMat.template stresses<strainTag>(c);
  auto tangent_muesli = muesliMat.template tangentModuli<strainTag>(c);

  auto energy_ikarus  = ikarusMat.template storedEnergy<strainTag>(c);
  auto stress_ikarus  = ikarusMat.template stresses<strainTag>(c);
  auto tangent_ikarus = ikarusMat.template tangentModuli<strainTag>(c);

  t.check(Dune::FloatCmp::eq(energy_muesli, energy_ikarus, 1e-14)) << testLocation();
  t.check(isApproxSame(stress_muesli, stress_ikarus, 1e-14)) << testLocation();
  t.check(isApproxSame(tangent_muesli, tangent_ikarus, 1e-14)) << testLocation();

  return t;
}

int main(int argc, char** argv) {
  Ikarus::init(argc, argv);
  TestSuite t;

  double Emod  = 1000;
  double nu    = 0.25; // Blatz Ko assumes nu = 0.25
  auto matPar_ = YoungsModulusAndPoissonsRatio{.emodul = Emod, .nu = nu};
  auto matPar  = toLamesFirstParameterAndShearModulus(matPar_);

  auto lin  = LinearElasticity(matPar);
  auto linm = MuesliElastic(matPar);

  t.subTest(testMaterials<StrainTags::linear>(linm, lin));

  auto nhm = Muesli::makeNeoHooke(matPar, false);
  auto nh  = NeoHooke(matPar);

  t.subTest(testMaterials<StrainTags::rightCauchyGreenTensor>(nhm, nh));
  t.subTest(testMaterials<StrainTags::deformationGradient>(nhm, nh));
  t.subTest(testMaterials<StrainTags::greenLagrangian>(nhm, nh));

  auto svk  = StVenantKirchhoff(matPar);
  auto svkm = Muesli::makeSVK(matPar);

  t.subTest(testMaterials<StrainTags::rightCauchyGreenTensor>(svkm, svk));
  t.subTest(testMaterials<StrainTags::deformationGradient>(svkm, svk));
  t.subTest(testMaterials<StrainTags::greenLagrangian>(svkm, svk));

  return t.exit();
}
