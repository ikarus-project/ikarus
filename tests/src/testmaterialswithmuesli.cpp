// SPDX-FileCopyrightText: 2021-2024 The Ikarus Developers mueller@ibb.uni-stuttgart.de
// SPDX-License-Identifier: LGPL-3.0-or-later

#include <config.h>

#include "testhelpers.hh"

#include <dune/common/test/testsuite.hh>

#include <Eigen/Eigenvalues>

#include <ikarus/finiteelements/mechanics/materials.hh>
#include <ikarus/finiteelements/mechanics/materials/muesli/mueslimaterials.hh>
#include <ikarus/finiteelements/physicshelper.hh>
#include <ikarus/utils/init.hh>

using namespace Ikarus;
using namespace Ikarus::Materials;
using Dune::TestSuite;

template <StrainTags strainTag>
double transformStrainAccordingToStrain(auto& e) {
  double strainDerivativeFactor = 1;

  if (strainTag == StrainTags::greenLagrangian or strainTag == StrainTags::linear) {
    e = ((e.transpose() + e + 3 * Eigen::Matrix3d::Identity()) / 10).eval();
    e /= e.array().maxCoeff();
    auto C = (2 * e + Eigen::Matrix3d::Identity()).eval();
    Eigen::EigenSolver<Eigen::Matrix3d> esC(C);
    e                      = 0.5 * (C / esC.eigenvalues().real().maxCoeff() - Eigen::Matrix3d::Identity());
    strainDerivativeFactor = 1;
  } else if (strainTag == StrainTags::rightCauchyGreenTensor) {
    e = (e.transpose() + e).eval();
    Eigen::EigenSolver<Eigen::Matrix3d> esC(e);
    e += (-esC.eigenvalues().real().minCoeff() + 1) * Eigen::Matrix3d::Identity();
    esC.compute(e);
    e /= esC.eigenvalues().real().maxCoeff();

    assert(esC.eigenvalues().real().minCoeff() > 0 &&
           " The smallest eigenvalue is negative this is unsuitable for the tests");

    strainDerivativeFactor = 0.5;
  } else if (strainTag == StrainTags::deformationGradient) {
    e = (e + 3 * Eigen::Matrix3d::Identity()).eval(); // create positive definite matrix
    e = e.sqrt();
  }
  return strainDerivativeFactor;
}

template <StrainTags strainTag, typename MuesliMAT, typename IkarusMAT>
requires(std::is_base_of_v<muesli::smallStrainMaterial, typename MuesliMAT::MaterialModel> or
         std::is_base_of_v<muesli::finiteStrainMaterial, typename MuesliMAT::MaterialModel>)
auto testMaterials(const MuesliMAT& muesliMat, const IkarusMAT& ikarusMat) {
  TestSuite t(MuesliMAT::name() + " vs " + IkarusMAT::name() + " InputStrainMeasure: " + toString(strainTag));

  Eigen::Matrix3d c;
  c.setRandom();
  transformStrainAccordingToStrain<strainTag>(c);

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

  auto nhm = MuesliFinite<Muesli::NeoHooke>(matPar, false);
  auto nh  = NeoHooke(matPar);

  t.subTest(testMaterials<StrainTags::rightCauchyGreenTensor>(nhm, nh));
  t.subTest(testMaterials<StrainTags::deformationGradient>(nhm, nh));
  t.subTest(testMaterials<StrainTags::greenLagrangian>(nhm, nh));



  // auto energy  = nhm.storedEnergy<StrainTags::rightCauchyGreenTensor>(c);
  // auto stress  = nhm.stresses<StrainTags::rightCauchyGreenTensor>(c);
  // auto tangent = nhm.tangentModuli<StrainTags::rightCauchyGreenTensor>(c);

  // std::cout << "Energy (NHM)\n" << energy << std::endl;
  // std::cout << "Stress (NHM)\n" << stress << std::endl;
  // std::cout << "Tangent (NHM)\n" << tangent << std::endl;

  // auto energyLin  = nh.storedEnergy<StrainTags::rightCauchyGreenTensor>(c);
  // auto stressLin  = nh.stresses<StrainTags::rightCauchyGreenTensor>(c);
  // auto tangentLin = nh.tangentModuli<StrainTags::rightCauchyGreenTensor>(c);

  // std::cout << "Energy (NH)\n" << energyLin << std::endl;
  // std::cout << "Stress (NH)\n" << stressLin << std::endl;
  // std::cout << "Tangent (NH)\n" << tangentLin << std::endl;

  return t.exit();
}
