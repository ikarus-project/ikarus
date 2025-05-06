// SPDX-FileCopyrightText: 2021-2025 The Ikarus Developers ikarus@ibb.uni-stuttgart.de
// SPDX-License-Identifier: LGPL-3.0-or-later

#include <config.h>

#include "testhelpers.hh"
#include "testhyperelasticity.hh"
#include "tests/src/resultcollection.hh"

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

namespace Testing {
auto testMatProp() {
  auto matPar_ = testMatPar();
  auto matPar  = toLamesFirstParameterAndShearModulus(matPar_);
  auto matProp = propertiesFromIkarusMaterialParameters(matPar);
  return std::make_pair(matPar, matProp);
}
} // namespace Testing

template <StrainTags strainTag, typename MuesliMAT, typename IkarusMAT>
auto compareIkarusAndMuesli(const MuesliMAT& muesliMat, const IkarusMAT& ikarusMat) {
  TestSuite t(MuesliMAT::name() + " vs " + IkarusMAT::name() + " InputStrainMeasure: " + toString(strainTag));

  Eigen::Matrix3d C = testMatrix();
  auto strain       = [&]() {
    if constexpr (strainTag == Ikarus::StrainTags::linear) // For linear we use the same as GreenLagrange
      return transformStrain<Ikarus::StrainTags::rightCauchyGreenTensor, Ikarus::StrainTags::greenLagrangian>(C).eval();
    else
      return transformStrain<Ikarus::StrainTags::rightCauchyGreenTensor, strainTag>(C).eval();
  }();

  auto tol = Testing::isPlaneStress<IkarusMAT> ? 1e-8 : 1e-14;

  auto energy_muesli = muesliMat.template storedEnergy<strainTag>(strain);
  auto stress_muesli = muesliMat.template stresses<strainTag>(strain);
  auto moduli_muesli = muesliMat.template tangentModuli<strainTag>(strain);

  auto energy_ikarus = ikarusMat.template storedEnergy<strainTag>(strain);
  auto stress_ikarus = ikarusMat.template stresses<strainTag>(strain);
  auto moduli_ikarus = ikarusMat.template tangentModuli<strainTag>(strain);

  checkScalars(t, energy_muesli, energy_ikarus, "Incorrect energy", tol);
  checkApproxMatrices(t, stress_muesli, stress_ikarus, "Incorrect stresses", tol);
  checkApproxMatrices(t, moduli_muesli, moduli_ikarus, "Incorrect tangentModuli", tol);

  return t;
}

// template <typename MAT>
auto checkConstructors() {
  TestSuite t("Check Constructors");
  auto test = [&]<typename MAT>(muesli::materialProperties matPar) {
    auto mm = MAT(matPar);

    auto mm1{mm};
    auto mm2 = mm1;

    auto mm3 = std::forward<MAT>(mm);

    t.check(mm1.material().check()) << testLocation();
    t.check(mm2.material().check()) << testLocation();
    t.check(mm3.material().check()) << testLocation();

    t.check(mm1.materialPoint().get() != NULL) << testLocation();
    t.check(mm2.materialPoint().get() != NULL) << testLocation();
    t.check(mm3.materialPoint().get() != NULL) << testLocation();
  };

  auto [matPar, matProp] = Testing::testMatProp();
  test.operator()<FiniteStrain<muesli::neohookeanMaterial>>(matProp);
  test.operator()<FiniteStrain<muesli::svkMaterial>>(matProp);
  test.operator()<SmallStrain<muesli::elasticIsotropicMaterial>>(matProp);

  makeMuesliNeoHooke(YoungsModulusAndBulkModulus{1000, 500});
  makeMuesliNeoHooke(YoungsModulusAndPoissonsRatio{1000, 0.2});
  makeMuesliSVK(YoungsModulusAndLamesFirstParameter{1000, 500});

  return t;
}

template <DeformationState def, typename DEV>
auto testDeviatoricMaterialResultMuesli(const DEV& dev) {
  Dune::TestSuite t("Test Deviatoric Function Results for the material model: " + dev.name() +
                    " with deformation type as " + toString(def));

  // This is more or less copied from tests/src/testhyperelasticity.hh
  auto [energyEx, firstDerivativesEx, secondDerivativesEx] = materialResults<DEV, def>();
  auto deformation                                         = Deformations{};
  constexpr double lambda                                  = 1.37;
  auto C = Impl::maybeFromVoigt(deformation.rightCauchyGreen<def>(lambda));

  auto W               = dev.template storedEnergy<StrainTags::rightCauchyGreenTensor>(C);
  auto stress          = dev.template stresses<StrainTags::rightCauchyGreenTensor, false>(C);
  constexpr double tol = 1e-14;

  // Transform firstDerivativesEx
  Eigen::SelfAdjointEigenSolver<decltype(C)> eigensolver{};
  eigensolver.compute(C);
  Eigen::Vector<double, 3> principalStretches = eigensolver.eigenvalues().array().sqrt().eval();
  auto N                                      = eigensolver.eigenvectors().eval();
  Eigen::Vector3d SDev                        = (firstDerivativesEx.array() / principalStretches.array()).eval();
  Eigen::Matrix3d stressEx                    = (N * SDev.asDiagonal() * N.transpose()).eval();

  // Check
  checkScalars(t, W, energyEx, testLocation() + dev.name() + ": Incorrect Energies.", tol);
  checkApproxMatrices(t, stress, stressEx, testLocation() + dev.name() + ": Incorrect stress.", tol);

  return t;
}

template <typename MAT>
auto testDeviatoricPart(TestSuite& t, const MAT& mat) {
  t.subTest(testDeviatoricMaterialResultMuesli<DeformationState::Undeformed>(mat));
  t.subTest(testDeviatoricMaterialResultMuesli<DeformationState::Uniaxial>(mat));
  t.subTest(testDeviatoricMaterialResultMuesli<DeformationState::Biaxial>(mat));
  t.subTest(testDeviatoricMaterialResultMuesli<DeformationState::PureShear>(mat));
  t.subTest(testDeviatoricMaterialResultMuesli<DeformationState::Random>(mat));
}

auto testMuesliAgainstIkarus() {
  TestSuite t("Test Muesli Materials against Ikarus Materials");
  using enum StrainTags;

  auto callMaterialComparisonTest = [&](const auto& mueslimat, const auto& ikarusmat) {
    t.subTest(compareIkarusAndMuesli<rightCauchyGreenTensor>(mueslimat, ikarusmat));
    t.subTest(compareIkarusAndMuesli<deformationGradient>(mueslimat, ikarusmat));
    t.subTest(compareIkarusAndMuesli<greenLagrangian>(mueslimat, ikarusmat));
  };

  auto [matPar, matProp] = Testing::testMatProp();
  auto mu                = matPar.mu;
  auto K                 = convertLameConstants(matPar).toBulkModulus();

  // Linear Elasticity
  auto lin  = LinearElasticity(matPar);
  auto linm = SmallStrain(matPar);

  t.subTest(compareIkarusAndMuesli<linear>(linm, lin));

  // NeoHooke
  auto nhm = makeMuesliNeoHooke(matPar, false);
  auto nh  = NeoHooke(matPar);
  callMaterialComparisonTest(nhm, nh);

  // SVK
  auto svkm = makeMuesliSVK(matPar);
  auto svk  = StVenantKirchhoff(matPar);
  callMaterialComparisonTest(svkm, svk);

  // ArrudaBoyce
  auto abm = makeMuesliArrudaBoyce(matPar.mu, 0.85, .0, true);
  auto ab  = makeArrudaBoyce({matPar.mu, 0.85});
  callMaterialComparisonTest(abm, ab);
  testDeviatoricPart(t, abm);

  // Yeoh
  auto yaohMatPar = std::array{mu / 2.0, mu / 6.0, mu / 3.0};
  auto yeohm      = makeMuesliYeoh(yaohMatPar, 0.0);
  auto yeoh       = makeYeoh(yaohMatPar);
  callMaterialComparisonTest(yeohm, yeoh);
  testDeviatoricPart(t, yeohm);

  // NeoHooke Deviatoric
  auto nhDevi  = makeOgden<1, PrincipalStretchTags::deviatoric>({mu}, {2.0}, K, VF3{});
  auto nhDevim = makeMuesliNeoHooke(matPar, true);
  callMaterialComparisonTest(nhDevim, nhDevi);

  // Test 2D Case
  auto nhplaneStrain  = planeStrain(nh);
  auto nhplaneStrainm = planeStrain(nhm);
  callMaterialComparisonTest(nhplaneStrainm, nhplaneStrain);

  auto nhplaneStress  = planeStress(nh);
  auto nhplaneStressm = planeStress(nhm);
  callMaterialComparisonTest(nhplaneStress, nhplaneStressm);

  // Test name Function
  t.check(svkm.name() == "FiniteStrain: SvkMaterial");
  t.check(nhm.name() == "FiniteStrain: NeohookeanMaterial");
  t.check(linm.name() == "SmallStrain: ElasticIsotropicMaterial");

  return t;
}

int main(int argc, char** argv) {
  Ikarus::init(argc, argv);
  TestSuite t("Muesli Materials Test");

  t.subTest(checkConstructors());
  t.subTest(testMuesliAgainstIkarus());

  return t.exit();
}
