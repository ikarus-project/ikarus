// SPDX-FileCopyrightText: 2021-2025 The Ikarus Developers mueller@ibb.uni-stuttgart.de
// SPDX-License-Identifier: LGPL-3.0-or-later

#include <config.h>

#include "tests/src/testhelpers.hh"
#include "tests/src/testhyperelasticity.hh"

#include <dune/common/test/testsuite.hh>

#include <ikarus/finiteelements/mechanics/materials.hh>
#include <ikarus/finiteelements/mechanics/materials/numericalmaterialinversion.hh>
#include <ikarus/utils/init.hh>
#include <ikarus/utils/lambertw.hh>

using Dune::TestSuite;
using namespace Ikarus;

namespace Testing {

struct StressStates
{
  using enum DeformationState;
  static constexpr int dim = 3;
  using MatrixType         = Eigen::Matrix<double, dim, dim>;

  static constexpr auto allStressStates = std::array{Undeformed, Uniaxial, Biaxial, PureShear, Random};

  static MatrixType PK2Stress(DeformationState ST, double lambda) {
    if (ST == Undeformed)
      return undeformed(lambda);
    else if (ST == Uniaxial)
      return uniaxialTensile(lambda);
    else if (ST == Biaxial)
      return biaxialTensile(lambda);
    else if (ST == PureShear)
      return pureShear(lambda);
    else
      return random(lambda);
  }

private:
  static MatrixType undeformed(double) { return MatrixType::Zero(); }

  static MatrixType uniaxialTensile(double lambda) {
    MatrixType m = MatrixType::Zero();
    m(0, 0)      = lambda;
    return m;
  }

  static MatrixType biaxialTensile(double lambda) {
    MatrixType m = MatrixType::Zero();
    m(0, 0)      = lambda;
    m(1, 1)      = lambda * 0.6;
    return m;
  }

  static MatrixType pureShear(double lambda) {
    MatrixType m = MatrixType::Zero();
    m(1, 0)      = lambda;
    m(0, 1)      = lambda;
    return m;
  }

  static MatrixType random(double) { return testMatrix(); }
};

auto commonErrorMessage(double tol, const std::string& stressStateName) {
  std::ostringstream oss;
  oss << "Stress (" << stressStateName
      << ") from inversion doesn't coincide with input stress (tol = " << std::scientific << std::setprecision(2) << tol
      << ")";
  std::string errorMessage = oss.str();
  return errorMessage;
}

template <typename Material, StrainTags tag, bool useVoigt = true>
void testMaterialInversionMethod(TestSuite& t, const Material& material, const Eigen::Matrix3d& stressState, double tol,
                                 const std::string& errorMessage) {
  auto [D, strain]    = material.template materialInversion<tag, useVoigt>(stressState);
  auto computedStress = material.template stresses<tag>(strain);
  checkApproxVectors(t, computedStress, toVoigt(stressState, false), errorMessage, tol);
}

template <typename Material>
void testNumericalMaterialInversion(TestSuite& t, const Material& material, const Eigen::Matrix3d& stressState,
                                    double tol, const std::string& errorMessage,
                                    const Eigen::Matrix3d& initialValue = Eigen::Matrix3d::Zero()) {
  auto [D, strain]    = numericalMaterialInversion(material, stressState, initialValue);
  auto computedStress = material.template stresses<Material::strainTag>(strain);
  checkApproxVectors(t, computedStress, toVoigt(stressState, false), errorMessage, tol);
}
} // namespace Testing

template <typename T>
auto testLambertW() {
  TestSuite t;

  // These are some numerical tests against https://github.com/JuliaMath/LambertW.jl and
  // https://github.com/DarkoVeberic/LambertW
  auto testNumbers = std::array<T, 6>{-1.0 / std::numbers::e_v<double>, 0.2, 2.0, 100.0, -0.1};
  auto expectedVals =
      std::array<T, 6>{-1.0, 0.16891597349910956, 0.8526055020137254, 3.38563014029005, -0.11183255915896297};

  for (auto i : Dune::range(testNumbers.size())) {
    auto val = util::lambertW0<T>(testNumbers[i]);
    checkScalars(t, expectedVals[i], val, "Wrong result for LambertW", 1e-10);
  }

  return t;
}

auto testMatPar(double nu, double Emod = 1000) {
  auto matPar = YoungsModulusAndPoissonsRatio{.emodul = Emod, .nu = nu};
  auto mu     = convertLameConstants(matPar).toShearModulus();
  auto Lambda = convertLameConstants(matPar).toLamesFirstParameter();
  return LamesFirstParameterAndShearModulus(Lambda, mu);
}

template <typename Material>
TestSuite testMaterial(double nu, const Eigen::Matrix3d& stressState, const std::string& stressStateName) {
  using namespace Ikarus::Materials;
  TestSuite t("Testing " + Material::name() + " with nu = " + std::to_string(nu));

  double tol        = nu > 0.45 ? (std::same_as<Material, NeoHooke> ? 1e-8 : 1e-9) : 1e-12;
  auto errorMessage = Testing::commonErrorMessage(tol, stressStateName);

  auto matPar = testMatPar(nu);
  Material material(matPar);

  // --- Test using the impl method ---
  {
    auto [D, strainMeasure] = material.materialInversionImpl(stressState);
    auto computedStress     = material.template stresses<Material::strainTag>(strainMeasure);
    checkApproxVectors(t, computedStress, toVoigt(stressState, false), errorMessage, tol);
  }
  // --- Test using the greenLagrangian strain tag ---
  Testing::testMaterialInversionMethod<Material, Ikarus::StrainTags::greenLagrangian>(t, material, stressState, tol,
                                                                                      errorMessage);
  // --- Test using the rightCauchyGreenTensor strain tag ---
  Testing::testMaterialInversionMethod<Material, Ikarus::StrainTags::rightCauchyGreenTensor>(t, material, stressState,
                                                                                             tol, errorMessage);
  // --- Test using the numerical approach ---
  { // Obtain Esave and pertube slightly
    auto [_, Esave] = material.template materialInversion<Ikarus::StrainTags::greenLagrangian, false>(stressState);
    Esave.diagonal().array() += 0.1;
    Testing::testNumericalMaterialInversion(t, material, stressState, tol, errorMessage, Esave);
  }
  // --- Test using the greenLagrangian strain tag and non voigt notation---
  Testing::testMaterialInversionMethod<Material, Ikarus::StrainTags::greenLagrangian, false>(t, material, stressState,
                                                                                             tol, errorMessage);

  return t;
}

template <Concepts::Material Material>
auto testHyperelasticMaterials(const Material& material, const Eigen::Matrix3d& stressState,
                               const std::string& stressStateName) {
  TestSuite t("Testing " + Material::name());
  double tol        = 1e-10;
  auto errorMessage = Testing::commonErrorMessage(tol, stressStateName);

  // --- Test using the numerical method directly ---
  Testing::testNumericalMaterialInversion(t, material, stressState, tol, errorMessage);
  // --- Test using the greenLagrangian strain tag ---
  Testing::testMaterialInversionMethod<Material, Ikarus::StrainTags::greenLagrangian>(t, material, stressState, tol,
                                                                                      errorMessage);
  // --- Test using the rightCauchyGreenTensor strain tag ---
  Testing::testMaterialInversionMethod<Material, Ikarus::StrainTags::rightCauchyGreenTensor>(t, material, stressState,
                                                                                             tol, errorMessage);
  // --- Test using the greenLagrangian strain tag and non voigt notation---
  Testing::testMaterialInversionMethod<Material, Ikarus::StrainTags::greenLagrangian, false>(t, material, stressState,
                                                                                             tol, errorMessage);

  return t;
}

auto testMaterialInversions3D() {
  using namespace Materials;
  using namespace Testing;
  TestSuite t("Material Inversions");

  for (auto nu : std::array{0.0, 0.25, 0.35, 0.4999})
    for (auto lambda : std::array{-0.5, 0.5, -100., 100., -1000.})
      for (auto stressState : StressStates::allStressStates) {
        if (std::fabs(lambda) > 100 && stressState == DeformationState::PureShear)
          break;
        t.subTest(
            testMaterial<StVenantKirchhoff>(nu, StressStates::PK2Stress(stressState, lambda), toString(stressState)));
        t.subTest(testMaterial<NeoHooke>(nu, StressStates::PK2Stress(stressState, lambda), toString(stressState)));
      }
  return t;
}

auto testNumericalInversions3d() {
  using namespace Materials;
  using namespace Testing;

  TestSuite t;

  auto matPar = testMatPar(0.3);
  auto K      = convertLameConstants(matPar).toBulkModulus();
  auto mu     = matPar.mu;
  auto lambda = matPar.lambda;

  std::array<double, 1> mu_og     = {mu};
  std::array<double, 1> alpha_og  = {2.0};
  std::array<double, 3> mu_og2    = {2.0 * mu / 3.0, mu / 6.0, mu / 6.0};
  std::array<double, 3> alpha_og2 = {1.23, 0.59, 0.18};

  auto bk        = makeBlatzKo(mu);
  auto ogden     = makeOgden<1, PrincipalStretchTags::total>(mu_og, alpha_og, lambda, VF3{});
  auto ogden2    = makeOgden<3, PrincipalStretchTags::total>(mu_og2, alpha_og2, lambda, VF2{});
  auto ogdenDev  = makeOgden<1, PrincipalStretchTags::deviatoric>(mu_og, alpha_og, K, VF3{});
  auto ogdenDev2 = makeOgden<3, PrincipalStretchTags::deviatoric>(mu_og2, alpha_og2, K, VF2{});

  auto mr   = makeMooneyRivlin({mu / 2.0, mu / 2.0}, K, VF3{});
  auto yeoh = makeYeoh({mu / 2.0, mu / 6.0, mu / 3.0}, K, VF3{});
  auto ab   = makeArrudaBoyce({mu, 0.85}, K, VF3{});
  auto gent = makeGent({mu, 2.5}, K, VF3{});

  auto matTuple = Dune::makeTupleVector(bk, ogden, ogden2, ogdenDev, ogdenDev2, mr, yeoh, ab, gent);

  Dune::Hybrid::forEach(matTuple, [](const auto& mat) {
    for (auto lambda : std::array{-0.5, 0.5})
      for (auto stressState : StressStates::allStressStates)
        testHyperelasticMaterials(mat, StressStates::PK2Stress(stressState, lambda), toString(stressState));
  });

  return t;
}

int main(int argc, char** argv) {
  Ikarus::init(argc, argv);
  TestSuite t;

  t.subTest(testLambertW<double>());
  t.subTest(testLambertW<long double>());
  t.subTest(testMaterialInversions3D());
  t.subTest(testNumericalInversions3d());

  return t.exit();
}
