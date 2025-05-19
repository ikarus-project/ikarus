// SPDX-FileCopyrightText: 2021-2025 The Ikarus Developers ikarus@ibb.uni-stuttgart.de
// SPDX-License-Identifier: LGPL-3.0-or-later

#include <config.h>

#include "tests/src/testhelpers.hh"
#include "tests/src/testhyperelasticity.hh"

#include <dune/common/test/testsuite.hh>

#include <autodiff/forward/dual.hpp>

#include <ikarus/finiteelements/mechanics/materials.hh>
#include <ikarus/finiteelements/mechanics/materials/numericalmaterialinversion.hh>
#include <ikarus/utils/init.hh>
#include <ikarus/utils/lambertw.hh>

using Dune::TestSuite;
using namespace Ikarus;

namespace Testing {

static constexpr auto allStressStates =
    std::array{DeformationState::Undeformed, DeformationState::Uniaxial, DeformationState::Biaxial,
               DeformationState::PureShear, DeformationState::Random};

template <int dim>
struct StressStates
{
  static Eigen::Matrix<double, 3, 3> randomMatrix() {
    return Eigen::Matrix<double, 3, 3>{
        {  0.600872, -0.179083, -0.1193886},
        { -0.179083,  0.859121,    -0.0358},
        {-0.1193886,   -0.0358,          1}
    };
  }

  using enum DeformationState;
  using MatrixType = Eigen::Matrix<double, dim, dim>;

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

  static MatrixType random(double)
  requires(dim == 3)
  {
    return randomMatrix();
  }
  static MatrixType random(double)
  requires(dim == 2)
  {
    return randomMatrix().block<2, 2>(0, 0);
  }
};

auto commonErrorMessage(double tol, const std::string& stressStateName) {
  std::ostringstream oss;
  oss << "Stress (" << stressStateName
      << ") from inversion doesn't coincide with input stress (tol = " << std::scientific << std::setprecision(2) << tol
      << ")";
  std::string errorMessage = oss.str();
  return errorMessage;
}

template <typename Material, StrainTags tag, int dim, bool useVoigt = true>
void testMaterialInversionMethod(TestSuite& t, const Material& material,
                                 const Eigen::Matrix<double, dim, dim>& stressState, double tol,
                                 const std::string& errorMessage) {
  // Input for 2d materials only as voigt
  auto [D, strain]    = material.template materialInversion<tag, useVoigt>(toVoigt(stressState, false));
  auto computedStress = material.template stresses<tag>(Impl::maybeToVoigt(strain, true));

  checkApproxVectors(t, computedStress, toVoigt(stressState, false), errorMessage, tol);
}

template <typename Material, int dim>
void testNumericalMaterialInversion(
    TestSuite& t, const Material& material, const Eigen::Matrix<double, dim, dim>& stressState, double tol,
    const std::string& errorMessage,
    const Eigen::Matrix<double, dim, dim>& initialValue = Eigen::Matrix<double, dim, dim>::Zero()) {
  auto [D, strain]    = numericalMaterialInversion(material, stressState, initialValue);
  auto computedStress = material.template stresses<Material::strainTag>(toVoigt(strain));
  checkApproxVectors(t, computedStress, toVoigt(stressState, false), errorMessage, tol);
}

template <auto n, Concepts::EigenVector6 Vec>
auto testIndicesToBeZero(TestSuite& t, const std::array<std::size_t, n>& indicesToBeZero, const Vec& vec,
                         const std::string& name) {
  std::vector<size_t> indicesNotZero;
  for (auto i : indicesToBeZero)
    if (Dune::FloatCmp::eq(vec[i], 1e-10))
      indicesNotZero.push_back(i);

  Eigen::VectorX<size_t> indicesNotZeroEigen =
      Eigen::Map<Eigen::VectorX<size_t>>(indicesNotZero.data(), indicesNotZero.size());

  t.check(indicesNotZero.empty()) << name << " not zero at " << indicesNotZeroEigen.transpose()
                                  << "\n Vector is: " << vec.transpose();
}

Eigen::Matrix3d planeStrainElasticity(double E, double nu) {
  Eigen::Matrix3d C;
  double factor = E / ((1 + nu) * (1 - 2 * nu));
  C << 1 - nu, nu, 0, nu, 1 - nu, 0, 0, 0, (1 - 2 * nu) / 2;
  return factor * C;
}

Eigen::Matrix3d planeStressElasticity(double E, double nu) {
  Eigen::Matrix3d C;
  double factor = E / (1 - nu * nu);
  C << 1, nu, 0, nu, 1, 0, 0, 0, (1 - nu) / 2;
  return factor * C;
}

Eigen::Matrix<double, 6, 6> elasticityTensor3D(double E, double nu) {
  Eigen::Matrix<double, 6, 6> C;
  double factor = E / ((1 + nu) * (1 - 2 * nu));
  // clang-format off
    C << 1 - nu,  nu,      nu,      0,       0,       0,
         nu,      1 - nu,  nu,      0,       0,       0,
         nu,      nu,      1 - nu,  0,       0,       0,
         0,       0,       0,       (1 - 2 * nu) / 2, 0, 0,
         0,       0,       0,       0,       (1 - 2 * nu) / 2, 0,
         0,       0,       0,       0,       0,       (1 - 2 * nu) / 2;
  // clang-format on
  return factor * C;
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

template <typename Material, int dim>
TestSuite testMaterial(const Material& material, double nu, const Eigen::Matrix<double, dim, dim>& stressState,
                       const std::string& stressStateName) {
  using namespace Ikarus::Materials;
  TestSuite t("Testing " + Material::name() + " with nu = " + std::to_string(nu));

  double tol = nu > 0.45 ? (std::same_as<Material, NeoHooke> ? 1e-8 : 1e-9) : 1e-12;
  // For plane stress we also increase the tolerance
  if (Ikarus::traits::isSpecializationNonTypeAndTypes<Ikarus::Materials::VanishingStress, Material>::value)
    tol = std::max(tol, 1e-9);

  const auto errorMessage = Testing::commonErrorMessage(tol, stressStateName);

  // --- Test using the impl method ---
  if constexpr (dim == 3) {
    auto [D, strainMeasure] = material.materialInversionImpl(stressState);
    auto computedStress     = material.template stresses<Material::strainTag>(strainMeasure);
    checkApproxVectors(t, computedStress, toVoigt(stressState, false), errorMessage, tol);
  }
  // --- Test using the greenLagrangian strain tag ---
  Testing::testMaterialInversionMethod<Material, Ikarus::StrainTags::greenLagrangian, dim>(t, material, stressState,
                                                                                           tol, errorMessage);
  // --- Test using the rightCauchyGreenTensor strain tag ---
  Testing::testMaterialInversionMethod<Material, Ikarus::StrainTags::rightCauchyGreenTensor, dim>(
      t, material, stressState, tol, errorMessage);
  // --- Test using the numerical approach ---
  { // Obtain Esave and pertube slightly
    auto [_, Esave] =
        material.template materialInversion<Ikarus::StrainTags::greenLagrangian, false>(toVoigt(stressState, false));
    Esave.diagonal().array() += 0.1;
    Testing::testNumericalMaterialInversion(t, material, stressState, tol, errorMessage, Esave);
  }
  // --- Test using the greenLagrangian strain tag and non voigt notation---
  Testing::testMaterialInversionMethod<Material, Ikarus::StrainTags::greenLagrangian, dim, false>(
      t, material, stressState, tol, errorMessage);

  return t;
}

template <Concepts::Material Material, int dim>
auto testHyperelasticMaterials(const Material& material, const Eigen::Matrix<double, dim, dim>& stressState,
                               const std::string& stressStateName, double tol) {
  TestSuite t("Testing " + Material::name());
  auto errorMessage = Testing::commonErrorMessage(tol, stressStateName);

  // --- Test using the numerical method directly ---
  Testing::testNumericalMaterialInversion(t, material, stressState, tol, errorMessage);
  // --- Test using the greenLagrangian strain tag ---
  Testing::testMaterialInversionMethod<Material, Ikarus::StrainTags::greenLagrangian, dim>(t, material, stressState,
                                                                                           tol, errorMessage);
  // --- Test using the rightCauchyGreenTensor strain tag ---
  Testing::testMaterialInversionMethod<Material, Ikarus::StrainTags::rightCauchyGreenTensor, dim>(
      t, material, stressState, tol, errorMessage);
  // --- Test using the greenLagrangian strain tag and non voigt notation---
  Testing::testMaterialInversionMethod<Material, Ikarus::StrainTags::greenLagrangian, dim, false>(
      t, material, stressState, tol, errorMessage);

  return t;
}

auto testMaterialInversionForSVKAndNH() {
  using namespace Materials;
  using namespace Testing;
  TestSuite t("Material Inversions");

  auto dimRange = Dune::Hybrid::integralRange(std::integral_constant<int, 2>(), std::integral_constant<int, 4>());

  for (auto nu : std::array{0.0, 0.25, 0.35, 0.4999})
    for (auto lambda : std::array{-0.5, 0.5, -100., 100., -1000.})
      for (auto stressState : allStressStates) {
        if (std::fabs(lambda) > 100 && stressState == DeformationState::PureShear)
          break;

        auto matPar = testMatPar(nu);
        auto svk    = StVenantKirchhoff(matPar);
        auto nh     = NeoHooke(matPar);

        /** For the VanishingStrain and VanishingStress case, by default `numericalMaterialInversion` is used.
         *  Since then a starting value is unavailable, the test is only for smaller stress states,
         *  such that convergence is obtained with the undeformed state being the initial value.
         */
        if (std::abs(lambda) < 100) {
          t.subTest(testMaterial(planeStrain(svk), nu, StressStates<2>::PK2Stress(stressState, lambda),
                                 toString(stressState)));
          t.subTest(testMaterial(planeStrain(nh), nu, StressStates<2>::PK2Stress(stressState, lambda),
                                 toString(stressState)));
          t.subTest(testMaterial(planeStress(svk, 1e-9), nu, StressStates<2>::PK2Stress(stressState, lambda),
                                 toString(stressState)));
          t.subTest(testMaterial(planeStress(nh, 1e-9), nu, StressStates<2>::PK2Stress(stressState, lambda),
                                 toString(stressState)));
        }
      }
  return t;
}

auto testHyperelasticMaterialInversion() {
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
      for (auto stressState : allStressStates) {
        testHyperelasticMaterials(mat, StressStates<3>::PK2Stress(stressState, lambda), toString(stressState), 1e-10);
        testHyperelasticMaterials(planeStrain(mat), StressStates<2>::PK2Stress(stressState, lambda),
                                  toString(stressState), 1e-10);
        testHyperelasticMaterials(planeStress(mat, 1e-9), StressStates<2>::PK2Stress(stressState, lambda),
                                  toString(stressState), 1e-9);
      }
  });

  return t;
}

template <Concepts::Material Material>
auto testPlaneStrainPlaneStressConditions() {
  using namespace Materials;
  using namespace Testing;

  auto matPar    = testMatPar(0.3);
  auto mat3d     = Material(matPar);
  auto mat2D_sig = planeStress(mat3d, 1e-9);
  auto mat2D_eps = planeStrain(mat3d);

  const auto GLStrain = StrainTags::greenLagrangian;

  TestSuite t("Testing PlaneStrain PlaneStress Conditions " + mat3d.name());

  // We start from a "random" 2d stress state
  auto stressState                = StressStates<2>::PK2Stress(DeformationState::Random, 0.0);
  Eigen::Vector<double, 3> Svoigt = toVoigt(stressState, false);

  // Plane Stress case
  auto [_, E_sig] = mat2D_sig.template materialInversion<GLStrain, true>(Svoigt);

  // Plane Strain case
  auto [__, E_eps] = mat2D_eps.template materialInversion<GLStrain, true>(Svoigt);

  auto E_sig_full = enlargeIfReduced<decltype(mat2D_sig)>(E_sig).eval();
  auto E_eps_full = enlargeIfReduced<decltype(mat2D_eps)>(E_eps).eval();

  auto sig_sig_final = mat3d.template stresses<GLStrain, true>(E_sig_full);
  auto sig_eps       = mat3d.template stresses<GLStrain, true>(E_eps_full);

  auto [D_eps, eps_eps_final] = mat3d.template materialInversion<GLStrain, true>(sig_eps);

  auto indicesToBeZero = std::array<size_t, 3>{2, 3, 4};
  testIndicesToBeZero(t, indicesToBeZero, sig_sig_final, "sig_sig_final");
  testIndicesToBeZero(t, indicesToBeZero, eps_eps_final, "eps_eps_final");

  return t;
}

auto checkSVK() {
  using namespace Materials;
  using namespace Testing;

  TestSuite t("Testing SVK 2D");

  double nu = 0.3, Emod = 1000;
  auto matPar    = testMatPar(nu, Emod);
  auto mat3d     = StVenantKirchhoff(matPar);
  auto mat2D_sig = planeStress(mat3d, 1e-9);
  auto mat2D_eps = planeStrain(mat3d);
  auto zero      = Eigen::Vector3d::Zero().eval();
  auto zero3d    = Eigen::Vector<double, 6>::Zero().eval();

  auto C_sig = mat2D_sig.tangentModuli<StrainTags::greenLagrangian, true>(zero);
  auto C_eps = mat2D_eps.tangentModuli<StrainTags::greenLagrangian, true>(zero);

  checkApproxMatrices(t, C_sig, planeStressElasticity(Emod, nu), "TangentModuli not same for plane stress");
  checkApproxMatrices(t, C_eps, planeStrainElasticity(Emod, nu), "TangentModuli not same for plane strain");

  auto [D_sig, _]  = mat2D_sig.template materialInversion<StrainTags::greenLagrangian, true>(zero);
  auto [D_eps, __] = mat2D_eps.template materialInversion<StrainTags::greenLagrangian, true>(zero);

  checkApproxMatrices(t, D_sig, planeStressElasticity(Emod, nu).inverse().eval(),
                      "Inverse tangentModuli not same for plane stress");

  checkApproxMatrices(t, D_eps, planeStrainElasticity(Emod, nu).inverse().eval(),
                      "Inverse tangentModuli not same for plane strain");

  auto [D_sig_numeric, ___] = mat2D_sig.template materialInversion<StrainTags::greenLagrangian, true, true>(zero);

  checkApproxMatrices(t, D_sig, D_sig_numeric, "Inverse tangentModuli should be the same for numeric and non-numeric");

  checkApproxMatrices(t, D_sig_numeric, planeStressElasticity(Emod, nu).inverse().eval(),
                      "Inverse tangentModuli not same for plane stress (numeric approach)");

  // 3D test
  auto C = mat3d.tangentModuli<StrainTags::greenLagrangian, true>(zero3d);
  checkApproxMatrices(t, C, elasticityTensor3D(Emod, nu), "TangentModuli not same for 3D");

  auto [D, ____] = mat3d.template materialInversion<StrainTags::greenLagrangian, true>(zero3d);
  checkApproxMatrices(t, D, elasticityTensor3D(Emod, nu).inverse().eval(), "Inverse tangentModuli not same for 3D");

  auto [D_numeric, _____] = mat3d.template materialInversion<StrainTags::greenLagrangian, true>(zero3d);
  checkApproxMatrices(t, D_numeric, elasticityTensor3D(Emod, nu).inverse().eval(),
                      "Inverse tangentModuli not same for 3D (numeric approach)");

  return t;
}

int main(int argc, char** argv) {
  Ikarus::init(argc, argv);
  TestSuite t;

  t.subTest(testLambertW<double>());
  t.subTest(testLambertW<long double>());
  t.subTest(testLambertW<autodiff::dual1st>());
  t.subTest(testLambertW<autodiff::dual2nd>());

  t.subTest(testMaterialInversionForSVKAndNH());
  t.subTest(testHyperelasticMaterialInversion());

  t.subTest(testPlaneStrainPlaneStressConditions<StVenantKirchhoff>());
  t.subTest(testPlaneStrainPlaneStressConditions<NeoHooke>());

  t.subTest(checkSVK());

  return t.exit();
}
