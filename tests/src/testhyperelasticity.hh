// SPDX-FileCopyrightText: 2021-2025 The Ikarus Developers ikarus@ibb.uni-stuttgart.de
// SPDX-License-Identifier: LGPL-3.0-or-later

#include <config.h>

#include "materialresultcollection.hh"
#include "testhelpers.hh"

#include <dune/common/test/testsuite.hh>

#include <Eigen/Eigenvalues>

#include <autodiff/forward/dual.hpp>
#include <autodiff/forward/dual/dual.hpp>
#include <autodiff/forward/dual/eigen.hpp>
#include <autodiff/forward/utils/derivative.hpp>

#include <ikarus/finiteelements/mechanics/materials.hh>
#include <ikarus/finiteelements/mechanics/materials/muesli/mueslifinite.hh>
#include <ikarus/finiteelements/physicshelper.hh>
#include <ikarus/utils/differentiablefunction.hh>
#include <ikarus/utils/functionsanitychecks.hh>

using namespace Ikarus;
using namespace Materials;
using Dune::TestSuite;

inline auto testMatrix() {
  return Eigen::Matrix3d{
      { 0.600872, -0.179083, 0},
      {-0.179083,  0.859121, 0},
      {        0,         0, 1}
  };
}

auto testMatPar() {
  double Emod = 1000;
  double nu   = 0.25; // Blatz Ko assumes nu = 0.25
  return YoungsModulusAndPoissonsRatio{.emodul = Emod, .nu = nu};
}

struct Deformations
{
  static constexpr int dim = 3;
  using MatrixType         = Eigen::Matrix<double, dim, dim>;
  Deformations()           = default;

  template <DeformationState DT>
  MatrixType rightCauchyGreen(double lambda) const {
    if constexpr (DT == DeformationState::Undeformed)
      return undeformed(lambda);
    else if constexpr (DT == DeformationState::Uniaxial)
      return uniaxialTensile(lambda);
    else if constexpr (DT == DeformationState::Biaxial)
      return biaxialTensile(lambda);
    else if constexpr (DT == DeformationState::PureShear)
      return pureShear(lambda);
    else
      return random(lambda);
  }

private:
  /// Convert deformation gradient to right Cauchy-Green tensor
  MatrixType toC(const MatrixType& F) const { return (F.transpose() * F).eval(); }

  MatrixType undeformed([[maybe_unused]] double lambda_) const {
    auto F = MatrixType::Identity().eval();
    return toC(F);
  }

  MatrixType uniaxialTensile(double lambda_) const {
    auto F = MatrixType::Zero().eval();
    F.diagonal() << lambda_, 1.0 / sqrt(lambda_), 1.0 / sqrt(lambda_);
    return toC(F);
  }

  MatrixType biaxialTensile(double lambda_) const {
    auto F = MatrixType::Zero().eval();
    F.diagonal() << lambda_, lambda_, 1.0 / (lambda_ * lambda_);
    return toC(F);
  }

  MatrixType pureShear(double lambda_) const {
    auto F = MatrixType::Zero().eval();
    F.diagonal() << lambda_, 1.0, 1.0 / lambda_;
    return toC(F);
  }

  MatrixType random([[maybe_unused]] double lambda_) const { return testMatrix(); }
};

auto testVolumetricFunctions() {
  const auto detC      = testMatrix().determinant();
  const auto J         = sqrt(detC);
  constexpr double tol = 1e-14;

  TestSuite t("Test volumetric functions by AD");

  auto checkFirstAndSecondDerivativeOfEnergy = [&]<typename VF>(const VF& vf) {
    autodiff::dual2nd x = J;
    auto f              = [&](const auto& xloc) { return vf.storedEnergyImpl(xloc); };
    auto derivs         = derivatives(f, autodiff::wrt(x), autodiff::at(x));
    auto uprime         = vf.firstDerivativeImpl(J);
    auto uprimeprime    = vf.secondDerivativeImpl(J);

    checkScalars(t, uprime, derivs[1], testLocation() + Dune::className<VF>(), tol);
    checkScalars(t, uprimeprime, derivs[2], testLocation() + Dune::className<VF>(), tol);
  };

  auto checkFirstDerivativeOfFirstDerivative = [&]<typename VF>(const VF& vf) {
    autodiff::dual x    = J;
    auto f              = [&](const auto& xloc) { return vf.firstDerivativeImpl(xloc); };
    auto uprimeprime_ad = derivative(f, autodiff::wrt(x), autodiff::at(x));
    auto uprimeprime    = vf.secondDerivativeImpl(J);

    checkScalars(t, uprimeprime, uprimeprime_ad, testLocation() + Dune::className<VF>(), tol);
  };

  auto checkVolumetricFunctionsByAD = [&]<typename VF>(const VF& vf) {
    checkFirstAndSecondDerivativeOfEnergy(vf);
    checkFirstDerivativeOfFirstDerivative(vf);
  };

  checkVolumetricFunctionsByAD(VF1{});
  checkVolumetricFunctionsByAD(VF2{});
  checkVolumetricFunctionsByAD(VF3{});
  checkVolumetricFunctionsByAD(VF4{0.5});
  checkVolumetricFunctionsByAD(VF5{});
  checkVolumetricFunctionsByAD(VF6{});
  checkVolumetricFunctionsByAD(VF7{0.5});
  checkVolumetricFunctionsByAD(VF8{});
  checkVolumetricFunctionsByAD(VF9{});
  checkVolumetricFunctionsByAD(VF10{0.4});
  checkVolumetricFunctionsByAD(VF11{});

  return t;
}

auto recoverNeoHookeTest() {
  using StrainTags::rightCauchyGreenTensor;

  TestSuite t("Recover NeoHooke material");

  auto matPar = testMatPar();
  auto mu     = convertLameConstants(matPar).toShearModulus();
  auto Lambda = convertLameConstants(matPar).toLamesFirstParameter();
  auto K      = convertLameConstants(matPar).toBulkModulus();

  std::array<double, 1> mu_og    = {mu};
  std::array<double, 1> alpha_og = {2.0};

  auto nhFromogdenTotal     = makeOgden<1, PrincipalStretchTags::total>(mu_og, alpha_og, Lambda, VF3{});
  auto nhFromogdenDevi      = makeOgden<1, PrincipalStretchTags::deviatoric>(mu_og, alpha_og);
  auto nhFromInvariantBased = makeInvariantBased<1>({mu / 2.0}, {1}, {0});
  auto nh                   = NeoHooke(toLamesFirstParameterAndShearModulus(matPar));

  auto deformation     = Deformations{};
  constexpr double tol = 1e-14;

  auto checkNHRecoveryImpl = [&]<DeformationState def, typename MAT1, typename MAT2>(const MAT1& mat1,
                                                                                     const MAT2& mat2) {
    auto c           = deformation.rightCauchyGreen<def>(1.37);
    auto energy_mat1 = mat1.template storedEnergy<rightCauchyGreenTensor>(c);
    auto stress_mat1 = mat1.template stresses<rightCauchyGreenTensor>(c);
    auto moduli_mat1 = mat1.template tangentModuli<rightCauchyGreenTensor>(c);

    auto energy_mat2 = mat2.template storedEnergy<rightCauchyGreenTensor>(c);
    auto stress_mat2 = mat2.template stresses<rightCauchyGreenTensor>(c);
    auto moduli_mat2 = mat2.template tangentModuli<rightCauchyGreenTensor>(c);

    const std::string matName = mat1.name() + " and " + mat2.name();

    checkScalars(t, energy_mat2, energy_mat1, testLocation() + matName + "Incorrect Energy.", tol);
    checkApproxVectors(t, stress_mat2, stress_mat1, testLocation() + matName + "Incorrect Stresses.", tol);
    checkApproxMatrices(t, moduli_mat2, moduli_mat1, testLocation() + matName + "Incorrect tangentModuli.", tol);
  };

  auto checkNHRecovery = [&]<DeformationState def>() {
    checkNHRecoveryImpl.operator()<def>(nh, nhFromogdenTotal);
    checkNHRecoveryImpl.operator()<def>(nhFromInvariantBased, nhFromogdenDevi);
  };

  // Checking for these deformation states will indirectly ensure correct transformation of quantities from principal
  // coordinate system to Cartesian coordinate system (even for duplicate principal stretches).
  checkNHRecovery.operator()<DeformationState::Undeformed>();
  checkNHRecovery.operator()<DeformationState::Uniaxial>();
  checkNHRecovery.operator()<DeformationState::Biaxial>();
  checkNHRecovery.operator()<DeformationState::PureShear>();
  checkNHRecovery.operator()<DeformationState::Random>();

  return t;
}

template <typename DEV, DeformationState def>
auto materialResults() {
  if constexpr (std::same_as<DEV, BlatzKoT<double>>) {
    return BlatzKoResults<def>();
  } else if constexpr (std::same_as<DEV, OgdenT<double, 3, PrincipalStretchTags::total>>) {
    return OgdenTotalResults<def>();
  } else if constexpr (std::same_as<DEV, OgdenT<double, 3, PrincipalStretchTags::deviatoric>>) {
    return OgdenDeviatoricResults<def>();
  } else if constexpr (std::same_as<DEV, InvariantBasedT<double, 2>>) {
    return MooneyRivlinResults<def>();
  } else if constexpr (std::same_as<DEV, InvariantBasedT<double, 3>> or
                       std::same_as<DEV, Ikarus::Materials::FiniteStrain<muesli::yeohMaterial>>) {
    return YeohResults<def>();
  } else if constexpr (std::same_as<DEV, ArrudaBoyceT<double>> or
                       std::same_as<DEV, Ikarus::Materials::FiniteStrain<muesli::arrudaboyceMaterial>>) {
    return ArrudaBoyceResults<def>();
  } else if constexpr (std::same_as<DEV, GentT<double>>) {
    return GentResults<def>();
  } else
    static_assert(Dune::AlwaysFalse<DEV>::value, "The requested deviatoric function is not implemented.");
}

template <DeformationState def, Concepts::DeviatoricFunction DEV>
auto testMaterialResult(const DEV& dev) {
  Dune::TestSuite t("Test Deviatoric Function Results for the material model: " + dev.name() +
                    " with deformation type as " + toString(def));
  auto [energyEx, firstDerivativesEx, secondDerivativesEx] = materialResults<DEV, def>();
  auto deformation                                         = Deformations{};
  constexpr double lambda                                  = 1.37;
  auto C = Impl::maybeFromVoigt(deformation.rightCauchyGreen<def>(lambda));
  Eigen::SelfAdjointEigenSolver<decltype(C)> eigensolver{};
  eigensolver.compute(C, Eigen::EigenvaluesOnly);
  Eigen::Vector<double, 3> principalStretches = eigensolver.eigenvalues().array().sqrt().eval();
  auto W                                      = dev.storedEnergyImpl(principalStretches);
  auto dWdLambda                              = dev.firstDerivativeImpl(principalStretches);
  auto ddWdLambda                             = dev.secondDerivativeImpl(principalStretches);
  constexpr double tol                        = 1e-14;

  checkScalars(t, W, energyEx, testLocation() + dev.name() + ": Incorrect Energies.", tol);
  checkApproxVectors(t, dWdLambda, firstDerivativesEx, testLocation() + dev.name() + ": Incorrect first derivatives.",
                     tol);
  checkApproxMatrices(t, ddWdLambda, secondDerivativesEx,
                      testLocation() + dev.name() + ": Incorrect second derivatives.", tol);

  return t;
}

auto testMaterialResults() {
  Dune::TestSuite t("Test Deviatoric Function Results ");
  auto matPar     = testMatPar();
  auto mu         = convertLameConstants(matPar).toShearModulus();
  auto K          = convertLameConstants(matPar).toBulkModulus();
  auto Lambda     = convertLameConstants(matPar).toLamesFirstParameter();
  auto lameMatPar = toLamesFirstParameterAndShearModulus(matPar);

  std::array<double, 3> mu_og    = {2.0 * mu / 3.0, mu / 6.0, mu / 6.0};
  std::array<double, 3> alpha_og = {1.23, 0.59, 0.18};

  auto bk         = BlatzKo(mu);
  auto ogdenTotal = Ogden<3, PrincipalStretchTags::total>(mu_og, alpha_og);
  auto ogdenDevi  = Ogden<3, PrincipalStretchTags::deviatoric>(mu_og, alpha_og);
  auto mr         = InvariantBased<2>({1, 0}, {0, 1}, {mu / 2.0, mu / 2.0});
  auto yeoh       = InvariantBased<3>({1, 2, 3}, {0, 0, 0}, {mu / 2.0, mu / 6.0, mu / 3.0});
  auto ab         = ArrudaBoyce({mu, 0.85});
  auto gent       = Gent({mu, 2.5});

  auto checkForDeformationType = [&]<DeformationState def>() {
    t.subTest(testMaterialResult<def>(bk));
    t.subTest(testMaterialResult<def>(ogdenTotal));
    t.subTest(testMaterialResult<def>(ogdenDevi));
    t.subTest(testMaterialResult<def>(mr));
    t.subTest(testMaterialResult<def>(yeoh));
    t.subTest(testMaterialResult<def>(ab));
    t.subTest(testMaterialResult<def>(gent));
  };

  checkForDeformationType.operator()<DeformationState::Undeformed>();
  checkForDeformationType.operator()<DeformationState::Uniaxial>();
  checkForDeformationType.operator()<DeformationState::Biaxial>();
  checkForDeformationType.operator()<DeformationState::PureShear>();
  checkForDeformationType.operator()<DeformationState::Random>();

  // Check malformed Invariant model constructor
  t.checkThrow([]() { InvariantBased<2>({0, 0}, {0, 1}, {500, 500}); });

  return t;
}
