// SPDX-FileCopyrightText: 2021-2024 The Ikarus Developers mueller@ibb.uni-stuttgart.de
// SPDX-License-Identifier: LGPL-3.0-or-later

#include <config.h>

#include "testhelpers.hh"

#include <dune/common/test/testsuite.hh>

#include <Eigen/Eigenvalues>

#include <ikarus/finiteelements/mechanics/materials.hh>
#include <ikarus/finiteelements/mechanics/materials/muesli/mueslimaterials.hh>

#include <ikarus/finiteelements/physicshelper.hh>
#include <ikarus/utils/functionsanitychecks.hh>
#include <ikarus/utils/init.hh>
#include <ikarus/utils/nonlinearoperator.hh>

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

template <StrainTags strainTag, typename MaterialImpl>
auto testMaterialWithStrain(const MaterialImpl& mat, const double tol = 1e-13) {
  TestSuite t(mat.name() + " InputStrainMeasure: " + toString(strainTag));
  std::cout << "Test: " << t.name() << " started\n";

  Eigen::Matrix3d e;
  e.setRandom();

  double strainDerivativeFactor = transformStrainAccordingToStrain<strainTag>(e);

  auto ev = toVoigtAndMaybeReduce(e, mat, true);
  static_assert(MaterialImpl::isReduced or
                (decltype(ev)::RowsAtCompileTime == 6 and decltype(ev)::ColsAtCompileTime == 1));

  auto energy    = mat.template storedEnergy<strainTag>(e);
  auto energyV   = mat.template storedEnergy<strainTag>(ev);
  auto stressesV = mat.template stresses<strainTag>(e);

  auto stressesVV = mat.template stresses<strainTag>(ev);

  auto moduliV  = mat.template tangentModuli<strainTag>(e);
  auto moduliVV = mat.template tangentModuli<strainTag>(ev);

  t.check(Dune::FloatCmp::le(std::abs(energyV - energy), tol))
      << std::string("Energy obtained from matrix and from Voigt representation do not coincide \n") << energy
      << "\n and \n"
      << energyV << "\n Diff: " << energy - energyV << " with tol: " << tol;
  if constexpr (requires { mat.impl().template stressesImpl<false>(e); }) {
    auto stresses   = mat.template stresses<strainTag, false>(e);
    auto stressesVM = mat.template stresses<strainTag, false>(ev);
    t.check(isApproxSame(toVoigt(stresses, false), enlargeIfReduced<MaterialImpl>(stressesV), tol))
        << std::string("Voigt representation of stresses does not coincide with matrix representation \n")
        << toVoigt(stresses, false) << "\n and \n"
        << enlargeIfReduced<MaterialImpl>(stressesV) << "\n Diff: \n"
        << toVoigt(stresses, false) - enlargeIfReduced<MaterialImpl>(stressesV);
    t.check(isApproxSame(toVoigt(stressesVM, false), enlargeIfReduced<MaterialImpl>(stressesV), tol))
        << std::string(" stresses obtained with strains from voigt, does not coincide with matrix representation \n")
        << stressesVM << "\n and \n"
        << stressesV;
  }
  t.check(isApproxSame(stressesVV, stressesV, tol))
      << std::string(
             "Voigt representation of stresses obtained with strains from "
             "voigt, does not coincide with matrix representation \n")
      << stressesVV << "\n and \n"
      << stressesV << "\n Diff: \n"
      << stressesVV - stressesV;

  if constexpr (requires { mat.impl().template tangentModuliImpl<false>(e); }) {
    auto moduli = mat.template tangentModuli<strainTag, false>(e);
    if constexpr (!MaterialImpl::isReduced) {
      t.check(isApproxSame(toVoigt(moduli), moduliV, tol))
          << std::string("Voigt representation of tangent moduli does not coincide with Tensor<4> representation \n")
          << toVoigt(moduli) << "\n and \n"
          << moduliV;
    }

    t.check(isApproxSame(moduliVV, moduliV, tol)) << std::string(
                                                         "Voigt representation of tangent moduli obtained with Voigt "
                                                         "object does not coincide with tangent moduli "
                                                         "obtained with matrix object  \n")
                                                  << moduliVV << "\n and \n"
                                                  << moduliV;
  }

  auto f  = [&](auto& xv) { return mat.template storedEnergy<strainTag>(xv); };
  auto df = [&](auto& xv) { return (mat.template stresses<strainTag>(xv) * strainDerivativeFactor).eval(); };

  auto ddf = [&](auto& xv) {
    return (mat.template tangentModuli<strainTag>(xv) * strainDerivativeFactor * strainDerivativeFactor).eval();
  };

  auto nonLinOp    = Ikarus::NonLinearOperator(functions(f, df, ddf), parameter(ev));
  auto subNonLinOp = nonLinOp.template subOperator<1, 2>();

  t.check(utils::checkGradient(nonLinOp, {.draw = false, .writeSlopeStatementIfFailed = true}))
      << std::string("checkGradient Failed");
  t.check(utils::checkHessian(nonLinOp, {.draw = false, .writeSlopeStatementIfFailed = true}))
      << std::string("checkHessian Failed");
  t.check(utils::checkJacobian(subNonLinOp, {.draw = false, .writeSlopeStatementIfFailed = true}))
      << std::string("checkJacobian Failed");

  return t;
}

template <typename Material>
auto testMaterial(Material mat) {
  TestSuite t("testMaterial");
  if constexpr (std::is_same_v<Material, LinearElasticity> or
                std::is_same_v<Material,
                               Ikarus::VanishingStress<std::array<Ikarus::Impl::MatrixIndexPair, 1>({{{2, 2}}}),
                                                       Ikarus::LinearElasticity>>) {
    t.subTest(testMaterialWithStrain<StrainTags::linear>(mat));
  } else if constexpr (std::is_same_v<Material,
                                      Ikarus::VanishingStress<std::array<Ikarus::Impl::MatrixIndexPair, 1>({{{2, 2}}}),
                                                              Ikarus::StVenantKirchhoff>>) {
    t.subTest(testMaterialWithStrain<StrainTags::greenLagrangian>(mat));
  } else {
    if constexpr (Material::isReduced) {
      t.subTest(testMaterialWithStrain<StrainTags::greenLagrangian>(mat, 1e-12));
      t.subTest(testMaterialWithStrain<StrainTags::rightCauchyGreenTensor>(mat, 1e-12));
    } else {
      t.subTest(testMaterialWithStrain<StrainTags::greenLagrangian>(mat));
      t.subTest(testMaterialWithStrain<StrainTags::rightCauchyGreenTensor>(mat));
    }
  }
  return t;
}

int main(int argc, char** argv) {
  Ikarus::init(argc, argv);
  TestSuite t;

  double Emod  = 1000;
  double nu    = 0.25; // Blatz Ko assumes nu = 0.25
  auto matPar_ = YoungsModulusAndPoissonsRatio{.emodul = Emod, .nu = nu};
  auto matPar  = toLamesFirstParameterAndShearModulus(matPar_);

  // Eigen::Matrix3d c;
  // c.setRandom();
  // transformStrainAccordingToStrain<StrainTags::rightCauchyGreenTensor>(c);
  //  auto c = Eigen::Matrix3d::Identity().eval();
  Eigen::Matrix3d c{
      { 0.600872, -0.179083, 0},
      {-0.179083,  0.859121, 0},
      {        0,         0, 1}
  };

  auto e = transformStrain<StrainTags::rightCauchyGreenTensor, StrainTags::greenLagrangian>(c);

  // auto muesliMaterial = MuesliElastic(matPar);
  // auto muesliMaterial2 = MuesliElastic(matPar_);

  // auto energy  = muesliMaterial.storedEnergy<StrainTags::linear>(e);
  // auto stress  = muesliMaterial.stresses<StrainTags::linear>(e);
  // auto tangent = muesliMaterial.tangentModuli<StrainTags::linear>(e);

  // std::cout << "Energy\n" << energy << std::endl;
  // std::cout << "Stress\n" << stress << std::endl;
  // std::cout << "Tangent\n" << tangent << std::endl;

  // auto linE = LinearElasticity(matPar);

  // auto energyLin  = linE.storedEnergy<StrainTags::linear>(e);
  // auto stressLin  = linE.stresses<StrainTags::linear>(e);
  // auto tangentLin = linE.tangentModuli<StrainTags::linear>(e);

  // std::cout << "Energy (LIN)\n" << energyLin << std::endl;
  // std::cout << "Stress (LIN)\n" << stressLin << std::endl;
  // std::cout << "Tangent (LIN)\n" << tangentLin << std::endl;

  auto nhm = MuesliFinite<Muesli::NeoHooke>(matPar, false);
  nhm.material().print(std::cout);


  auto energy  = nhm.storedEnergy<StrainTags::rightCauchyGreenTensor>(c);
  auto stress  = nhm.stresses<StrainTags::rightCauchyGreenTensor>(c);
  auto tangent = nhm.tangentModuli<StrainTags::rightCauchyGreenTensor>(c);

  std::cout << "Energy (NHM)\n" << energy << std::endl;
  std::cout << "Stress (NHM)\n" << stress << std::endl;
  std::cout << "Tangent (NHM)\n" << tangent << std::endl;

  auto nh = NeoHooke(matPar);

  auto energyLin  = nh.storedEnergy<StrainTags::rightCauchyGreenTensor>(c);
  auto stressLin  = nh.stresses<StrainTags::rightCauchyGreenTensor>(c);
  auto tangentLin = nh.tangentModuli<StrainTags::rightCauchyGreenTensor>(c);

  std::cout << "Energy (NH)\n" << energyLin << std::endl;
  std::cout << "Stress (NH)\n" << stressLin << std::endl;
  std::cout << "Tangent (NH)\n" << tangentLin << std::endl;

  return t.exit();
}
