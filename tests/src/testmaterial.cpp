// SPDX-FileCopyrightText: 2021-2023 The Ikarus Developers mueller@ibb.uni-stuttgart.de
// SPDX-License-Identifier: LGPL-3.0-or-later

#include <config.h>

#include "testhelpers.hh"

#include <dune/common/test/testsuite.hh>

#include <Eigen/Eigenvalues>

#include <ikarus/finiteelements/mechanics/materials.hh>
#include <ikarus/finiteelements/physicshelper.hh>
#include <ikarus/linearalgebra/nonlinearoperator.hh>
#include <ikarus/utils/functionsanitychecks.hh>
#include <ikarus/utils/init.hh>

using namespace Ikarus;
using Dune::TestSuite;

template <StrainTags strainTag, typename MaterialImpl>
auto testMaterialWithStrain(const MaterialImpl& mat, const double tol = 1e-14) {
  TestSuite t(mat.name() + " InputStrainMeasure: " + toString(strainTag));
  std::cout << "Test: " << t.name() << " started\n";
  Eigen::Matrix3d e;
  e.setRandom();
  double strainDerivativeFactor = 1;
  if (strainTag == StrainTags::greenLagrangian or strainTag == StrainTags::linear) {
    e = ((e.transpose() + e + 3 * Eigen::Matrix3d::Identity()) / 10).eval();
    e /= e.array().maxCoeff();
    auto C = (2 * e + Eigen::Matrix3d::Identity()).eval();
    Eigen::EigenSolver<Eigen::Matrix3d> esC(C);
    e                      = 0.5 * (C / esC.eigenvalues().real().maxCoeff() - Eigen::Matrix3d::Identity());
    strainDerivativeFactor = 1;
  } else if (strainTag == StrainTags::rightCauchyGreenTensor) {
    e = ((e.transpose() + e + 3 * Eigen::Matrix3d::Identity()) / 10).eval();  // create positive definite matrix
    Eigen::EigenSolver<Eigen::Matrix3d> esC(e);
    e /= esC.eigenvalues().real().maxCoeff();
    strainDerivativeFactor = 0.5;
  } else if (strainTag == StrainTags::deformationGradient) {
    e = (e + 3 * Eigen::Matrix3d::Identity()).eval();  // create positive definite matrix
    e = e.sqrt();
  }
  auto ev = toVoigtAndMaybeReduce(e, mat, true);
  static_assert(MaterialImpl::isReduced
                or (decltype(ev)::RowsAtCompileTime == 6 and decltype(ev)::ColsAtCompileTime == 1));

  auto energy     = mat.template storedEnergy<strainTag>(e);
  auto energyV    = mat.template storedEnergy<strainTag>(ev);
  auto stressesV  = mat.template stresses<strainTag>(e);
  auto stressesVV = mat.template stresses<strainTag>(ev);

  auto moduliV  = mat.template tangentModuli<strainTag>(e);
  auto moduliVV = mat.template tangentModuli<strainTag>(ev);

  t.check(Dune::FloatCmp::le(std::abs(energyV - energy), tol))
      << "Energy obtained from matrix and from Voigt representation do not coincide \n"
      << energy << "\n and \n"
      << energyV << "\n Diff: " << energy - energyV << " with tol: " << tol;
  if constexpr (requires { mat.impl().template stressesImpl<false>(e); }) {
    auto stresses   = mat.template stresses<strainTag, false>(e);
    auto stressesVM = mat.template stresses<strainTag, false>(ev);
    t.check(isApproxSame(toVoigt(stresses, false), enlargeIfReduced<MaterialImpl>(stressesV), tol))
        << "Voigt representation of stresses does not coincide with matrix representation \n"
        << toVoigt(stresses, false) << "\n and \n"
        << enlargeIfReduced<MaterialImpl>(stressesV) << "\n Diff: \n"
        << toVoigt(stresses, false) - enlargeIfReduced<MaterialImpl>(stressesV);
    t.check(isApproxSame(toVoigt(stressesVM, false), enlargeIfReduced<MaterialImpl>(stressesV), tol))
        << " stresses obtained with strains from voigt, does not coincide with matrix representation \n"
        << stressesVM << "\n and \n"
        << stressesV;
  }
  t.check(isApproxSame(stressesVV, stressesV, tol)) << "Voigt representation of stresses obtained with strains from "
                                                       "voigt, does not coincide with matrix representation \n"
                                                    << stressesVV << "\n and \n"
                                                    << stressesV << "\n Diff: \n"
                                                    << stressesVV - stressesV;

  if constexpr (requires { mat.impl().template tangentModuliImpl<false>(e); }) {
    auto moduli = mat.template tangentModuli<strainTag, false>(e);
    if constexpr (!MaterialImpl::isReduced) {
      t.check(isApproxSame(toVoigt(moduli), moduliV, tol))
          << "Voigt representation of tangent moduli does not coincide with Tensor<4> representation \n"
          << toVoigt(moduli) << "\n and \n"
          << moduliV;
    }

    t.check(isApproxSame(moduliVV, moduliV, tol))
        << "Voigt representation of tangent moduli obtained with Voigt object does not coincide with tangent moduli "
           "obtained with matrix object  \n"
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

  t.check(checkGradient(nonLinOp, {.draw = false, .writeSlopeStatementIfFailed = true})) << "checkGradient Failed";
  t.check(checkHessian(nonLinOp, {.draw = false, .writeSlopeStatementIfFailed = true})) << "checkHessian Failed";
  t.check(checkJacobian(subNonLinOp, {.draw = false, .writeSlopeStatementIfFailed = true})) << "checkJacobian Failed";

  return t;
}

template <typename Material>
auto testMaterial(Material mat) {
  TestSuite t("testMaterial");
  if constexpr (
      std::is_same_v<
          Material,
          LinearElasticity> or std::is_same_v<Material, Ikarus::VanishingStress<std::array<Ikarus::Impl::StressIndexPair, 1>({{{2, 2}}}), Ikarus::LinearElasticity>>) {
    t.subTest(testMaterialWithStrain<StrainTags::linear>(mat));
  } else if constexpr (std::is_same_v<Material,
                                      Ikarus::VanishingStress<std::array<Ikarus::Impl::StressIndexPair, 1>({{{2, 2}}}),
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

  LamesFirstParameterAndShearModulus matPar{.lambda = 1000, .mu = 500};

  auto svk = StVenantKirchhoff(matPar);
  t.subTest(testMaterial(svk));

  auto nh = NeoHooke(matPar);
  t.subTest(testMaterial(nh));

  auto le = LinearElasticity(matPar);
  t.subTest(testMaterial(le));

  auto leRed = makeVanishingStress<Impl::StressIndexPair{2, 2}>(le, 1e-12);
  t.subTest(testMaterial(leRed));

  auto svkRed = makeVanishingStress<Impl::StressIndexPair{2, 2}>(svk, 1e-12);
  t.subTest(testMaterial(svkRed));

  auto nhRed = makeVanishingStress<Impl::StressIndexPair{2, 2}>(nh, 1e-12);
  t.subTest(testMaterial(nhRed));

  auto nhRed2 = makeVanishingStress<Impl::StressIndexPair{1, 1}, Impl::StressIndexPair{2, 2}>(nh, 1e-12);
  t.subTest(testMaterial(nhRed2));

  auto nhRed3
      = makeVanishingStress<Impl::StressIndexPair{2, 1}, Impl::StressIndexPair{2, 0}, Impl::StressIndexPair{2, 2}>(
          nh, 1e-12);
  t.subTest(testMaterial(nhRed3));

  auto nhRed4 = shellMaterial(nh, 1e-12);
  t.subTest(testMaterial(nhRed4));

  auto nhRed5 = beamMaterial(nh, 1e-12);
  t.subTest(testMaterial(nhRed5));

  auto nhRed6 = planeStress(nh, 1e-12);
  t.subTest(testMaterial(nhRed6));

  return t.exit();
}
