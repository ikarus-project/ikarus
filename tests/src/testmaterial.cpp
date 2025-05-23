// SPDX-FileCopyrightText: 2021-2025 The Ikarus Developers ikarus@ibb.uni-stuttgart.de
// SPDX-License-Identifier: LGPL-3.0-or-later

#include <config.h>

#include "testhelpers.hh"

#include <dune/common/test/testsuite.hh>

#include <Eigen/Eigenvalues>

#include <ikarus/finiteelements/mechanics/materials.hh>
#if ENABLE_MUESLI
  #include <ikarus/finiteelements/mechanics/materials/muesli/mueslimaterials.hh>
#endif
#include <ikarus/finiteelements/physicshelper.hh>
#include <ikarus/utils/differentiablefunction.hh>
#include <ikarus/utils/functionsanitychecks.hh>
#include <ikarus/utils/init.hh>

using namespace Ikarus;
using namespace Ikarus::Materials;
using Dune::TestSuite;

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

  auto fl  = [&](auto& xv) { return mat.template storedEnergy<strainTag>(xv); };
  auto dfl = [&](auto& xv) { return (mat.template stresses<strainTag>(xv) * strainDerivativeFactor).eval(); };

  auto ddfl = [&](auto& xv) {
    return (mat.template tangentModuli<strainTag>(xv) * strainDerivativeFactor * strainDerivativeFactor).eval();
  };

  auto f  = Ikarus::makeDifferentiableFunction(functions(fl, dfl, ddfl), ev);
  auto df = derivative(f);

  t.check(utils::checkGradient(f, ev, {.draw = false, .writeSlopeStatementIfFailed = true}))
      << std::string("checkGradient Failed");
  t.check(utils::checkHessian(f, ev, {.draw = false, .writeSlopeStatementIfFailed = true}))
      << std::string("checkHessian Failed");
  t.check(utils::checkJacobian(df, ev, {.draw = false, .writeSlopeStatementIfFailed = true}))
      << std::string("checkJacobian Failed");

  return t;
}

template <typename Material>
auto testMaterial(Material mat) {
  TestSuite t("testMaterial");
  if constexpr (Material::isLinear) {
    t.subTest(testMaterialWithStrain<StrainTags::linear>(mat));
  } else if constexpr (std::is_same_v<Material, VanishingStress<std::array<MatrixIndexPair, 1>({{{2, 2}}}),
                                                                Ikarus::Materials::StVenantKirchhoff>>) {
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

template <StrainTags strainTag, typename MaterialImpl>
auto testPlaneStrainAgainstPlaneStress(const double tol = 1e-10) {
  TestSuite t(MaterialImpl::name() + " InputStrainMeasure: " + toString(strainTag));
  std::cout << "TestPlaneStrinAgainstPlaneStress: " << t.name() << " started\n";

  Eigen::Matrix3d e;
  e.setRandom();

  transformStrainAccordingToStrain<strainTag>(e);

  // instantiate material models
  LamesFirstParameterAndShearModulus matPar{.lambda = 0, .mu = 1000}; // nu = 0

  auto mat            = MaterialImpl{matPar};
  auto planeStrainMat = planeStrain(mat);
  auto planeStressMat = planeStress(mat);

  // energy should be the same for plane stress and plane strain for nu = 0
  auto energies = std::array<double, 2>{planeStrainMat.template storedEnergy<strainTag>(e),
                                        planeStressMat.template storedEnergy<strainTag>(e)};

  t.check(Dune::FloatCmp::le(std::abs(energies[0] - energies[1]), tol))
      << "Energies for plane strain and plane stress should be the same but are"
      << "\n"
      << energies[0] << " and " << energies[1] << "\n Diff: " << energies[0] - energies[1] << " with tol: " << tol;

  // Stresses should be the same for nu = 0
  auto stressPlaneStrain = planeStrainMat.template stresses<strainTag>(e);
  auto stressPlaneStress = planeStressMat.template stresses<strainTag>(e);

  t.check(isApproxSame(stressPlaneStrain, stressPlaneStress, tol))
      << "Stresses for plane strain and plane stress should be the same but are"
      << "\n"
      << stressPlaneStrain << "\nand\n " << stressPlaneStress << "\nDiff:\n"
      << stressPlaneStrain - stressPlaneStress << " with tol: " << tol;

  // If we compare the plain stress material tensor with plain strain material tensor it should be the same for nu = 0
  auto matTangentPlaneStrain = planeStrainMat.template tangentModuli<strainTag>(e);
  auto matTangentPlaneStress = planeStressMat.template tangentModuli<strainTag>(e);

  t.check(isApproxSame(matTangentPlaneStrain, matTangentPlaneStress, tol))
      << "Material Tangent for plane strain and plane stress should be the same but are"
      << "\n"
      << matTangentPlaneStrain << "\nand\n"
      << matTangentPlaneStress << "\n Diff:\n"
      << matTangentPlaneStrain - matTangentPlaneStress << " with tol: " << tol;

  // Sanity check, for nu != 0 the quantities should differ
  matPar         = {.lambda = 1000, .mu = 50}; // nu = 0
  mat            = MaterialImpl{matPar};
  planeStrainMat = planeStrain(mat);
  planeStressMat = planeStress(mat);

  energies = std::array<double, 2>{planeStrainMat.template storedEnergy<strainTag>(e),
                                   planeStressMat.template storedEnergy<strainTag>(e)};

  t.check(not Dune::FloatCmp::eq(energies[0], energies[1], tol))
      << "Energies for plane strain and plane stress should not be the same but are"
      << "\n"
      << energies[0] << " and " << energies[1] << " with tol: " << tol;

  stressPlaneStrain = planeStrainMat.template stresses<strainTag>(e);
  stressPlaneStress = planeStressMat.template stresses<strainTag>(e);

  t.check(not isApproxSame(stressPlaneStrain, stressPlaneStress, tol))
      << "Stresses for plane strain and plane stress should not be the same but are"
      << "\n"
      << stressPlaneStrain << "\nand\n " << stressPlaneStress << " with tol: " << tol;

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

  auto leRed = makeVanishingStress<MatrixIndexPair{2, 2}>(le, 1e-12);
  t.subTest(testMaterial(leRed));

  auto svkRed = makeVanishingStress<MatrixIndexPair{2, 2}>(svk, 1e-12);
  t.subTest(testMaterial(svkRed));

  auto nhRed = makeVanishingStress<MatrixIndexPair{2, 2}>(nh, 1e-12);
  t.subTest(testMaterial(nhRed));

  auto nhRed2 = makeVanishingStress<MatrixIndexPair{1, 1}, MatrixIndexPair{2, 2}>(nh, 1e-12);
  t.subTest(testMaterial(nhRed2));

  auto nhRed3 = makeVanishingStress<MatrixIndexPair{2, 1}, MatrixIndexPair{2, 0}, MatrixIndexPair{2, 2}>(nh, 1e-12);
  t.subTest(testMaterial(nhRed3));

  auto nhRed4 = shellMaterial(nh, 1e-12);
  t.subTest(testMaterial(nhRed4));

  auto nhRed5 = beamMaterial(nh, 1e-12);
  t.subTest(testMaterial(nhRed5));

  auto nhRed6 = planeStress(nh, 1e-12);
  t.subTest(testMaterial(nhRed6));

  auto svkPlaneStrain = planeStrain(svk);
  t.subTest(testMaterial(svkPlaneStrain));

  auto nhPlaneStrain = planeStrain(nh);
  t.subTest(testMaterial(nhPlaneStrain));

  auto linPlaneStrain = planeStrain(le);
  t.subTest(testMaterialWithStrain<StrainTags::linear>(linPlaneStrain));

  auto nhRed8 = makeVanishingStrain<MatrixIndexPair{1, 1}, MatrixIndexPair{2, 2}>(nh);
  t.subTest(testMaterial(nhRed8));

  // Hyperelasticity
  double mu     = matPar.mu;
  double lambda = matPar.lambda;
  double K      = convertLameConstants(matPar).toBulkModulus();

  auto bk = makeBlatzKo(mu);
  t.subTest(testMaterial(bk));

  // Material parameters (example values)
  std::array<double, 1> mu_og    = {mu};
  std::array<double, 1> alpha_og = {2.0};

  std::array<double, 3> mu_og2    = {2.0 * mu / 3.0, mu / 6.0, mu / 6.0};
  std::array<double, 3> alpha_og2 = {1.23, 0.59, 0.18};

  auto ogden     = makeOgden<1, PrincipalStretchTags::total>(mu_og, alpha_og, lambda, VF3{});
  auto ogden2    = makeOgden<3, PrincipalStretchTags::total>(mu_og2, alpha_og2, lambda, VF2{});
  auto ogdenDev  = makeOgden<1, PrincipalStretchTags::deviatoric>(mu_og, alpha_og, K, VF3{});
  auto ogdenDev2 = makeOgden<3, PrincipalStretchTags::deviatoric>(mu_og2, alpha_og2, K, VF2{});

  t.subTest(testMaterial(ogden));
  t.subTest(testMaterial(ogden2));
  t.subTest(testMaterial(ogdenDev));
  t.subTest(testMaterial(ogdenDev2));

  auto mr   = makeMooneyRivlin({mu / 2.0, mu / 2.0});
  auto yeoh = makeYeoh({mu / 2.0, mu / 6.0, mu / 3.0});
  auto ab   = makeArrudaBoyce({mu, 0.85});
  auto gent = makeGent({mu, 2.5});

  t.subTest(testMaterial(mr));
  t.subTest(testMaterial(yeoh));
  t.subTest(testMaterial(ab));
  t.subTest(testMaterial(gent));

  auto psOgden  = planeStrain(ogden);
  auto pstOgden = planeStress(ogden2, 1e-12);
  auto psmr     = planeStrain(mr);

  t.subTest(testMaterial(psOgden));
  t.subTest(testMaterial(pstOgden));
  t.subTest(testMaterial(psmr));

  t.subTest(testPlaneStrainAgainstPlaneStress<StrainTags::linear, LinearElasticity>());
  t.subTest(testPlaneStrainAgainstPlaneStress<StrainTags::greenLagrangian, StVenantKirchhoff>());
  t.subTest(testPlaneStrainAgainstPlaneStress<StrainTags::rightCauchyGreenTensor, NeoHooke>());

#if ENABLE_MUESLI
  // Muesli
  auto muesliLin = makeMuesliLinearElasticity(matPar);
  t.subTest(testMaterial(muesliLin));

  auto muesliSVK = makeMuesliSVK(matPar);
  t.subTest(testMaterial(muesliSVK));

  auto muesliNeoHooke = makeMuesliNeoHooke(matPar, false);
  t.subTest(testMaterial(muesliNeoHooke));

  auto muesliNeoHookeReg = makeMuesliNeoHooke(matPar, true);
  t.subTest(testMaterial(muesliNeoHookeReg));

  auto muesliMR = makeMooneyRivlin({K, mu / 2, mu / 2});
  t.subTest(testMaterial(muesliMR));

  auto muesliYeoh = makeMuesliYeoh({mu / 2, mu / 6, mu / 3}, K);
  t.subTest(testMaterial(muesliYeoh));

  auto muesliAB = makeMuesliArrudaBoyce(mu, 0.85, K, false);
  t.subTest(testMaterial(muesliAB));

  auto muesliSVKPlaneStrain = planeStrain(muesliSVK);
  t.subTest(testMaterial(muesliSVKPlaneStrain));

  auto muesliNHPlaneStress = planeStress(muesliNeoHooke, 1e-11);
  t.subTest(testMaterial(muesliNHPlaneStress));
#endif
  return t.exit();
}
