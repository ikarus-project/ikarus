// SPDX-FileCopyrightText: 2021-2025 The Ikarus Developers ikarus@ibb.uni-stuttgart.de
// SPDX-License-Identifier: LGPL-3.0-or-later

#include <config.h>

#include "testcantileverbeam.hh"

#include <dune/common/test/testsuite.hh>

#include <Eigen/Core>

#include <ikarus/finiteelements/fefactory.hh>
#include <ikarus/finiteelements/mechanics/displacementpressure.hh>
#include <ikarus/finiteelements/mechanics/materials.hh>
#include <ikarus/finiteelements/mixin.hh>
#include <ikarus/solver/eigenvaluesolver/generalizedeigensolverfactory.hh>

using namespace Ikarus;
using Dune::TestSuite;

template <typename FEType>
auto testAutoDiff(const FEType& fe, double tol = 1e-10) {
  Dune::TestSuite t("Check calculateMatrixImpl() and calculateVectorImpl() by Automatic Differentiation");

  auto n = fe.size();
  Eigen::VectorXd d(n);
  d << 0.0, 0.1, 0.2, 0.3, 0.0, 0.1, 0.4, -0.1, 1e-5;
  double lambda                 = 0.0;
  auto req                      = typename FEType::Requirement(d, lambda);
  const std::string feClassName = Dune::className(fe);

  using AutoDiffBasedFE = Ikarus::AutoDiffFE<FEType, true>;
  AutoDiffBasedFE feAutoDiff(fe);

  Eigen::MatrixXd K, KAutoDiff;
  K.setZero(n, n);
  KAutoDiff.setZero(n, n);

  Eigen::VectorXd R, RAutoDiff;
  R.setZero(n);
  RAutoDiff.setZero(n);

  calculateMatrix(fe, req, Ikarus::MatrixAffordance::stiffness, K);
  calculateMatrix(feAutoDiff, req, Ikarus::MatrixAffordance::stiffness, KAutoDiff);

  checkApproxMatrices(t, K, KAutoDiff, testLocation() + "\nIncorrect stiffness matrices." + feClassName, tol);
  checkSymmetricMatrix(t, K, tol, "K");
  checkSymmetricMatrix(t, KAutoDiff, tol, "KAutoDiff");

  return t;
}

template <typename FEType>
auto testEigenValues(const FEType& fe) {
  TestSuite t;

  auto n = fe.size();

  Eigen::VectorXd d;
  d.setZero(n);
  double lambda = 0.0;
  auto req      = typename FEType::Requirement(d, lambda);

  auto testEVs = [&](const auto& d, const auto& evsExpected, const auto& normfexpected, const auto& expectedStress) {
    req.insertGlobalSolution(d);
    Eigen::MatrixXd k;
    k.setZero(n, n);

    Eigen::VectorXd f;
    f.setZero(n);

    calculateMatrix(fe, req, Ikarus::MatrixAffordance::stiffness, k);
    calculateVector(fe, req, Ikarus::VectorAffordance::forces, f);

    auto essaK = makeIdentitySymEigenSolver<EigenValueSolverType::Eigen>(k);
    essaK.compute();
    auto eigenValuesComputed = essaK.eigenvalues();

    checkApproxVectors(t, eigenValuesComputed, evsExpected, " Incorrect eigenvalues of K ", 1e-10);
    checkScalars(t, f.norm(), normfexpected, " Incorrect norm of f", 1e-10);

    Dune::FieldVector<double, 2> firstQP{0.2113248654052, 0.21132486540520};
    auto S = fe.template calculateAt<ResultTypes::PK2Stress>(req, firstQP);

    checkApproxMatrices(t, S.asMat(), expectedStress, "Incorrect Stress", 1e-8);
  };

  auto expectedEigenValues0 = Eigen::Vector<double, 9>{
      -0.00899964037263862, -2.1729667857278618e-14, 2.946341739174961e-14, 1.1436002929121896e-13, 223.72258617281094,
      260.9992542878448,    260.9992542878448,       671.1409395973144,     671.1409395973158};
  auto expectedEigenValuesD = Eigen::Vector<double, 9>{
      -35.283103786777886, -0.01060548937017789, -5.267913592314902e-15, 4.55693342945812e-14, 142.21651420720613,
      198.0477505746059,   382.0279746326626,    715.9234881273067,      1102.7590489808263};

  auto expectedStress0 = Eigen::Matrix<double, 2, 2>::Zero().eval();
  auto expectedStressD = Eigen::Matrix<double, 2, 2>{
      {79.83415515,   43.14730543},
      {43.14730543, -115.61657004}
  };

  testEVs(d, expectedEigenValues0, 8.038873388460925e-14, expectedStress0);

  d << 0.0, 0.1, 0.2, 0.3, 0.0, 0.1, 0.4, -0.1, 1e-5;
  testEVs(d, expectedEigenValuesD, 278.35552561152656, expectedStressD);

  return t;
};

auto testStuff() {
  TestSuite t("Eigenvalue Test");

  // SETUP
  auto matParameter = toLamesFirstParameterAndShearModulus({.emodul = 1000.0, .nu = 0.49});
  auto kappa        = convertLameConstants(matParameter).toBulkModulus();
  auto matDEV       = Materials::makeOgden<1, Ikarus::PrincipalStretchTags::deviatoric>({matParameter.mu}, {2.0});
  auto matVOL       = Materials::makeMaterialLawFromPenaltyFunction(Materials::PVF1());

  using Grid     = Dune::YaspGrid<2>;
  const double L = 1;

  Dune::FieldVector<double, 2> bbox       = {L, L};
  std::array<int, 2> elementsPerDirection = {1, 1};
  auto grid                               = std::make_shared<Grid>(bbox, elementsPerDirection);
  auto gridView                           = grid->leafGridView();

  using namespace Dune::Functions::BasisFactory;
  auto basis = Ikarus::makeBasis(
      gridView, composite(power<2>(lagrange<1>(), FlatInterleaved{}), lagrange<0>(), BlockedLexicographic{}));

  auto sk      = skills(displacementPressure(planeStrain(matDEV), planeStrain(matVOL), kappa));
  using FEType = decltype(makeFE(basis, sk));
  std::vector<FEType> fes;

  for (auto&& ge : elements(gridView)) {
    fes.emplace_back(makeFE(basis, sk));
    fes.back().bind(ge);
  }

  auto& fe = fes.front();

  auto nDOF = basis.flat().size();
  auto n    = fe.size();

  t.check(n == nDOF);

  // TESTS
  t.subTest(testEigenValues(fe));
  t.subTest(testAutoDiff(fe));

  return t;
}

int main(int argc, char** argv) {
  Ikarus::init(argc, argv);
  Dune::TestSuite t;

  t.subTest(testStuff());

  return t.exit();
}
