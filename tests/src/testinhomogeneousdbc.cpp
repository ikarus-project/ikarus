// SPDX-FileCopyrightText: 2021-2025 The Ikarus Developers ikarus@ibb.uni-stuttgart.de
// SPDX-License-Identifier: LGPL-3.0-or-later

#include <config.h>

#include "testelasticstrip.hh"

#include <dune/common/test/testsuite.hh>
#include <dune/grid/uggrid.hh>

#include <Eigen/Core>

#include <ikarus/finiteelements/fefactory.hh>
#include <ikarus/finiteelements/mechanics/enhancedassumedstrains.hh>
#include <ikarus/finiteelements/mechanics/linearelastic.hh>
#include <ikarus/finiteelements/mechanics/materials.hh>
#include <ikarus/finiteelements/mixin.hh>

using namespace Ikarus;
using Dune::TestSuite;

template <bool useArcLength>
auto elasticStripTestResults() {
  if constexpr (not useArcLength) {
    std::array<std::tuple<int, double, double>, 4> expectedResultsSVK = {
        std::make_tuple(6, 1.814746879163122, 1.0), std::make_tuple(6, 1.850732157345016, 1.0),
        std::make_tuple(10, 1.817539174273711, 1.0), std::make_tuple(11, 1.8408528524451402, 1.0)};

    std::array<std::tuple<int, double, double>, 4> expectedResultsNH = {
        std::make_tuple(7, 2.207111977584091, 1.0), std::make_tuple(7, 2.1944518710582974, 1.0),
        std::make_tuple(12, 2.2070109913128926, 1.0), std::make_tuple(10, 2.2078938377725192, 1.0)};

    std::array<decltype(expectedResultsSVK), 2> expectedResults = {expectedResultsSVK, expectedResultsNH};
    return expectedResults;
  } else {
    std::array<std::tuple<int, double, double>, 4> expectedResultsSVK = {
        std::make_tuple(4, 1.023987990175817, 0.481292374796350),
        std::make_tuple(4, 0.5603762787930263, 0.2547863399026589),
        std::make_tuple(8, 1.026515156370454, 0.481209660350352),
        std::make_tuple(8, 1.032752731266185, 0.4811653863224171)};

    std::array<std::tuple<int, double, double>, 4> expectedResultsNH = {
        std::make_tuple(4, 1.319709650556954, 0.4670817259366557),
        std::make_tuple(4, 0.7788676577352444, 0.246974652921659),
        std::make_tuple(8, 1.320939735347231, 0.4668834888190049),
        std::make_tuple(7, 1.32134857972011, 0.466884919792626)};

    std::array<decltype(expectedResultsSVK), 2> expectedResults = {expectedResultsSVK, expectedResultsNH};
    return expectedResults;
  }
}

/** Adapted from Macneal and Harder 1985 (https://doi.org/10.1016/0168-874X(85)90003-4) */
auto linearPatchTestWithIDBC(DBCOption dbcOption) {
  TestSuite t("Patch test for a geometrically linear problem involving inhomogeneous Dirichlet BCs with dbcOption = " +
              toString(dbcOption));

  std::cout << "Started: " << t.name() << std::endl;

  constexpr int gridDim     = 2;
  constexpr int basis_order = 1;

  const double E  = 1000;
  const double nu = 0.25;
  double Dhat     = 0.001;

  Dune::GridFactory<Dune::UGGrid<gridDim>> gridFactory;
  gridFactory.insertVertex({0.0, 0.0});   // 0
  gridFactory.insertVertex({0.24, 0.0});  // 1
  gridFactory.insertVertex({0.0, 0.12});  // 2
  gridFactory.insertVertex({0.24, 0.12}); // 3
  gridFactory.insertVertex({0.04, 0.02}); // 4
  gridFactory.insertVertex({0.18, 0.03}); // 5
  gridFactory.insertVertex({0.08, 0.08}); // 6
  gridFactory.insertVertex({0.16, 0.08}); // 7

  gridFactory.insertElement(Dune::GeometryTypes::quadrilateral, {0, 1, 4, 5});
  gridFactory.insertElement(Dune::GeometryTypes::quadrilateral, {5, 1, 7, 3});
  gridFactory.insertElement(Dune::GeometryTypes::quadrilateral, {6, 7, 2, 3});
  gridFactory.insertElement(Dune::GeometryTypes::quadrilateral, {0, 4, 2, 6});
  gridFactory.insertElement(Dune::GeometryTypes::quadrilateral, {4, 5, 6, 7});
  auto grid = gridFactory.createGrid();

  auto gridView = grid->leafGridView();

  using namespace Dune::Functions::BasisFactory;
  auto basis = Ikarus::makeBasis(gridView, power<2>(lagrange<basis_order>(), FlatInterleaved()));

  auto matParameter = toLamesFirstParameterAndShearModulus({.emodul = E, .nu = nu});
  Materials::LinearElasticity linMat(matParameter);
  auto reducedMat = Materials::planeStress(linMat);

  auto sk      = skills(linearElastic(reducedMat));
  using FEType = decltype(makeFE(basis, sk));
  std::vector<FEType> fes;

  for (auto&& ge : elements(gridView)) {
    fes.emplace_back(makeFE(basis, sk));
    fes.back().bind(ge);
  }

  Ikarus::DirichletValues dirichletValues(basis.flat());

  Dune::FieldVector<double, gridDim> bottomLeftPos{0.0, 0.0};
  Dune::FieldVector<double, gridDim> topLeftPos{0.0, 0.12};
  const auto bottomLeftIndices = utils::globalIndexFromGlobalPosition(basis.flat(), bottomLeftPos);
  const auto topLeftIndices    = utils::globalIndexFromGlobalPosition(basis.flat(), topLeftPos);
  dirichletValues.setSingleDOF(bottomLeftIndices[0], true); // ux at (0, 0) = 0
  dirichletValues.setSingleDOF(bottomLeftIndices[1], true); // uy at (0, 0) = 0
  dirichletValues.setSingleDOF(topLeftIndices[0], true);    // ux at (0, 0.12) = 0
  t.check(dirichletValues.fixedDOFsize() == 3, "Number of fixed DOFs is not equal to 3");

  // Inhomogeneous Boundary Conditions
  auto inhomogeneousDisplacement = [Dhat]<typename T>(const auto& globalCoord, const T& lambda) {
    Eigen::Vector<T, 2> localInhomogeneous;
    if (Dune::FloatCmp::eq(globalCoord[0], 0.24)) {
      localInhomogeneous[0] = Dhat * lambda;
      localInhomogeneous[1] = 0;
    } else
      localInhomogeneous.setZero();
    return localInhomogeneous;
  };

  dirichletValues.storeInhomogeneousBoundaryCondition(inhomogeneousDisplacement);
  t.check(dirichletValues.fixedDOFsize() == 5, "Number of fixed DOFs is not equal to 5");

  auto denseFlatAssembler = makeDenseFlatAssembler(fes, dirichletValues);

  auto req     = typename FEType::Requirement(basis);
  auto& d      = req.globalSolution();
  auto& lambda = req.parameter();
  lambda       = 1.0; // linear case
  auto dRed    = d;
  denseFlatAssembler->bind(req, Ikarus::AffordanceCollections::elastoStatics, dbcOption);
  auto linSolver = LinearSolver(SolverTypeTag::d_LDLT);

  if (dbcOption == DBCOption::Reduced)
    dRed = denseFlatAssembler->createReducedVector(dRed);

  const auto& K = denseFlatAssembler->matrix();
  auto R        = denseFlatAssembler->vector();

  const auto F_dirichlet = utils::obtainForcesDueToIDBC(denseFlatAssembler);
  R += F_dirichlet;

  linSolver.compute(K);
  linSolver.solve(dRed, -R);

  if (dbcOption == DBCOption::Reduced)
    d = denseFlatAssembler->createFullVector(dRed);
  else
    d = dRed;

  Eigen::VectorXd inhomogeneousDisp(basis.flat().dimension());
  dirichletValues.evaluateInhomogeneousBoundaryCondition(inhomogeneousDisp, 1.0);
  for (int i = 0; i < basis.flat().dimension(); ++i)
    if (Dune::FloatCmp::ne(inhomogeneousDisp[i], 0.0))
      d[i] = inhomogeneousDisp[i];

  double expectedSigmaXX = 4.1666666666666667;
  Eigen::VectorXd expectedDisplacement;
  expectedDisplacement.setZero(d.size());
  expectedDisplacement << 0.0, 0.0, 0.001, 0.0, 0.0, -0.000125, 0.001, -0.000125, 0.0001666666666666667,
      -0.0000208333333333333, 0.000750, -0.00003125, 0.0003333333333333333, -0.0000833333333333333,
      0.000666666666666667, -0.0000833333333333333;

  constexpr double tol = 1e-10;
  for (const auto i : Dune::range(d.size()))
    if (std::abs(expectedDisplacement[i]) > tol)
      checkScalars(t, d[i], expectedDisplacement[i], " Incorrect displacement at i = " + std::to_string(i), tol);

  // constant stress states for patch test
  for (const auto& fe : fes) {
    const auto sigma = fe.calculateAt<ResultTypes::linearStress>(req, {0.5, 0.5}).asVec();
    checkScalars(t, sigma[0], expectedSigmaXX, " Incorrect sigma_xx", tol);
  }

  return t;
}

int main(int argc, char** argv) {
  Ikarus::init(argc, argv);
  Dune::TestSuite t;

  auto matParameterSVK = toLamesFirstParameterAndShearModulus({.emodul = 100.0, .nu = 0.3});
  Materials::StVenantKirchhoff matSVK(matParameterSVK);

  LamesFirstParameterAndShearModulus matParameterNH = {.lambda = 24.0, .mu = 6.0};
  Materials::NeoHooke matNH(matParameterNH);

  auto reducedMats = Dune::makeTupleVector(planeStrain(matSVK), planeStrain(matNH));

  auto testRange = Dune::Hybrid::integralRange(std::integral_constant<int, 0>(), std::integral_constant<int, 2>());

  auto testFunctor = [&]<bool useArcLength = false>(DBCOption dbcOption) {
    auto expectedResults = elasticStripTestResults<useArcLength>();
    Dune::Hybrid::forEach(testRange, [&](auto i) {
      t.subTest(elasticStripTest<useArcLength>(dbcOption, reducedMats[i], skills(), 1, expectedResults[i][0]));
      t.subTest(elasticStripTest<useArcLength>(dbcOption, reducedMats[i], skills(), 1, expectedResults[i][1], 2));
      if (dbcOption == DBCOption::Reduced) {
        t.checkThrow<Dune::NotImplemented>(
            [&]() {
              elasticStripTest<useArcLength>(dbcOption, reducedMats[i], skills(eas<EAS::DisplacementGradient>(4)), 2,
                                             expectedResults[i][2]);
            },
            "elasticStripTest with EAS should throw a Dune::NotImplemented for DBCOption::Reduced.");
      } else {
        t.subTest(elasticStripTest<useArcLength>(dbcOption, reducedMats[i], skills(eas<EAS::DisplacementGradient>(4)),
                                                 2, expectedResults[i][2]));
        t.subTest(elasticStripTest<useArcLength>(
            dbcOption, reducedMats[i], skills(eas<EAS::DisplacementGradientTransposed>(4)), 2, expectedResults[i][3]));
      }
    });
  };

  testFunctor(DBCOption::Full);
  testFunctor(DBCOption::Reduced);

  testFunctor.operator()<true>(DBCOption::Full);
  testFunctor.operator()<true>(DBCOption::Reduced);

  t.subTest(linearPatchTestWithIDBC(DBCOption::Full));
  t.subTest(linearPatchTestWithIDBC(DBCOption::Reduced));

  return t.exit();
}
