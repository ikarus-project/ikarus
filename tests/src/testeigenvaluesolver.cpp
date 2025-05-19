// SPDX-FileCopyrightText: 2021-2025 The Ikarus Developers ikarus@ibb.uni-stuttgart.de
// SPDX-License-Identifier: LGPL-3.0-or-later
#include <config.h>

#include "dummyproblem.hh"
#include "testhelpers.hh"

#include <dune/common/float_cmp.hh>
#include <dune/common/test/testsuite.hh>
#include <dune/functions/functionspacebases/lagrangebasis.hh>
#include <dune/functions/functionspacebases/powerbasis.hh>
#include <dune/grid/uggrid.hh>
#include <dune/grid/yaspgrid.hh>

#include <ikarus/assembler/simpleassemblers.hh>
#include <ikarus/finiteelements/fefactory.hh>
#include <ikarus/finiteelements/mechanics/linearelastic.hh>
#include <ikarus/finiteelements/mechanics/loads/volume.hh>
#include <ikarus/finiteelements/mixin.hh>
#include <ikarus/solver/eigenvaluesolver/generalizedeigensolverfactory.hh>
#include <ikarus/utils/basis.hh>
#include <ikarus/utils/dirichletvalues.hh>
#include <ikarus/utils/init.hh>

using Dune::TestSuite;

namespace Testing {

auto consistentRandomNumber(int i) {
  std::mt19937 gen(i);
  std::uniform_real_distribution<double> dist(1.0, 100.0); // Random numbers in [1, 100)
  return dist(gen);
}

auto generateRandomMatrices(int n) {
  Eigen::SparseMatrix<double> M(n, n);
  Eigen::MatrixX<double> MD(n, n);
  MD.setZero();

  for (auto i : Dune::range(n)) {
    double randomNumber = consistentRandomNumber(i);
    M.coeffRef(i, i)    = randomNumber;
    MD(i, i)            = randomNumber;
  }
  return std::make_pair(M, MD);
}
} // namespace Testing

template <Ikarus::Concepts::EigenValueSolver SOL1, Ikarus::Concepts::EigenValueSolver SOL2>
auto testEigenValues(const SOL1& solver1, const SOL2& solver2, std::optional<Eigen::Index> nev_ = std::nullopt,
                     const std::string& message = "") {
  TestSuite t("Test Eigenvalues" + Dune::className<SOL1>() + ", " + Dune::className<SOL2>());
  auto eigenvalues1 = solver1.eigenvalues();
  auto eigenvalues2 = solver2.eigenvalues();

  if (not nev_.has_value())
    checkApproxVectors(t, eigenvalues1, eigenvalues2, message);
  else
    checkApproxVectors(t, eigenvalues1.head(*nev_).eval(), eigenvalues2.head(*nev_).eval(), message);

  return t;
}

// as1 corresponds to solver1, which is tested against solver2
template <Ikarus::Concepts::EigenValueSolver SOL1, Ikarus::Concepts::EigenValueSolver SOL2>
auto testEigenVectors(const SOL1& solver1, const SOL2& solver2, const auto& K, const std::string& message = "") {
  TestSuite t("Test Eigenvalues " + Dune::className<SOL1>() + ", " + Dune::className<SOL2>());

  Eigen::MatrixXd evecs = solver1.eigenvectors();
  Eigen::VectorXd evals = solver2.eigenvalues();

  Eigen::VectorXd T = (evecs.transpose() * K * evecs).diagonal();

  checkApproxVectors(t, T, evals, message, 1e-10);

  return t;
}

auto testEigenSolvers() {
  TestSuite t("EigenValueSolver Test", Dune::TestSuite::ThrowPolicy::ThrowOnRequired);
  using Grid     = Dune::YaspGrid<2>;
  using GridView = Grid::LeafGridView;

  DummyProblem<Grid, true> testCase({5, 5});
  auto& req             = testCase.requirement();
  auto& fes             = testCase.finiteElements();
  auto& dirichletValues = testCase.dirichletValues();

  auto assK  = Ikarus::makeSparseFlatAssembler(fes, dirichletValues);
  auto assKD = Ikarus::makeDenseFlatAssembler(fes, dirichletValues);

  assK->bind(req, Ikarus::AffordanceCollections::elastoStatics, Ikarus::DBCOption::Reduced);
  assKD->bind(req, Ikarus::AffordanceCollections::elastoStatics, Ikarus::DBCOption::Reduced);

  const auto K  = assK->matrix();
  const auto KD = assKD->matrix();

  const auto [M, MD] = Testing::generateRandomMatrices(KD.rows());

  int nev = 10; // number of requested eigenvalues
  using Ikarus::EigenValueSolverType;

  auto partialSolver = Ikarus::PartialGeneralizedSymEigenSolver(K, M, nev);
  t.checkThrow([&]() { partialSolver.eigenvalues(); }) << testLocation();
  t.check(partialSolver.compute()) << testLocation();

  auto partialSolverDense = Ikarus::PartialGeneralizedSymEigenSolver(KD, MD, nev);
  t.require(partialSolverDense.compute()) << testLocation();

  auto solver1 = Ikarus::makeGeneralizedSymEigenSolver<EigenValueSolverType::Eigen>(KD, MD);
  t.require(solver1.compute()) << testLocation();

  auto solver2 = Ikarus::makeGeneralizedSymEigenSolver<EigenValueSolverType::Spectra>(K, M);
  t.require(solver2.compute()) << testLocation();

  auto solver3 = Ikarus::makeGeneralizedSymEigenSolver<EigenValueSolverType::Spectra>(KD, MD);
  t.require(solver3.compute()) << testLocation();

  t.subTest(testEigenVectors(solver2, solver3, K, "Testing solver2 against solver3 with K"));
  t.subTest(testEigenVectors(solver3, solver2, KD, "Testing solver3 against solver2 with KD"));
  t.subTest(testEigenVectors(solver1, solver2, KD, "Testing solver1 against solver2 with KD"));
  t.subTest(testEigenVectors(solver2, solver1, K, "Testing solver2 against solver1 with K"));

  t.subTest(testEigenValues(partialSolver, solver1, nev, "Testing partialSolver against solver1"));
  t.subTest(testEigenValues(partialSolver, solver2, nev, "Testing partialSolver against solver2"));
  t.subTest(testEigenValues(partialSolver, solver3, nev, "Testing partialSolver against solver3"));
  t.subTest(testEigenValues(solver1, solver2, {}, "Testing solver1 against solver2"));
  t.subTest(testEigenValues(solver1, solver3, {}, "Testing solver1 against solver3"));
  t.subTest(testEigenValues(solver2, solver3, {}, "Testing solver2 against solver3"));
  t.subTest(testEigenValues(partialSolver, partialSolverDense, {}, "Testing partialSolver against partialSolverDense"));

  static_assert(Ikarus::Concepts::EigenValueSolver<decltype(partialSolver)>);
  static_assert(Ikarus::Concepts::EigenValueSolver<decltype(solver1)>);
  static_assert(Ikarus::Concepts::EigenValueSolver<decltype(solver2)>);
  static_assert(Ikarus::Concepts::EigenValueSolver<decltype(solver3)>);

  auto partialIdentitySolver1 = Ikarus::makePartialIdentitySymEigenSolver(K, 10);
  auto partialIdentitySolver2 = Ikarus::makePartialIdentitySymEigenSolver(KD, 10);

  t.require(partialIdentitySolver1.compute()) << testLocation();
  t.require(partialIdentitySolver2.compute()) << testLocation();

  t.subTest(testEigenValues(partialIdentitySolver1, partialIdentitySolver2));

  auto solver1Identity = Ikarus::makeIdentitySymEigenSolver<EigenValueSolverType::Eigen>(assKD);
  t.require(solver1Identity.compute()) << testLocation();

  auto solver2Identity = Ikarus::makeIdentitySymEigenSolver<EigenValueSolverType::Spectra>(assK);
  t.require(solver2Identity.compute()) << testLocation();

  auto solver3Identity = Ikarus::makeIdentitySymEigenSolver<EigenValueSolverType::Spectra>(KD);
  t.require(solver3Identity.compute()) << testLocation();

  t.subTest(
      testEigenVectors(solver2Identity, solver3Identity, K, "Testing solver2Identity against solver3Identity with K"));
  t.subTest(testEigenVectors(solver3Identity, solver2Identity, KD,
                             "Testing solver3Identity against solver2Identity with KD"));
  t.subTest(testEigenVectors(solver1Identity, solver2Identity, KD,
                             "Testing solver1Identity against solver2Identity with KD"));
  t.subTest(
      testEigenVectors(solver2Identity, solver1Identity, K, "Testing solver2Identity against solver1Identity with K"));

  t.subTest(testEigenValues(solver1Identity, solver2Identity, {}, "Testing solver1Identity against solver2Identity"));
  t.subTest(testEigenValues(solver1Identity, solver3Identity, {}, "Testing solver1Identity against solver3Identity"));
  t.subTest(testEigenValues(solver2Identity, solver3Identity, {}, "Testing solver2Identity against solver3Identity"));

  return t;
}

auto checkThrows() {
  TestSuite t;
  using Ikarus::EigenValueSolverType;

  auto testMatrix1 = Eigen::MatrixX<double>::Random(10, 10).eval();
  auto testMatrix2 = Eigen::MatrixX<double>::Random(9, 9).eval();
  auto testMatrix3 = Eigen::MatrixX<double>::Random(9, 10).eval();

  auto sparseTestMatrix1 = Eigen::SparseMatrix<double>(testMatrix1.sparseView());
  auto sparseTestMatrix2 = Eigen::SparseMatrix<double>(testMatrix2.sparseView());
  auto sparseTestMatrix3 = Eigen::SparseMatrix<double>(testMatrix3.sparseView());

  t.checkThrow([&]() {
    Ikarus::GeneralizedSymEigenSolver<EigenValueSolverType::Spectra, Eigen::SparseMatrix<double>>(sparseTestMatrix1,
                                                                                                  sparseTestMatrix2);
  }) << testLocation();
  t.checkThrow([&]() {
    Ikarus::GeneralizedSymEigenSolver<EigenValueSolverType::Spectra, Eigen::SparseMatrix<double>>(sparseTestMatrix2,
                                                                                                  sparseTestMatrix1);
  }) << testLocation();
  t.checkThrow([&]() {
    Ikarus::GeneralizedSymEigenSolver<EigenValueSolverType::Spectra, Eigen::SparseMatrix<double>>(sparseTestMatrix3,
                                                                                                  sparseTestMatrix1);
  }) << testLocation();

  t.checkThrow([&]() {
    Ikarus::GeneralizedSymEigenSolver<EigenValueSolverType::Eigen, Eigen::MatrixX<double>>(testMatrix1, testMatrix2);
  }) << testLocation();
  t.checkThrow([&]() {
    Ikarus::GeneralizedSymEigenSolver<EigenValueSolverType::Spectra, Eigen::MatrixX<double>>(testMatrix2, testMatrix1);
  }) << testLocation();
  t.checkThrow([&]() {
    Ikarus::GeneralizedSymEigenSolver<EigenValueSolverType::Spectra, Eigen::MatrixX<double>>(testMatrix3, testMatrix1);
  }) << testLocation();

  t.checkThrow([&]() { Ikarus::PartialGeneralizedSymEigenSolver<Eigen::MatrixX<double>>(testMatrix1, testMatrix2, 7); })
      << testLocation();
  t.checkThrow([&]() { Ikarus::PartialGeneralizedSymEigenSolver<Eigen::MatrixX<double>>(testMatrix2, testMatrix1, 7); })
      << testLocation();
  t.checkThrow([&]() { Ikarus::PartialGeneralizedSymEigenSolver<Eigen::MatrixX<double>>(testMatrix3, testMatrix1, 7); })
      << testLocation();

  t.checkThrow([&]() {
    Ikarus::PartialGeneralizedSymEigenSolver<Eigen::MatrixX<double>>(testMatrix1, testMatrix1, 15);
  }) << testLocation();

  return t;
}

int main(int argc, char** argv) {
  Ikarus::init(argc, argv);
  TestSuite t;

  t.subTest(testEigenSolvers());
  t.subTest(checkThrows());

  return t.exit();
}
