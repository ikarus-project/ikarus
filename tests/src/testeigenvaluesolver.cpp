// SPDX-FileCopyrightText: 2021-2024 The Ikarus Developers mueller@ibb.uni-stuttgart.de
// SPDX-License-Identifier: LGPL-3.0-or-later
#include <config.h>

#include "tests/src/testhelpers.hh"

#include <vector>

#include <dune/common/float_cmp.hh>
#include <dune/common/test/testsuite.hh>
#include <dune/fufem/dunepython.hh>
#include <dune/functions/functionspacebases/basistags.hh>
#include <dune/functions/functionspacebases/boundarydofs.hh>
#include <dune/functions/functionspacebases/interpolate.hh>
#include <dune/functions/functionspacebases/lagrangebasis.hh>
#include <dune/functions/functionspacebases/powerbasis.hh>
#include <dune/grid/uggrid.hh>
#include <dune/grid/yaspgrid.hh>
#include <dune/localfefunctions/eigenDuneTransformations.hh>
#include <dune/vtk/datacollectors/discontinuousdatacollector.hh>

#include <ikarus/assembler/simpleassemblers.hh>
#include <ikarus/finiteelements/fefactory.hh>
#include <ikarus/finiteelements/mechanics/linearelastic.hh>
#include <ikarus/finiteelements/mechanics/loads/volume.hh>
#include <ikarus/finiteelements/mixin.hh>
#include <ikarus/io/resultfunction.hh>
#include <ikarus/io/vtkwriter.hh>
#include <ikarus/solver/eigenvaluesolver/generaleigensolver.hh>
#include <ikarus/utils/basis.hh>
#include <ikarus/utils/dirichletvalues.hh>
#include <ikarus/utils/dynamics/dynamics.hh>
#include <ikarus/utils/init.hh>

using Dune::TestSuite;

template <Ikarus::Concepts::EigenValueSolver SOL1, Ikarus::Concepts::EigenValueSolver SOL2>
auto testSolversAgainstEachOther(const SOL1& solver1, const SOL2& solver2,
                                 std::optional<Eigen::Index> nev_ = std::nullopt) {
  TestSuite t("Test " + Dune::className<SOL1>() + ", " + Dune::className<SOL2>());
  auto eigenvalues1  = solver1.eigenvalues();
  auto eigenvectors1 = solver1.normalizedEigenvectors().cwiseAbs().eval();

  auto eigenvalues2  = solver2.eigenvalues();
  auto eigenvectors2 = solver2.normalizedEigenvectors().cwiseAbs().eval();

  if (not nev_.has_value()) {
    t.check(isApproxSame(eigenvalues1, eigenvalues2, 1e-10)) << testLocation();
    t.check(isApproxSame(eigenvectors1, eigenvectors2, 1e-10)) << testLocation();
  } else {
    t.check(isApproxSame(eigenvalues1.head(*nev_).eval(), eigenvalues2.head(*nev_).eval(), 1e-10))
        << testLocation() << "\n"
        << eigenvalues1.head(*nev_).eval() << "\n\n"
        << eigenvalues2.head(*nev_).eval();
    ;
    t.check(isApproxSame(eigenvectors1.leftCols(*nev_).eval(), eigenvectors2.leftCols(*nev_).eval(), 1e-10))
        << testLocation() << "\n"
        << eigenvectors1.leftCols(*nev_).eval() << "\n\n"
        << eigenvectors2.leftCols(*nev_).eval();
  }

  return t;
}

static auto testRealWorldProblem() {
  TestSuite t("EigenvalueTest");
  using Grid = Dune::YaspGrid<2>;

  const double Lx                         = 4.0;
  const double Ly                         = 4.0;
  Dune::FieldVector<double, 2> bbox       = {Lx, Ly};
  std::array<int, 2> elementsPerDirection = {4, 4};
  auto grid                               = std::make_shared<Grid>(bbox, elementsPerDirection);

  auto gridView = grid->leafGridView();

  using namespace Dune::Functions::BasisFactory;
  auto basis = Ikarus::makeBasis(gridView, power<2>(lagrange<1>()));

  auto vL      = []([[maybe_unused]] auto& globalCoord, auto& lamb) { return Eigen::Vector2d{0, -1}; };
  auto linMat  = Ikarus::LinearElasticity(Ikarus::toLamesFirstParameterAndShearModulus({.emodul = 100, .nu = 0.2}));
  auto skills_ = Ikarus::skills(Ikarus::linearElastic(Ikarus::planeStress(linMat)), Ikarus::volumeLoad<2>(vL));

  using FEType = decltype(Ikarus::makeFE(basis, skills_));
  std::vector<FEType> fes;
  for (auto&& element : elements(gridView)) {
    fes.emplace_back(Ikarus::makeFE(basis, skills_));
    fes.back().bind(element);
  }

  auto& fe = fes[0];
  auto req = FEType::Requirement();

  typename FEType::Requirement::SolutionVectorType d;
  d.setZero(basis.flat().size());
  req.insertGlobalSolution(d);

  auto dirichletValues = Ikarus::DirichletValues(basis.flat());
  dirichletValues.fixBoundaryDOFs([](auto& dirichFlags, auto&& indexGlobal) { dirichFlags[indexGlobal] = true; });

  auto assM = Ikarus::makeSparseFlatAssembler(fes, dirichletValues);
  auto assK = Ikarus::makeSparseFlatAssembler(fes, dirichletValues);

  auto assMD = Ikarus::makeDenseFlatAssembler(fes, dirichletValues);
  auto assKD = Ikarus::makeDenseFlatAssembler(fes, dirichletValues);

  assM->bind(req, Ikarus::AffordanceCollections::dynamics, Ikarus::DBCOption::Reduced);
  assK->bind(req, Ikarus::AffordanceCollections::elastoStatics, Ikarus::DBCOption::Reduced);

  assMD->bind(req, Ikarus::AffordanceCollections::dynamics, Ikarus::DBCOption::Reduced);
  assKD->bind(req, Ikarus::AffordanceCollections::elastoStatics, Ikarus::DBCOption::Reduced);

  int nev = 10; // number of requested eigenvalues
  using Ikarus::EigenSolverTypeTag;

  auto solver = Ikarus::PartialGeneralSymEigenSolver(assK, assM, nev);
  t.checkThrow([&]() { solver.eigenvalues(); }) << testLocation();
  bool success = solver.compute();
  t.check(success) << testLocation();

  auto solver3 = Ikarus::makeGeneralSymEigenSolver<EigenSolverTypeTag::Eigen>(assKD, assMD);
  solver3.compute();

  auto solver4 = Ikarus::makeGeneralSymEigenSolver<EigenSolverTypeTag::Spectra>(assK, assM);
  solver4.compute();

  auto solver5 = Ikarus::makeGeneralSymEigenSolver<EigenSolverTypeTag::Spectra>(assKD, assMD);
  solver5.compute();

  t.subTest(testSolversAgainstEachOther(solver, solver3, nev));
  t.subTest(testSolversAgainstEachOther(solver, solver4, nev));
  t.subTest(testSolversAgainstEachOther(solver, solver5, nev));
  t.subTest(testSolversAgainstEachOther(solver3, solver4));
  t.subTest(testSolversAgainstEachOther(solver4, solver5));

  static_assert(Ikarus::Concepts::EigenValueSolver<decltype(solver)>);
  static_assert(Ikarus::Concepts::EigenValueSolver<decltype(solver3)>);
  static_assert(Ikarus::Concepts::EigenValueSolver<decltype(solver4)>);
  static_assert(Ikarus::Concepts::EigenValueSolver<decltype(solver5)>);

  return t;
}

int main(int argc, char** argv) {
  Ikarus::init(argc, argv);
  TestSuite t;

  t.subTest(testRealWorldProblem());
  return t.exit();
}
