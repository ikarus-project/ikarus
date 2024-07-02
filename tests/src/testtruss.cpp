// SPDX-FileCopyrightText: 2021-2024 The Ikarus Developers mueller@ibb.uni-stuttgart.de
// SPDX-License-Identifier: LGPL-3.0-or-later

#include <config.h>

#include "checkfebyautodiff.hh"
#include "testcommon.hh"
#include "testhelpers.hh"

#include <matplot/matplot.h>

#include <dune/foamgrid/foamgrid.hh>
#include <dune/functions/functionspacebases/lagrangebasis.hh>
#include <dune/functions/functionspacebases/powerbasis.hh>
#include <dune/functions/gridfunctions/discreteglobalbasisfunction.hh>
#include <dune/grid/io/file/vtk/vtkwriter.hh>

#include <Eigen/Core>

#include <ikarus/assembler/simpleassemblers.hh>
#include <ikarus/controlroutines/loadcontrol.hh>
#include <ikarus/finiteelements/fefactory.hh>
#include <ikarus/finiteelements/mechanics/truss.hh>
#include <ikarus/solver/linearsolver/linearsolver.hh>
#include <ikarus/solver/nonlinearsolver/newtonraphson.hh>
#include <ikarus/solver/nonlinearsolver/nonlinearsolverfactory.hh>
#include <ikarus/utils/basis.hh>
#include <ikarus/utils/dirichletvalues.hh>
#include <ikarus/utils/init.hh>
#include <ikarus/utils/observer/controlvtkwriter.hh>
#include <ikarus/utils/observer/genericobserver.hh>
#include <ikarus/utils/observer/nonlinearsolverlogger.hh>

using namespace Ikarus;
using Dune::TestSuite;

static auto vonMisesTrussTest() {
  TestSuite t("vonMisesTrussTest");
  /// Construct grid
  Dune::GridFactory<Dune::FoamGrid<1, 2, double>> gridFactory;
  constexpr double h = 1.0;
  constexpr double L = 10.0;
  gridFactory.insertVertex({0, 0});
  gridFactory.insertVertex({L, h});
  gridFactory.insertVertex({2 * L, 0});
  gridFactory.insertElement(Dune::GeometryTypes::line, {0, 1});
  gridFactory.insertElement(Dune::GeometryTypes::line, {1, 2});
  auto grid     = gridFactory.createGrid();
  auto gridView = grid->leafGridView();

  /// Construct basis
  using namespace Dune::Functions::BasisFactory;
  auto basis = Ikarus::makeBasis(gridView, power<2>(lagrange<1>()));

  /// Create finite elements
  const double E = 30000;
  const double A = 0.1;
  auto sk        = skills(truss(E, A));
  using FEType   = decltype(makeFE(basis, sk));
  std::vector<FEType> fes;
  for (auto&& ge : elements(gridView)) {
    fes.emplace_back(makeFE(basis, sk));
    fes.back().bind(ge);
  }

  /// Collect dirichlet nodes
  auto basisP = std::make_shared<const decltype(basis)>(basis);
  DirichletValues dirichletValues(basisP->flat());
  dirichletValues.fixBoundaryDOFs(
      [&](auto& dirichletFlags, auto&& globalIndex) { dirichletFlags[globalIndex] = true; });
  t.check(dirichletValues.fixedDOFsize() == 4, "Number of fixed DOFs is not equal to 4");

  /// Create assembler
  auto denseFlatAssembler = makeDenseFlatAssembler(fes, dirichletValues);

  /// Create non-linear operator
  double lambda = 0.0;
  Eigen::VectorXd d;
  d.setZero(basis.flat().size());

  auto req = FEType::Requirement();
  req.insertGlobalSolution(d).insertParameter(lambda);
  denseFlatAssembler->bind(req, AffordanceCollections::elastoStatics, DBCOption::Full);

  /// Choose linear solver
  auto linSolver = LinearSolver(SolverTypeTag::d_LDLT);

  NewtonRaphsonConfig<decltype(linSolver)> nrConfig{.linearSolver = linSolver};

  NonlinearSolverFactory nrFactory(nrConfig);
  auto nr = nrFactory.create(denseFlatAssembler);

  /// Create Observer to write information of the non-linear solver
  auto nonLinearSolverObserver = std::make_shared<NonLinearSolverLogger>();

  constexpr int loadSteps = 10;

  Eigen::Matrix3Xd lambdaAndDisp;
  lambdaAndDisp.setZero(Eigen::NoChange, loadSteps + 1);
  /// Create Observer which executes when control routines messages
  auto lvkObserver =
      std::make_shared<GenericObserver<ControlMessages>>(ControlMessages::SOLUTION_CHANGED, [&](int step) {
        lambdaAndDisp(0, step) = lambda; // load factor
        lambdaAndDisp(1, step) = d[2];   // horizontal displacement at center node
        lambdaAndDisp(2, step) = d[3];   // vertical displacement at center node
      });

  /// Create Observer which writes vtk files when control routines messages
  auto vtkWriter = std::make_shared<ControlSubsamplingVertexVTKWriter<std::remove_cvref_t<decltype(basis.flat())>>>(
      basis.flat(), d, 2);
  vtkWriter->setFieldInfo("displacement", Dune::VTK::FieldInfo::Type::vector, 2);
  vtkWriter->setFileNamePrefix("vonMisesTruss");

  /// Create loadcontrol
  auto lc = LoadControl(nr, loadSteps, {0, 0.5});
  lc.nonlinearSolver().subscribeAll(nonLinearSolverObserver);
  lc.subscribeAll({vtkWriter, lvkObserver});

  /// Execute!
  auto controlInfo = lc.run();
  t.check(controlInfo.success == true) << "Load control failed to converge";

  Eigen::VectorXd lambdaVec = lambdaAndDisp.row(0);
  Eigen::VectorXd uVec      = lambdaAndDisp.row(1);
  Eigen::VectorXd vVec      = -lambdaAndDisp.row(2);

  // return the load factor as a function of the vertical displacement
  auto analyticalLoadDisplacementCurve = [&](auto& v) {
    const double Ltruss = std::sqrt(h * h + L * L);
    return 2.0 * E * A * Dune::power(h, 3) / Dune::power(Ltruss, 3) *
           (v / h - 1.5 * Dune::power(v / h, 2) + 0.5 * Dune::power(v / h, 3));
  };

  constexpr double tol = 1e-14;

  for (std::size_t i = 0; i < lambdaAndDisp.cols(); ++i) {
    checkScalars(t, lambdaVec[i], analyticalLoadDisplacementCurve(vVec[i]), " Incorrect load factor", tol);
    checkScalars(t, uVec[i], 0.0, " Incorrect horizontal displacement", tol);
  }

  double expectedAxialForce = -2.785092479363964262;

  // due to the symmetry of the problem, same axial forces are expected in both elements
  for (const auto& fe : fes) {
    const auto N = fe.calculateAt<ResultTypes::axialForce>(req, {0.5, 0.5}).asVec()[0];
    checkScalars(t, N, expectedAxialForce, " Incorrect Axial force", tol);
  }

  return t;
}

auto singleElementTest() {
  TestSuite t("Truss autodiff");
  using namespace Dune::Functions::BasisFactory;

  {
    Dune::GridFactory<Dune::FoamGrid<1, 2, double>> gridFactory;
    const double h = 1.7;
    const double L = 2.0;
    gridFactory.insertVertex({0, 0});
    gridFactory.insertVertex({L, h});
    gridFactory.insertVertex({2 * L, 0});
    gridFactory.insertElement(Dune::GeometryTypes::line, {0, 1});
    gridFactory.insertElement(Dune::GeometryTypes::line, {1, 2});
    auto grid     = gridFactory.createGrid();
    auto gridView = grid->leafGridView();

    t.subTest(checkFESByAutoDiff(gridView, power<2>(lagrange<1>()), skills(truss(1000.0, 0.0956)),
                                 AffordanceCollections::elastoStatics));
  }
  return t;
}

int main(int argc, char** argv) {
  Ikarus::init(argc, argv);
  TestSuite t("TrussTest");
  t.subTest(singleElementTest());
  t.subTest(vonMisesTrussTest());
  return t.exit();
}
