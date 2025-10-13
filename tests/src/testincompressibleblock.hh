// SPDX-FileCopyrightText: 2021-2025 The Ikarus Developers ikarus@ibb.uni-stuttgart.de
// SPDX-License-Identifier: LGPL-3.0-or-later

#pragma once

#include "testcommon.hh"
#include "testhelpers.hh"

#include <chrono>

#include <dune/common/test/testsuite.hh>
#include <dune/functions/functionspacebases/basistags.hh>
#include <dune/functions/functionspacebases/boundarydofs.hh>
#include <dune/functions/functionspacebases/lagrangebasis.hh>
#include <dune/functions/functionspacebases/lagrangedgbasis.hh>
#include <dune/functions/functionspacebases/powerbasis.hh>
#include <dune/vtk/datacollectors/discontinuousdatacollector.hh>

#include <Eigen/Core>

#include <ikarus/assembler/assemblermanipulatorfuser.hh>
#include <ikarus/assembler/simpleassemblers.hh>
#include <ikarus/controlroutines/loadcontrol.hh>
#include <ikarus/finiteelements/fefactory.hh>
#include <ikarus/finiteelements/mechanics/enhancedassumedstrains.hh>
#include <ikarus/finiteelements/mechanics/materials.hh>
#include <ikarus/finiteelements/mechanics/nonlinearelastic.hh>
#include <ikarus/solver/linearsolver/linearsolver.hh>
#include <ikarus/solver/nonlinearsolver/newtonraphson.hh>
#include <ikarus/solver/nonlinearsolver/nonlinearsolverfactory.hh>
#include <ikarus/utils/basis.hh>
#include <ikarus/utils/dirichletvalues.hh>
#include <ikarus/utils/functionhelper.hh>
#include <ikarus/utils/init.hh>
#include <ikarus/utils/listener/controllogger.hh>
#include <ikarus/utils/listener/controlvtkwriter.hh>
#include <ikarus/utils/listener/nonlinearsolverlogger.hh>

using namespace Ikarus;
using Dune::TestSuite;

namespace Testing {

template <typename Grid>
auto makeGrid(double L, int nele) {
  Dune::FieldVector<double, 3> bbox       = {L, L, L};
  std::array<int, 3> elementsPerDirection = {nele, nele, nele};
  return std::make_shared<Grid>(bbox, elementsPerDirection);
}
} // namespace Testing

template <int pD, int pP, bool continuous, typename Skills>
auto incompressibelBlockTest(Skills&& elePre, std::pair<int, double> testResults, int nele, bool logToConsole = false,
                             bool writeVTK = false) {
  TestSuite t("Incompressible Block Test for Displacement Pressure Element " + Dune::className(elePre));

  spdlog::info("Simulating Incompressible Block Test with pD = " + std::to_string(pD) +
               " and pP = " + std::to_string(pP));

  constexpr double tol  = 1e-10;
  constexpr int gridDim = 3;
  using Grid            = Dune::YaspGrid<gridDim>;
  const double L        = 50;

  auto grid     = Testing::makeGrid<Grid>(L, nele);
  auto gridView = grid->leafGridView();
  using namespace Dune::Functions::BasisFactory;
  auto basis = [&]() {
    if constexpr (continuous)
      return Ikarus::makeBasis(
          gridView, composite(power<3>(lagrange<pD>(), FlatInterleaved{}), lagrange<pP>(), BlockedLexicographic{}));
    else
      return Ikarus::makeBasis(
          gridView, composite(power<3>(lagrange<pD>(), FlatInterleaved{}), lagrangeDG<pP>(), BlockedLexicographic{}));
  }();
  auto vL = [](auto& globalCoord, auto& lambda) {
    Eigen::Vector<double, 3> fext{0.0, 0.0, -9 * lambda};
    return fext;
  };

  Dune::BitSetVector<1> neumannVertices(gridView.size(2), false);
  const auto& indexSet = gridView.indexSet();

  for (auto&& vertex : vertices(gridView)) {
    auto pos       = vertex.geometry().center();
    bool isNeumann = Dune::FloatCmp::eq(pos[2], L) and (pos[1] < (L / 2) + 1e-8) and (pos[0] < (L / 2) + 1e-8);
    neumannVertices[indexSet.index(vertex)] = isNeumann;
  }

  BoundaryPatch neumannBoundary(gridView, neumannVertices);

  auto sk      = skills(elePre, neumannBoundaryLoad(&neumannBoundary, vL));
  using FEType = decltype(makeFE(basis, sk));
  std::vector<FEType> fes;

  for (auto&& ge : elements(gridView)) {
    fes.emplace_back(makeFE(basis, sk));
    fes.back().bind(ge);
  }

  Ikarus::DirichletValues dirichletValues(basis.flat());

  // fix u_x at top and X_x = 0
  dirichletValues.fixBoundaryDOFs(
      [&](auto& dirichletFlags, auto&& localIndex, auto&& localView, auto&& intersection) {
        if (Dune::FloatCmp::eq(intersection.geometry().center()[2], L, 1e-8) or
            std::abs(intersection.geometry().center()[0]) < 1e-8)
          dirichletFlags[localView.index(localIndex)] = true;
      },
      Dune::TypeTree::treePath(Dune::Indices::_0, 0));

  // fix u_y at top and X_y = 0
  dirichletValues.fixBoundaryDOFs(
      [&](auto& dirichletFlags, auto&& localIndex, auto&& localView, auto&& intersection) {
        if (Dune::FloatCmp::eq(intersection.geometry().center()[2], L, 1e-8) or
            std::abs(intersection.geometry().center()[1]) < 1e-8)
          dirichletFlags[localView.index(localIndex)] = true;
      },
      Dune::TypeTree::treePath(Dune::Indices::_0, 1));

  // fix u_z bottom
  dirichletValues.fixBoundaryDOFs(
      [&](auto& dirichletFlags, auto&& localIndex, auto&& localView, auto&& intersection) {
        if (std::abs(intersection.geometry().center()[2]) < 1e-8)
          dirichletFlags[localView.index(localIndex)] = true;
      },
      Dune::TypeTree::treePath(Dune::Indices::_0, 2));

  auto sparseFlatAssembler = makeSparseFlatAssembler(fes, dirichletValues);

  Eigen::VectorXd d;
  d.setZero(basis.flat().size());
  double lambda = 0.0;
  auto req      = typename FEType::Requirement(d, lambda);

  sparseFlatAssembler->bind(req, Ikarus::AffordanceCollections::elastoStatics, DBCOption::Full);

  auto linSolver = LinearSolver(SolverTypeTag::sd_SparseQR);
  AffordanceCollection elastoStaticsNoScalar(VectorAffordance::forces, MatrixAffordance::stiffness);
  auto nonOp =
      DifferentiableFunctionFactory::op(sparseFlatAssembler, elastoStaticsNoScalar, sparseFlatAssembler->dBCOption());

  auto nrConfig =
      Ikarus::NewtonRaphsonConfig<decltype(linSolver)>{.parameters = {.tol = tol}, .linearSolver = linSolver};
  auto nrFactory = NonlinearSolverFactory(nrConfig);
  auto nr        = nrFactory.create(sparseFlatAssembler);

  auto lc = ControlRoutineFactory::create(LoadControlConfig{3, 0.0, 1.0}, nr, sparseFlatAssembler);

  auto nonLinearSolverLogger = NonLinearSolverLogger();
  auto controlLogger         = ControlLogger();

  if (logToConsole) {
    nonLinearSolverLogger.subscribeTo(lc.nonLinearSolver());
    controlLogger.subscribeTo(lc);
  }

  auto start                            = std::chrono::high_resolution_clock::now();
  const auto controlInfo                = lc.run(req);
  auto end                              = std::chrono::high_resolution_clock::now();
  std::chrono::duration<double> elapsed = end - start;
  std::cout << "Elapsed time: " << elapsed.count() << " seconds (nele = " << nele << ")" << std::endl;

  d      = req.globalSolution();
  lambda = req.parameter();

  if (writeVTK) {
    Ikarus::Vtk::Writer writer(sparseFlatAssembler);
    auto dBasis = Dune::Functions::subspaceBasis(basis.flat(), Dune::Indices::_0);
    auto pBasis = Dune::Functions::subspaceBasis(basis.flat(), Dune::Indices::_1);

    writer.addInterpolation(d, dBasis, "displacement");
    writer.addInterpolation(d, pBasis, "pressure");

    writer.write("block_" + std::to_string(nele) + "_" + std::to_string(pD) + "_" + std::to_string(continuous));
  }

  double expectedLambda                            = 1.0;
  const auto [expectedIterations, expectedMaxDisp] = testResults;

  t.check(controlInfo.success);
  checkScalars(t, controlInfo.totalIterations, expectedIterations, " Number of iterations");
  const auto maxDisp = std::ranges::max(d.cwiseAbs());

  checkScalars(t, maxDisp, expectedMaxDisp, " Max. displacement", tol);
  checkScalars(t, lambda, expectedLambda, " Load factor", tol);

  return t;
}
