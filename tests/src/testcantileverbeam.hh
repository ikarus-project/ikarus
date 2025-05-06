// SPDX-FileCopyrightText: 2021-2025 The Ikarus Developers ikarus@ibb.uni-stuttgart.de
// SPDX-License-Identifier: LGPL-3.0-or-later

#pragma once

#include "testcommon.hh"
#include "testhelpers.hh"

#include <dune/common/test/testsuite.hh>
#include <dune/functions/functionspacebases/basistags.hh>
#include <dune/functions/functionspacebases/boundarydofs.hh>
#include <dune/functions/functionspacebases/lagrangebasis.hh>
#include <dune/functions/functionspacebases/powerbasis.hh>

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

template <typename Grid, int gridDim>
auto makeGrid(double L, double h) {
  if constexpr (gridDim == 2) {
    Dune::FieldVector<double, 2> bbox             = {L, h};
    std::array<int, gridDim> elementsPerDirection = {10, 1};
    return std::make_shared<Grid>(bbox, elementsPerDirection);
  } else {
    Dune::FieldVector<double, 3> bbox             = {L, h, h};
    std::array<int, gridDim> elementsPerDirection = {10, 1, 1};
    return std::make_shared<Grid>(bbox, elementsPerDirection);
  }
}
template <int gridDim>
requires(gridDim == 2)
auto makePointLoad(const auto& basis, double L, double h) {
  Dune::FieldVector<double, gridDim> topRightPos{L, h};
  const auto globalIndices = utils::globalIndexFromGlobalPosition(basis.flat(), topRightPos);
  auto pointLoad           = [=](const auto&, const auto& par, auto, auto, Eigen::VectorXd& vec) -> void {
    auto loadFactor = par.parameter();
    vec[globalIndices[1]] -= -loadFactor * 1.0;
  };
  return pointLoad;
}

template <int gridDim>
requires(gridDim == 3)
auto makePointLoad(const auto& basis, double L, double h) {
  Dune::FieldVector<double, gridDim> topRightPos1{L, h, h};
  Dune::FieldVector<double, gridDim> topRightPos2{L, h, 0};
  const auto globalIndices1 = utils::globalIndexFromGlobalPosition(basis.flat(), topRightPos1);
  const auto globalIndices2 = utils::globalIndexFromGlobalPosition(basis.flat(), topRightPos2);

  auto pointLoad = [=](const auto&, const auto& par, auto, auto, Eigen::VectorXd& vec) -> void {
    auto loadFactor = par.parameter();
    vec[globalIndices1[1]] -= -loadFactor * 1.0;
    vec[globalIndices2[1]] -= -loadFactor * 1.0;
  };
  return pointLoad;
}
} // namespace Testing

template <int dim, typename MAT, typename Skills>
requires(dim == 2 or dim == 3)
auto cantileverBeamTest(const MAT& material, Skills&& additionalSkills, std::pair<int, double> testResults,
                        bool logToConsole = false, bool writeVTK = false) {
  if constexpr (dim == 2)
    static_assert(MAT::isReduced,
                  "cantileverBeamTest (3d) is only valid for a reduced material (planeStress or planeStrain).");
  if constexpr (dim == 3)
    static_assert(not MAT::isReduced, "cantileverBeamTest (3d) is only valid for 3d material");

  TestSuite t("Cantilever Beam for Nonlinear element with additional skills " + Dune::className(additionalSkills) +
              " and material type as " + Dune::className<MAT>());

  constexpr double tol  = 1e-10;
  constexpr int gridDim = dim;
  using Grid            = Dune::YaspGrid<gridDim>;
  const double L        = 10;
  const double h        = 2;

  auto grid     = Testing::makeGrid<Grid, gridDim>(L, h);
  auto gridView = grid->leafGridView();
  using namespace Dune::Functions::BasisFactory;
  auto basis = Ikarus::makeBasis(gridView, power<gridDim>(lagrange<1>(), FlatInterleaved()));

  auto sk      = merge(skills(nonLinearElastic(material)), std::move(additionalSkills));
  using FEType = decltype(makeFE(basis, sk));
  std::vector<FEType> fes;

  for (auto&& ge : elements(gridView)) {
    fes.emplace_back(makeFE(basis, sk));
    fes.back().bind(ge);
  }

  Ikarus::DirichletValues dirichletValues(basis.flat());

  // fix left edge (x=0)
  dirichletValues.fixBoundaryDOFs([&](auto& dirichletFlags, auto&& localIndex, auto&& localView, auto&& intersection) {
    if (std::abs(intersection.geometry().center()[0]) < 1e-8)
      dirichletFlags[localView.index(localIndex)] = true;
  });

  auto sparseFlatAssembler = makeSparseFlatAssembler(fes, dirichletValues);
  auto sparseAssemblerAM   = makeAssemblerManipulator(*sparseFlatAssembler);

  Eigen::VectorXd d;
  d.setZero(basis.flat().size());
  double lambda = 0.0;
  auto req      = typename FEType::Requirement(d, lambda);

  sparseAssemblerAM->bind(req, Ikarus::AffordanceCollections::elastoStatics, DBCOption::Full);
  // Apply constant point load at the top right corner
  sparseAssemblerAM->bind(Testing::makePointLoad<gridDim>(basis, L, h));

  auto linSolver = LinearSolver(SolverTypeTag::sd_UmfPackLU);
  AffordanceCollection elastoStaticsNoScalar(VectorAffordance::forces, MatrixAffordance::stiffness);
  auto nonOp =
      DifferentiableFunctionFactory::op(sparseAssemblerAM, elastoStaticsNoScalar, sparseAssemblerAM->dBCOption());

  auto nrConfig =
      Ikarus::NewtonRaphsonConfig<decltype(linSolver)>{.parameters = {.tol = tol}, .linearSolver = linSolver};
  NonlinearSolverFactory nrFactory(nrConfig);
  auto nr = nrFactory.create(sparseAssemblerAM);

  // Only when creating the control routine via the Factory, the elements get registered for correction update
  // automatically.
  auto lc = ControlRoutineFactory::create(LoadControlConfig{20, 0.0, 1.0}, nr, sparseFlatAssembler);

  auto nonLinearSolverObserver = NonLinearSolverLogger();
  auto controlLogger           = ControlLogger();
  auto vtkWriter = ControlSubsamplingVertexVTKWriter<std::remove_cvref_t<decltype(basis.flat())>>(basis.flat(), d);

  if (logToConsole) {
    nonLinearSolverObserver.subscribeTo(lc.nonLinearSolver());
    controlLogger.subscribeTo(lc);
  }
  if (writeVTK) {
    vtkWriter.setFileNamePrefix("CantileverNonlinearEAS");
    vtkWriter.setFieldInfo("Displacement", Dune::VTK::FieldInfo::Type::vector, 2);
    vtkWriter.subscribeTo(lc);
  }

  const auto controlInfo = lc.run(req);
  d                      = req.globalSolution();
  lambda                 = req.parameter();

  double expectedLambda                            = 1.0;
  const auto [expectedIterations, expectedMaxDisp] = testResults;

  t.check(controlInfo.success);
  checkScalars(t, controlInfo.totalIterations, expectedIterations, " Number of iterations");
  const auto maxDisp = std::ranges::max(d.cwiseAbs());

  checkScalars(t, maxDisp, expectedMaxDisp, " Max. displacement", tol);
  checkScalars(t, lambda, expectedLambda, " Load factor", tol);

  return t;
}
