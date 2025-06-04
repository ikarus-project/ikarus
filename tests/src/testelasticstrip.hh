// SPDX-FileCopyrightText: 2021-2025 The Ikarus Developers ikarus@ibb.uni-stuttgart.de
// SPDX-License-Identifier: LGPL-3.0-or-later

#pragma once

#include "testhelpers.hh"

#include <dune/common/test/testsuite.hh>
#include <dune/functions/functionspacebases/basistags.hh>
#include <dune/functions/functionspacebases/boundarydofs.hh>
#include <dune/functions/functionspacebases/lagrangebasis.hh>
#include <dune/functions/functionspacebases/powerbasis.hh>
#include <dune/grid/yaspgrid.hh>

#include <Eigen/Core>

#include <ikarus/assembler/assemblermanipulatorfuser.hh>
#include <ikarus/assembler/simpleassemblers.hh>
#include <ikarus/controlroutines/loadcontrol.hh>
#include <ikarus/controlroutines/pathfollowing.hh>
#include <ikarus/finiteelements/fefactory.hh>
#include <ikarus/finiteelements/mechanics/enhancedassumedstrains.hh>
#include <ikarus/finiteelements/mechanics/materials.hh>
#include <ikarus/finiteelements/mechanics/nonlinearelastic.hh>
#include <ikarus/solver/linearsolver/linearsolver.hh>
#include <ikarus/solver/nonlinearsolver/newtonraphson.hh>
#include <ikarus/solver/nonlinearsolver/newtonraphsonwithscalarsubsidiaryfunction.hh>
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

/** Adapted from Pfefferkorn et al. 2021 (https://doi.org/10.1002/nme.6605) */
template <bool useArcLength, typename MAT, typename Skills>
auto elasticStripTest(DBCOption dbcOption, const MAT& material, Skills&& additionalSkills, int loadSteps,
                      std::tuple<int, double, double> testResults, int order = 1, bool logToConsole = false,
                      bool writeVTK = false) {
  static_assert(MAT::isReduced, "elasticStripTest is only valid for a reduced material (planeStress or planeStrain).");

  TestSuite t("Elastic Strip Test for nonlinear solid element with additional skills " +
              Dune::className(additionalSkills) + ", material type as " + Dune::className<MAT>() +
              ", dbcOption = " + toString(dbcOption) + " and order = " + std::to_string(order));

  std::cout << "Started: " << t.name() << std::endl;

  constexpr double tol  = 1e-8;
  constexpr int gridDim = 2;
  using Grid            = Dune::YaspGrid<gridDim>;
  const double L        = 10;
  const double h        = L;

  Dune::FieldVector<double, 2> bbox       = {L, h};
  std::array<int, 2> elementsPerDirection = {10, 10};
  auto grid                               = std::make_shared<Grid>(bbox, elementsPerDirection);
  auto gridView                           = grid->leafGridView();
  using namespace Dune::Functions::BasisFactory;
  auto basis = Ikarus::makeBasis(gridView, power<gridDim>(lagrange(order), FlatInterleaved()));

  auto sk      = merge(skills(nonLinearElastic(material)), std::move(additionalSkills));
  using FEType = decltype(makeFE(basis, sk));
  std::vector<FEType> fes;

  for (auto&& ge : elements(gridView)) {
    fes.emplace_back(makeFE(basis, sk));
    fes.back().bind(ge);
  }

  Ikarus::DirichletValues dirichletValues(basis.flat());

  // fix ux and uy at left edge (x=0)
  dirichletValues.fixBoundaryDOFs([&](auto& dirichletFlags, auto&& localIndex, auto&& localView, auto&& intersection) {
    if (std::abs(intersection.geometry().center()[0]) < 1e-8)
      dirichletFlags[localView.index(localIndex)] = true;
  });

  /// fix uy at the right edge (x=L)
  dirichletValues.fixDOFs([&](auto& basis_, auto&& dirichletFlags) {
    Dune::Functions::forEachBoundaryDOF(Dune::Functions::subspaceBasis(basis_, 1),
                                        [&](auto&& localIndex, auto&& localView, auto&& intersection) {
                                          if (std::abs(intersection.geometry().center()[0] - L) < 1e-8)
                                            dirichletFlags[localView.index(localIndex)] = true;
                                        });
  });

  int numBoundaryNodesX = order * elementsPerDirection[0] + 1;
  int numBoundaryNodesY = order * elementsPerDirection[1] + 1;
  t.check(dirichletValues.fixedDOFsize() == 3 * numBoundaryNodesY,
          "Number of constrained DOFs is not equal to " + std::to_string(3 * numBoundaryNodesY));

  // Inhomogeneous Boundary Conditions
  auto inhomogeneousDisplacement = [L, h]<typename T>(const auto& globalCoord, const T& lambda) {
    Eigen::Vector<T, 2> localInhomogeneous;
    if (Dune::FloatCmp::eq(globalCoord[0], L)) {
      localInhomogeneous[0] = 10.0 * lambda;
      localInhomogeneous[1] = 0;
    } else
      localInhomogeneous.setZero();
    return localInhomogeneous;
  };

  dirichletValues.storeInhomogeneousBoundaryCondition(inhomogeneousDisplacement);

  t.check(dirichletValues.fixedDOFsize() == 4 * numBoundaryNodesY,
          "Number of constrained DOFs is not equal to " + std::to_string(4 * numBoundaryNodesY));

  auto sparseFlatAssembler = makeSparseFlatAssembler(fes, dirichletValues);

  auto req     = typename FEType::Requirement(basis);
  auto& d      = req.globalSolution();
  auto& lambda = req.parameter();

  sparseFlatAssembler->bind(req, Ikarus::AffordanceCollections::elastoStatics, dbcOption);

  auto linSolver = LinearSolver(SolverTypeTag::sd_UmfPackLU);
  AffordanceCollection elastoStaticsNoScalar(VectorAffordance::forces, MatrixAffordance::stiffness);
  auto nonOp =
      DifferentiableFunctionFactory::op(sparseFlatAssembler, elastoStaticsNoScalar, sparseFlatAssembler->dBCOption());

  // Only when creating the control routine via the Factory, the elements get registered for correction update
  // automatically.
  auto controlRoutine = [&]() -> decltype(auto) {
    if constexpr (useArcLength) { /// Create ArcLengthControl
      auto nrConfig = Ikarus::NewtonRaphsonWithSubsidiaryFunctionConfig<decltype(linSolver)>{.parameters = {.tol = tol},
                                                                                             .linearSolver = linSolver};
      auto nrFactory = NonlinearSolverFactory(nrConfig).withIDBCForceFunction(sparseFlatAssembler);
      auto nr        = nrFactory.create(sparseFlatAssembler);
      auto pft       = Ikarus::ArcLength{}; // Type of path following technique
      auto alc =
          ControlRoutineFactory::create(PathFollowingConfig{loadSteps, 45.0 / loadSteps, pft}, nr, sparseFlatAssembler);
      return alc;
    } else { /// Create loadcontrol
      auto nrConfig =
          Ikarus::NewtonRaphsonConfig<decltype(linSolver)>{.parameters = {.tol = tol}, .linearSolver = linSolver};
      auto nrFactory = NonlinearSolverFactory(nrConfig).withIDBCForceFunction(sparseFlatAssembler);
      auto nr        = nrFactory.create(sparseFlatAssembler);
      auto lc        = ControlRoutineFactory::create(LoadControlConfig{loadSteps, 0.0, 1.0}, nr, sparseFlatAssembler);
      return lc;
    }
  };

  auto cr = controlRoutine();

  auto nonLinearSolverLogger = NonLinearSolverLogger();
  auto controlLogger         = ControlLogger();
  auto vtkWriter             = ControlSubsamplingVertexVTKWriter(basis.flat());

  if (logToConsole) {
    nonLinearSolverLogger.subscribeTo(cr.nonLinearSolver());
    controlLogger.subscribeTo(cr);
  }
  if (writeVTK) {
    vtkWriter.setFileNamePrefix("ElasticStrip");
    vtkWriter.setFieldInfo("Displacement", Dune::VTK::FieldInfo::Type::vector, 2);
    vtkWriter.subscribeTo(cr);
  }

  const auto controlInfo = cr.run(req);
  d                      = req.globalSolution();
  lambda                 = req.parameter();

  // check that inhomogeneous boundary conditions are correctly applied in the final state
  Eigen::VectorXd inhomogeneousDisplacementExpected(basis.flat().dimension());
  dirichletValues.evaluateInhomogeneousBoundaryCondition(inhomogeneousDisplacementExpected, lambda);
  for (int i = 0; i < basis.flat().dimension(); ++i)
    if (std::abs(inhomogeneousDisplacementExpected[i]) > 1e-8)
      checkScalars(t, d[i], inhomogeneousDisplacementExpected[i],
                   " Inhomogeneous boundary condition not correctly applied at i = " + std::to_string(i), tol);

  if constexpr (FEType::hasEAS())
    for (const auto& fe : fes)
      t.check(not fe.internalVariable().isZero(),
              "Internal variables are zero and have not been updated by the nonlinear solver.");

  const auto [expectedIterations, expectedMaxDisp, expectedLambda] = testResults;

  t.check(controlInfo.success, "Failed to converge.");
  checkScalars(t, controlInfo.totalIterations, expectedIterations, " Number of iterations");

  Dune::FieldVector<double, gridDim> bottomCenterPos{L / 2.0, 0.0};
  const auto globalIndicesAtBottomCenterPos = utils::globalIndexFromGlobalPosition(basis.flat(), bottomCenterPos);
  const auto maxVerticalDisp                = d[globalIndicesAtBottomCenterPos[1]];

  checkScalars(t, maxVerticalDisp, expectedMaxDisp, " Max. vertical displacement", tol);
  checkScalars(t, lambda, expectedLambda, " Load factor", tol);

  return t;
}
