// SPDX-FileCopyrightText: 2021-2023 The Ikarus Developers mueller@ibb.uni-stuttgart.de
// SPDX-License-Identifier: LGPL-3.0-or-later

#include <config.h>

#include "testcommon.hh"
#include "testhelpers.hh"

#include <dune/common/test/testsuite.hh>
#include <dune/functions/functionspacebases/boundarydofs.hh>
#include <dune/functions/functionspacebases/powerbasis.hh>
#include <dune/functions/functionspacebases/subspacebasis.hh>
#include <dune/iga/nurbsbasis.hh>

#include <spdlog/spdlog.h>

#include <Eigen/Core>

#include <ikarus/assembler/simpleassemblers.hh>
#include <ikarus/controlroutines/pathfollowingtechnique.hh>
#include <ikarus/finiteelements/mechanics/kirchhoffloveshell.hh>
#include <ikarus/io/resultfunction.hh>
#include <ikarus/linearalgebra/dirichletvalues.hh>
#include <ikarus/linearalgebra/nonlinearoperator.hh>
#include <ikarus/solver/nonlinearsolver/newtonraphson.hh>
#include <ikarus/utils/algorithms.hh>
#include <ikarus/utils/basis.hh>
#include <ikarus/utils/init.hh>
#include <ikarus/utils/observer/controllogger.hh>
#include <ikarus/utils/observer/controlvtkwriter.hh>
#include <ikarus/utils/observer/nonlinearsolverlogger.hh>

using Dune::TestSuite;

template <typename PathFollowingType>
auto KLShellAndAdaptiveStepSizing(const PathFollowingType& pft, const std::vector<std::vector<int>>& expectedIterations,
                                  const std::vector<std::vector<double>>& expectedResults, const int targetIterations,
                                  const double stepSize) {
  TestSuite t("KLShellAndAdaptiveStepSizing --> " + pft.name);
  constexpr auto dimWorld        = 3;
  const std::array<int, 2> order = {2, 1};

  const std::array<std::vector<double>, 2> knotSpans = {{{0, 0, 0, 1, 1, 1}, {0, 0, 1, 1}}};

  using ControlPoint = Dune::IGA::NURBSPatchData<2, dimWorld>::ControlPointType;

  const std::vector<std::vector<ControlPoint>> controlPoints
      = {{{.p = {-0.5, 0, 0.0}, .w = 10}, {.p = {-0.5, 1.0, 0.0}, .w = 1}},
         {{.p = {0, 0, 0.2}, .w = 0.766044}, {.p = {0, 1.0, 0.2}, .w = 0.766044}},
         {{.p = {0.5, 0, 0.0}, .w = 1}, {.p = {0.5, 1.0, 0.0}, .w = 1}}};

  std::array<int, 2> dimsize = {(int)(controlPoints.size()), (int)(controlPoints[0].size())};

  auto controlNet = Dune::IGA::NURBSPatchData<2, dimWorld>::ControlPointNetType(dimsize, controlPoints);
  using Grid      = Dune::IGA::NURBSGrid<2, dimWorld>;

  Dune::IGA::NURBSPatchData<2, dimWorld> patchData;
  patchData.knotSpans     = knotSpans;
  patchData.degree        = order;
  patchData.controlPoints = controlNet;
  using namespace Dune::IGA;
  patchData = degreeElevate(patchData, 1, 1);

  auto grid = std::make_shared<Grid>(patchData);

  grid->globalRefineInDirection(0, 2);
  grid->globalRefineInDirection(1, 2);
  auto gridView        = grid->leafGridView();
  const auto& indexSet = gridView.indexSet();

  using GridView = decltype(gridView);
  using namespace Ikarus;
  using namespace Dune::Functions::BasisFactory;
  const double E         = 4.32e4;
  const double nu        = 0.3;
  const double thickness = 0.00001;
  auto basis             = Ikarus::makeBasis(gridView, power<3>(nurbs()));

  /// Neumann boundary condition at right edge (x=0.5)
  Dune::BitSetVector<1> neumannVertices(gridView.size(2), false);
  for (auto&& vertex : vertices(gridView))
    if (std::abs(vertex.geometry().corner(0)[0] - 0.5) < 1e-8) neumannVertices[indexSet.index(vertex)] = true;

  BoundaryPatch<decltype(gridView)> neumannBoundary(gridView, neumannVertices);

  auto neumannBoundaryLoad = [&](auto& globalCoord, auto& lamb) {
    Eigen::Vector3d fext;
    fext.setZero();
    fext[0] = -lamb;
    return fext;
  };

  using ElementType = Ikarus::AutoDiffFE<Ikarus::KirchhoffLoveShell<decltype(basis)>>;
  std::vector<ElementType> fes;

  for (auto& element : elements(gridView))
    fes.emplace_back(basis, element, E, nu, thickness, LoadDefault(), &neumannBoundary, neumannBoundaryLoad);

  auto basisP = std::make_shared<const decltype(basis)>(basis);
  Ikarus::DirichletValues dirichletValues(basisP->flat());

  /// fix all DOFs at the left edge (x=-0.5)
  dirichletValues.fixDOFs([&](auto& basis_, auto&& dirichletFlags) {
    Dune::Functions::forEachBoundaryDOF(basis_, [&](auto&& localIndex, auto&& localView, auto&& intersection) {
      if (std::abs(intersection.geometry().center()[0] + 0.5) < 1e-8)
        dirichletFlags[localView.index(localIndex)] = true;
    });
  });

  /// fix y and z DOFs at the right edge (x=0.5)
  dirichletValues.fixDOFs([&](auto& basis_, auto&& dirichletFlags) {
    for (auto fixedDirection = 1; fixedDirection < 3; ++fixedDirection) {
      Dune::Functions::forEachBoundaryDOF(Dune::Functions::subspaceBasis(basis_, fixedDirection),
                                          [&](auto&& localIndex, auto&& localView, auto&& intersection) {
                                            if (std::abs(intersection.geometry().center()[0] - 0.5) < 1e-8)
                                              dirichletFlags[localView.index(localIndex)] = true;
                                          });
    }
  });

  auto sparseAssembler = SparseFlatAssembler(fes, dirichletValues);

  Eigen::VectorXd d;
  d.setZero(basis.flat().size());
  double lambda = 0.0;

  auto req = FErequirements().addAffordance(Ikarus::AffordanceCollections::elastoStatics);

  auto residualFunction = [&](auto&& disp_, auto&& lambdaLocal) -> auto& {
    req.insertGlobalSolution(Ikarus::FESolutions::displacement, disp_)
        .insertParameter(Ikarus::FEParameter::loadfactor, lambdaLocal);
    return sparseAssembler.getVector(req);
  };

  auto KFunction = [&](auto&& disp_, auto&& lambdaLocal) -> auto& {
    req.insertGlobalSolution(Ikarus::FESolutions::displacement, disp_)
        .insertParameter(Ikarus::FEParameter::loadfactor, lambdaLocal);
    return sparseAssembler.getMatrix(req);
  };

  auto nonLinOp  = Ikarus::NonLinearOperator(functions(residualFunction, KFunction), parameter(d, lambda));
  auto linSolver = Ikarus::LinearSolver(Ikarus::SolverTypeTag::sd_SimplicialLDLT);

  int loadSteps = 6;

  auto nr   = Ikarus::makeNewtonRaphsonWithSubsidiaryFunction(nonLinOp, std::move(linSolver));
  auto dass = Ikarus::AdaptiveStepSizing::IterationBased{};
  auto nass = Ikarus::AdaptiveStepSizing::NoOp{};
  dass.setTargetIterations(targetIterations);

  /// control routine with and without step sizing
  auto crWSS  = Ikarus::PathFollowing(nr, loadSteps, stepSize, pft, dass);
  auto crWoSS = Ikarus::PathFollowing(nr, loadSteps, stepSize, pft, nass);

  auto nonLinearSolverObserver = std::make_shared<Ikarus::NonLinearSolverLogger>();
  auto pathFollowingObserver   = std::make_shared<Ikarus::ControlLogger>();
  nr->subscribeAll(nonLinearSolverObserver);

  t.checkThrow<Dune::InvalidStateException>(
      [&]() { nonLinearSolverObserver->update(Ikarus::NonLinearSolverMessages::BEGIN); },
      "nonLinearSolverObserver should have failed for the BEGIN message");

  t.checkThrow<Dune::InvalidStateException>(
      [&]() { nonLinearSolverObserver->update(Ikarus::NonLinearSolverMessages::END); },
      "nonLinearSolverObserver should have failed for the END message");

  /// Create Observer which writes vtk files when control routines messages
  /// SOLUTION_CHANGED
  auto vtkWriter = std::make_shared<ControlSubsamplingVertexVTKWriter<std::remove_cvref_t<decltype(basis.flat())>>>(
      basis.flat(), d, 2);
  vtkWriter->setFieldInfo("displacement", Dune::VTK::FieldInfo::Type::vector, 3);
  vtkWriter->setFileNamePrefix("testAdaptiveStepSizing" + pft.name);

  crWoSS.subscribeAll({vtkWriter, pathFollowingObserver});
  crWoSS.unSubscribeAll(vtkWriter);
  crWSS.subscribeAll({pathFollowingObserver});
  crWSS.subscribe(Ikarus::ControlMessages::SOLUTION_CHANGED, vtkWriter);

  const std::string& message1 = " --> " + pft.name + " with default adaptive step sizing";
  const std::string& message2 = " --> " + pft.name + " without default adaptive step sizing";

  t.checkThrow<Dune::InvalidStateException>(
      [&]() {
        auto dass2 = Ikarus::AdaptiveStepSizing::IterationBased{};
        dass2.setTargetIterations(0);
        auto crWSS2                = Ikarus::PathFollowing(nr, loadSteps, stepSize, pft, dass2);
        const auto controlInfoWSS2 = crWSS2.run();
      },
      "IterationBased should fail for targetIterations being 0");

  resetNonLinearOperatorParametersToZero(nonLinOp);
  const auto controlInfoWSS = crWSS.run();
  checkScalars(t, std::ranges::max(d), expectedResults[0][0], message1 + "<Max Displacement>");
  checkScalars(t, lambda, expectedResults[0][1], message1 + "<Lambda>");
  resetNonLinearOperatorParametersToZero(nonLinOp);

  const auto controlInfoWoSS = crWoSS.run();
  checkScalars(t, std::ranges::max(d), expectedResults[1][0], message2 + "<Max Displacement>");
  checkScalars(t, lambda, expectedResults[1][1], message2 + "<Lambda>");
  resetNonLinearOperatorParametersToZero(nonLinOp);

  const int controlInfoWSSIterations
      = std::accumulate(controlInfoWSS.solverInfos.begin(), controlInfoWSS.solverInfos.end(), 0,
                        [](int a, auto& b) { return b.iterations + a; });

  t.check(controlInfoWSS.success, "No convergence" + message1);
  t.check(controlInfoWoSS.success, "No convergence" + message2);

  t.check(controlInfoWSS.totalIterations < controlInfoWoSS.totalIterations)
      << "Total iterations should be less --> " << controlInfoWSS.totalIterations << " > "
      << std::to_string(controlInfoWoSS.totalIterations) + " --> " + pft.name;

  t.check(controlInfoWSSIterations == controlInfoWSS.totalIterations)
      << "Total number of iterations is wrong --> " << controlInfoWSS.totalIterations << " is not equal to "
      << std::to_string(controlInfoWSSIterations) + " --> " + pft.name;

  checkSolverInfos(t, expectedIterations[0], controlInfoWSS, loadSteps, message1);
  checkSolverInfos(t, expectedIterations[1], controlInfoWoSS, loadSteps, message2);

  return t;
}

int main(int argc, char** argv) {
  Ikarus::init(argc, argv);

  auto alc = Ikarus::StandardArcLength{};
  auto lc  = Ikarus::LoadControlWithSubsidiaryFunction{};

  /// expected iterations for each step for a path following type with and without step sizing
  const std::vector<std::vector<int>> expectedIterationsALC = {{9, 9, 6, 5, 4, 4}, {9, 8, 6, 5, 5, 5}};
  const std::vector<std::vector<int>> expectedIterationsLC  = {{10, 4, 4, 3, 3, 3}, {10, 6, 5, 4, 4, 4}};

  /// expected results(max(displacement) and lambda) for a path following type with and without step sizing
  const std::vector<std::vector<double>> expectedResultsALC
      = {{0.1032139637288574, 0.0003103004514250302}, {0.162759603260405, 0.0007765975850229621}};
  const std::vector<std::vector<double>> expectedResultsLC
      = {{0.08741028329554587, 0.0002318693543601816}, {0.144353999993308, 6e-4}};

  KLShellAndAdaptiveStepSizing(alc, expectedIterationsALC, expectedResultsALC, 3, 0.4);
  KLShellAndAdaptiveStepSizing(lc, expectedIterationsLC, expectedResultsLC, 2, 1e-4);
}
