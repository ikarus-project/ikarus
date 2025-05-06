// SPDX-FileCopyrightText: 2021-2025 The Ikarus Developers ikarus@ibb.uni-stuttgart.de
// SPDX-License-Identifier: LGPL-3.0-or-later

#include <config.h>

#include "testhelpers.hh"

#include <dune/common/test/testsuite.hh>

#include <Eigen/Core>

#include <ikarus/controlroutines/pathfollowing.hh>
#include <ikarus/utils/differentiablefunction.hh>
#include <ikarus/utils/init.hh>
#include <ikarus/utils/listener/controllogger.hh>
#include <ikarus/utils/listener/nonlinearsolverlogger.hh>

using namespace Ikarus::Concepts;
using Dune::TestSuite;

using DummyFERequirements = Ikarus::FERequirements<Ikarus::FESolutions::displacement, Ikarus::FEParameter::loadfactor>;

static auto residual(const Eigen::VectorXd& D, double lambda) {
  Eigen::VectorXd vec;
  vec.resize(2);
  auto& w = D[1];
  auto& u = D[0];
  vec << 16.0 * u - 3 * w - w * w * w - (2.0 / 3.0) * lambda, -3 * u + 4 * w - 3 * u * w * w - (9.0 / 4.0) * lambda;
  return vec;
}
static auto stiffnessMatrix(const Eigen::VectorXd& D, [[maybe_unused]] double lambda) {
  Eigen::MatrixXd mat;
  mat.setZero(2, 2);
  auto& w = D[1];
  auto& u = D[0];
  mat << 16.0, -3.0 * w * w - 3, -3 * w * w - 3, -6 * u * w + 4;
  return mat;
}

template <typename DifferentiableFunction>
static auto simple2DOperatorArcLengthTest(DifferentiableFunction& f, typename DifferentiableFunction::Domain& req,
                                          double stepSize, int loadSteps) {
  req.globalSolution().setZero();
  req.parameter() = 0.0;
  auto linSolver  = Ikarus::LinearSolver(Ikarus::SolverTypeTag::d_LDLT);
  auto pft        = Ikarus::ArcLength{}; // Type of path following technique

  auto nrSettings              = Ikarus::NewtonRaphsonWithSubsidiaryFunctionConfig({}, linSolver);
  auto nr                      = Ikarus::createNonlinearSolver(nrSettings, f);
  auto alc                     = Ikarus::PathFollowing(nr, loadSteps, stepSize, pft);
  auto nonLinearSolverObserver = Ikarus::NonLinearSolverLogger().subscribeTo(nr);
  auto pathFollowingObserver   = Ikarus::ControlLogger().subscribeTo(alc);

  const auto controlInfo              = alc.run(req);
  std::vector<int> expectedIterations = {1, 3, 3, 3, 3};
  Eigen::Vector2d expectedDisplacement;
  expectedDisplacement << 0.0883524725970593, 0.3486891582376427;
  double expectedLambda = 0.4877655288280236;

  TestSuite t("Arc Length with Subsidiary function");
  t.check(controlInfo.success, "No convergence");
  for (auto i = 0; i < 2; ++i)
    checkScalars(t, req.globalSolution()[i], expectedDisplacement[i], " --> " + pft.name());
  checkScalars(t, req.parameter(), expectedLambda, " --> " + pft.name());
  checkSolverInfos(t, expectedIterations, controlInfo, loadSteps);
  return t;
}

template <typename DifferentiableFunction>
static auto simple2DOperatorArcLengthTestAsDefault(DifferentiableFunction& f,
                                                   typename DifferentiableFunction::Domain& req, double stepSize,
                                                   int loadSteps) {
  req.globalSolution().setZero();
  req.parameter() = 0.0;
  auto linSolver  = Ikarus::LinearSolver(Ikarus::SolverTypeTag::d_LDLT);
  auto nrSettings = Ikarus::NewtonRaphsonWithSubsidiaryFunctionConfig<decltype(linSolver)>{.linearSolver = linSolver};
  auto nr         = Ikarus::createNonlinearSolver(nrSettings, f);
  auto alc        = Ikarus::PathFollowing(nr, loadSteps, stepSize);
  auto nonLinearSolverObserver = Ikarus::NonLinearSolverLogger().subscribeTo(nr);
  auto pathFollowingObserver   = Ikarus::ControlLogger().subscribeTo(alc);

  const auto controlInfo              = alc.run(req);
  std::vector<int> expectedIterations = {1, 3, 3, 3, 3};
  Eigen::Vector2d expectedDisplacement;
  expectedDisplacement << 0.0883524725970593, 0.3486891582376427;
  double expectedLambda = 0.4877655288280236;

  TestSuite t("Arc Length as Default Test");
  t.check(controlInfo.success, "No convergence");
  for (auto i = 0; i < 2; ++i)
    checkScalars(t, req.globalSolution()[i], expectedDisplacement[i]);
  checkScalars(t, req.parameter(), expectedLambda);
  checkSolverInfos(t, expectedIterations, controlInfo, loadSteps);
  return t;
}

template <typename DifferentiableFunction>
static auto simple2DOperatorLoadControlTest(DifferentiableFunction& f, typename DifferentiableFunction::Domain& req,
                                            double stepSize, int loadSteps) {
  req.globalSolution().setZero();
  req.parameter() = 0.0;
  auto linSolver  = Ikarus::LinearSolver(Ikarus::SolverTypeTag::d_LDLT);
  auto pft        = Ikarus::LoadControlSubsidiaryFunction{}; // Type of path following technique

  auto nrSettings              = Ikarus::NewtonRaphsonWithSubsidiaryFunctionConfig({}, linSolver);
  auto nr                      = Ikarus::createNonlinearSolver(nrSettings, f);
  auto lc                      = Ikarus::PathFollowing(nr, loadSteps, stepSize, pft);
  auto nonLinearSolverObserver = Ikarus::NonLinearSolverLogger().subscribeTo(nr);
  auto pathFollowingObserver   = Ikarus::ControlLogger().subscribeTo(lc);

  const auto controlInfo              = lc.run(req);
  std::vector<int> expectedIterations = {2, 3, 3, 3, 3};
  Eigen::Vector2d expectedDisplacement;
  expectedDisplacement << 0.0908533884835060, 0.3581294588381901;
  double expectedLambda = 0.5;

  TestSuite t("Load Control with Subsidiary function");
  t.check(controlInfo.success, "No convergence");
  for (auto i = 0; i < 2; ++i)
    checkScalars(t, req.globalSolution()[i], expectedDisplacement[i], " --> " + pft.name());
  checkScalars(t, req.parameter(), expectedLambda, " --> " + pft.name());
  checkSolverInfos(t, expectedIterations, controlInfo, loadSteps);
  return t;
}

template <typename DifferentiableFunction>
static auto simple2DOperatorDisplacementControlTest(DifferentiableFunction& f,
                                                    typename DifferentiableFunction::Domain& req, double stepSize,
                                                    int loadSteps) {
  req.globalSolution().setZero();
  req.parameter()                    = 0.0;
  auto linSolver                     = Ikarus::LinearSolver(Ikarus::SolverTypeTag::d_LDLT);
  std::vector<int> controlledIndices = {0};

  auto pft = Ikarus::DisplacementControl{controlledIndices}; // Type of path following technique

  auto nrSettings              = Ikarus::NewtonRaphsonWithSubsidiaryFunctionConfig({}, linSolver);
  auto nr                      = Ikarus::createNonlinearSolver(nrSettings, f);
  auto dc                      = Ikarus::PathFollowing(nr, loadSteps, stepSize, pft);
  auto nonLinearSolverObserver = Ikarus::NonLinearSolverLogger().subscribeTo(nr);
  auto pathFollowingObserver   = Ikarus::ControlLogger().subscribeTo(dc);

  const auto controlInfo              = dc.run(req);
  std::vector<int> expectedIterations = {3, 3, 3, 3, 3};
  Eigen::Vector2d expectedDisplacement;
  expectedDisplacement << 0.5, 1.4781013410920430;
  double expectedLambda = 0.5045466678049050;

  TestSuite t("Displacement Control with Subsidiary function");
  t.check(controlInfo.success, "No convergence");
  for (auto i = 0; i < 2; ++i)
    checkScalars(t, req.globalSolution()[i], expectedDisplacement[i], " --> " + pft.name());
  checkScalars(t, req.parameter(), expectedLambda, " --> " + pft.name());
  checkSolverInfos(t, expectedIterations, controlInfo, loadSteps);
  return t;
}

int main(int argc, char** argv) {
  Ikarus::init(argc, argv);
  TestSuite t;

  double lambda = 0;
  Eigen::VectorXd D;
  D.setZero(2);

  auto fvLambda  = [&](auto&& req_) { return residual(req_.globalSolution(), req_.parameter()); };
  auto dfvLambda = [&](auto&& req_) { return stiffnessMatrix(req_.globalSolution(), req_.parameter()); };
  DummyFERequirements req;
  req.insertGlobalSolution(D);
  req.insertParameter(lambda);

  auto f = Ikarus::makeDifferentiableFunction(Ikarus::functions(fvLambda, dfvLambda), req);

  double stepSize = 0.1;
  int loadSteps   = 5;

  t.subTest(simple2DOperatorArcLengthTest(f, req, stepSize, loadSteps));
  t.subTest(simple2DOperatorArcLengthTestAsDefault(f, req, stepSize, loadSteps));
  t.subTest(simple2DOperatorLoadControlTest(f, req, stepSize, loadSteps));
  t.subTest(simple2DOperatorDisplacementControlTest(f, req, stepSize, loadSteps));

  return t.exit();
}
