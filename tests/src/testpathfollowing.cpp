// SPDX-FileCopyrightText: 2021-2024 The Ikarus Developers mueller@ibb.uni-stuttgart.de
// SPDX-License-Identifier: LGPL-3.0-or-later

#include <config.h>

#include "testcommon.hh"
#include "testhelpers.hh"

#include <dune/common/test/testsuite.hh>

#include <Eigen/Core>

#include <ikarus/controlroutines/loadcontrol.hh>
#include <ikarus/controlroutines/pathfollowing.hh>
#include <ikarus/solver/nonlinearsolver/newtonraphson.hh>
#include <ikarus/utils/init.hh>
#include <ikarus/utils/nonlinearoperator.hh>
#include <ikarus/utils/observer/controllogger.hh>
#include <ikarus/utils/observer/genericobserver.hh>
#include <ikarus/utils/observer/nonlinearsolverlogger.hh>
#include <ikarus/utils/observer/observable.hh>

using namespace Ikarus::Concepts;
using Dune::TestSuite;

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

class OurDispObserver : public Ikarus::IObserver<Ikarus::ControlObservable>
{
public:
  void updateImpl(MessageType message, const StateType& state) override {
    if (message == Ikarus::ControlMessages::CONTROL_STARTED or message == Ikarus::ControlMessages::SOLUTION_CHANGED)
      std::cout << std::setprecision(16) << "The displacement is \t" << state.sol->transpose() << std::endl;
  }
};

template <typename NonLinearOperator>
static auto simple2DOperatorArcLengthTest(NonLinearOperator& nonLinOp, double stepSize, size_t loadSteps) {
  resetNonLinearOperatorParametersToZero(nonLinOp);
  auto linSolver = Ikarus::LinearSolver(Ikarus::SolverTypeTag::d_LDLT);
  auto pft       = Ikarus::ArcLength{}; // Type of path following technique

  auto nrSettings = Ikarus::NewtonRaphsonWithSubsidiaryFunctionConfig<decltype(linSolver)>{.linearSolver = linSolver};
  auto nr         = Ikarus::createNonlinearSolver(nrSettings, nonLinOp);

  auto alc = Ikarus::PathFollowing(nr, loadSteps, stepSize, pft);

  auto nonLinearSolverObserver = std::make_shared<Ikarus::NonLinearSolverLogger>();
  auto pathFollowingObserver   = std::make_shared<Ikarus::ControlLogger>();
  nr->subscribeAll(nonLinearSolverObserver);
  alc.subscribeAll(pathFollowingObserver);
  const auto controlState             = alc.run();
  std::vector<int> expectedIterations = {1, 3, 3, 3, 3};
  Eigen::Vector2d expectedDisplacement;
  expectedDisplacement << 0.0883524725970593, 0.3486891582376427;
  double expectedLambda = 0.4877655288280236;

  TestSuite t("Arc Length with Subsidiary function");
  t.check(controlState.success, "No convergence");
  for (auto i = 0; i < 2; ++i)
    checkScalars(t, nonLinOp.firstParameter()[i], expectedDisplacement[i], " --> " + alc.name());
  checkScalars(t, nonLinOp.lastParameter(), expectedLambda, " --> " + alc.name());

  checkSolverInfos(t, expectedIterations, controlState, loadSteps);
  return t;
}

template <typename NonLinearOperator>
static auto simple2DOperatorArcLengthTestAsDefault(NonLinearOperator& nonLinOp, double stepSize, size_t loadSteps) {
  resetNonLinearOperatorParametersToZero(nonLinOp);
  auto linSolver  = Ikarus::LinearSolver(Ikarus::SolverTypeTag::d_LDLT);
  auto nrSettings = Ikarus::NewtonRaphsonWithSubsidiaryFunctionConfig<decltype(linSolver)>{.linearSolver = linSolver};
  auto nr         = Ikarus::createNonlinearSolver(nrSettings, nonLinOp);
  auto alc        = Ikarus::PathFollowing(nr, loadSteps, stepSize);
  auto nonLinearSolverObserver = std::make_shared<Ikarus::NonLinearSolverLogger>();
  auto pathFollowingObserver   = std::make_shared<Ikarus::ControlLogger>();
  nr->subscribeAll(nonLinearSolverObserver);
  alc.subscribeAll(pathFollowingObserver);
  const auto controlState             = alc.run();
  std::vector<int> expectedIterations = {1, 3, 3, 3, 3};
  Eigen::Vector2d expectedDisplacement;
  expectedDisplacement << 0.0883524725970593, 0.3486891582376427;
  double expectedLambda = 0.4877655288280236;

  TestSuite t("Arc Length as Default Test");
  t.check(controlState.success, "No convergence");
  for (auto i = 0; i < 2; ++i)
    checkScalars(t, nonLinOp.firstParameter()[i], expectedDisplacement[i]);
  checkScalars(t, nonLinOp.lastParameter(), expectedLambda);

  checkSolverInfos(t, expectedIterations, controlState, loadSteps);
  return t;
}

template <typename NonLinearOperator>
static auto simple2DOperatorLoadControlAsSubsidiaryFunctionTest(NonLinearOperator& nonLinOp, double stepSize,
                                                                size_t loadSteps) {
  resetNonLinearOperatorParametersToZero(nonLinOp);
  TestSuite t("Load Control as Subsidiary function");
  std::vector<int> expectedIterations = {2, 3, 3, 3, 3};
  Eigen::Matrix2Xd expectedDisplacement;
  expectedDisplacement.setZero(Eigen::NoChange, loadSteps);
  expectedDisplacement << 0.01715872957844366, 0.0345464428730192, 0.0524126112865617, 0.0710534689402604,
      0.0908533884835060, 0.0691806374841585, 0.1389097864303651, 0.2097895325120464, 0.2825443193976919,
      0.3581294588381901;
  double expectedLambda = 0.5;
  auto linSolver        = Ikarus::LinearSolver(Ikarus::SolverTypeTag::d_LDLT);
  auto pft              = Ikarus::LoadControlSubsidiaryFunction{}; // Type of path following technique

  auto nrSettings = Ikarus::NewtonRaphsonWithSubsidiaryFunctionConfig<decltype(linSolver)>{.linearSolver = linSolver};
  auto nr         = Ikarus::createNonlinearSolver(nrSettings, nonLinOp);
  auto lc         = Ikarus::PathFollowing(nr, loadSteps, stepSize, pft);
  auto nonLinearSolverObserver = std::make_shared<Ikarus::NonLinearSolverLogger>();
  auto pathFollowingObserver   = std::make_shared<Ikarus::ControlLogger>();

  // create an observer to print only the solution vector
  auto ourDispObserver = std::make_shared<OurDispObserver>();

  /// Create GenericObserver which executes when control routines messages
  Eigen::Matrix2Xd dispMat;
  dispMat.setZero(Eigen::NoChange, loadSteps);
  auto dispObserver = std::make_shared<Ikarus::GenericObserver<Ikarus::ControlObservable>>(
      Ikarus::ControlMessages::SOLUTION_CHANGED, [&](int step) {
        const auto& d    = nonLinOp.firstParameter();
        dispMat(0, step) = d[0];
        dispMat(1, step) = d[1];
      });

  nr->subscribeAll(nonLinearSolverObserver);
  lc.subscribeAll({pathFollowingObserver, dispObserver, ourDispObserver});
  const auto controlState = lc.run();
  t.check(controlState.success, "No convergence");
  for (auto i = 0; i < 2; ++i)
    for (auto j = 0; j < loadSteps; ++j)
      checkScalars(t, dispMat(i, j), expectedDisplacement(i, j), " --> " + lc.name());
  checkScalars(t, nonLinOp.lastParameter(), expectedLambda, " --> " + lc.name());

  checkSolverInfos(t, expectedIterations, controlState, loadSteps);
  return t;
}

template <typename NonLinearOperator>
static auto simple2DOperatorLoadControlTest(NonLinearOperator& nonLinOp, double stepSize, size_t loadSteps) {
  resetNonLinearOperatorParametersToZero(nonLinOp);
  TestSuite t("Load Control Method");
  std::vector<int> expectedIterations = {2, 3, 3, 3, 3};
  Eigen::Matrix2Xd expectedDisplacement;
  expectedDisplacement.setZero(Eigen::NoChange, loadSteps);
  expectedDisplacement << 0.01715872957844366, 0.0345464428730192, 0.0524126112865617, 0.0710534689402604,
      0.0908533884835060, 0.0691806374841585, 0.1389097864303651, 0.2097895325120464, 0.2825443193976919,
      0.3581294588381901;
  double expectedLambda = 0.5;

  auto linSolver = Ikarus::LinearSolver(Ikarus::SolverTypeTag::d_LDLT);
  auto nr        = Ikarus::makeNewtonRaphson(nonLinOp, std::move(linSolver));
  auto lc        = Ikarus::LoadControl(nr, loadSteps, {0, 0.5});

  auto nonLinearSolverObserver = std::make_shared<Ikarus::NonLinearSolverLogger>();
  auto pathFollowingObserver   = std::make_shared<Ikarus::ControlLogger>();

  // create an observer to print only the solution vector
  auto ourDispObserver = std::make_shared<OurDispObserver>();

  /// Create GenericObserver which executes when control routines messages
  Eigen::Matrix2Xd dispMat;
  dispMat.setZero(Eigen::NoChange, loadSteps);
  auto dispObserver = std::make_shared<Ikarus::GenericObserver<Ikarus::ControlObservable>>(
      Ikarus::ControlMessages::SOLUTION_CHANGED, [&](int step) {
        const auto& d    = nonLinOp.firstParameter();
        dispMat(0, step) = d[0];
        dispMat(1, step) = d[1];
      });

  nr->subscribeAll(nonLinearSolverObserver);
  lc.subscribeAll({pathFollowingObserver, dispObserver, ourDispObserver});
  const auto controlState = lc.run();
  t.check(controlState.success, "No convergence");
  for (auto i = 0; i < 2; ++i)
    for (auto j = 0; j < loadSteps; ++j)
      checkScalars(t, dispMat(i, j), expectedDisplacement(i, j), " --> " + lc.name());
  checkScalars(t, nonLinOp.lastParameter(), expectedLambda, " --> " + lc.name());

  checkSolverInfos(t, expectedIterations, controlState, loadSteps);

  t.checkThrow<Dune::InvalidStateException>(
      [&]() { auto lc1 = Ikarus::LoadControl(nr, 0, {0, 1}); },
      "An object of LoadControl should not have been constructed with loadSteps = 0.");
  return t;
}

template <typename NonLinearOperator>
static auto simple2DOperatorDisplacementControlTest(NonLinearOperator& nonLinOp, double stepSize, size_t loadSteps) {
  resetNonLinearOperatorParametersToZero(nonLinOp);
  auto linSolver                     = Ikarus::LinearSolver(Ikarus::SolverTypeTag::d_LDLT);
  std::vector<int> controlledIndices = {0};

  auto pft = Ikarus::DisplacementControl{controlledIndices}; // Type of path following technique

  auto nrSettings = Ikarus::NewtonRaphsonWithSubsidiaryFunctionConfig<decltype(linSolver)>{.linearSolver = linSolver};
  auto nr         = Ikarus::createNonlinearSolver(nrSettings, nonLinOp);
  auto dc         = Ikarus::PathFollowing(nr, loadSteps, stepSize, pft);
  auto nonLinearSolverObserver = std::make_shared<Ikarus::NonLinearSolverLogger>();
  auto pathFollowingObserver   = std::make_shared<Ikarus::ControlLogger>();
  nr->subscribeAll(nonLinearSolverObserver);
  dc.subscribeAll(pathFollowingObserver);
  const auto controlState             = dc.run();
  std::vector<int> expectedIterations = {3, 3, 3, 3, 3};
  Eigen::Vector2d expectedDisplacement;
  expectedDisplacement << 0.5, 1.4781013410920430;
  double expectedLambda = 0.5045466678049050;

  TestSuite t("Displacement Control with Subsidiary function");
  t.check(controlState.success, "No convergence");
  for (auto i = 0; i < 2; ++i)
    checkScalars(t, nonLinOp.firstParameter()[i], expectedDisplacement[i], " --> " + dc.name());
  checkScalars(t, nonLinOp.lastParameter(), expectedLambda, " --> " + dc.name());

  checkSolverInfos(t, expectedIterations, controlState, loadSteps);
  return t;
}

int main(int argc, char** argv) {
  Ikarus::init(argc, argv);
  TestSuite t;

  double lambda = 0;
  Eigen::VectorXd D;
  D.setZero(2);

  auto fvLambda  = [&](auto&& D_, auto&& lambda_) { return residual(D_, lambda_); };
  auto dfvLambda = [&](auto&& D_, auto&& lambda_) { return stiffnessMatrix(D_, lambda_); };

  auto nonLinOp = Ikarus::NonLinearOperator(Ikarus::functions(fvLambda, dfvLambda), Ikarus::parameter(D, lambda));

  constexpr double stepSize  = 0.1;
  constexpr size_t loadSteps = 5;

  t.subTest(simple2DOperatorArcLengthTest(nonLinOp, stepSize, loadSteps));
  t.subTest(simple2DOperatorArcLengthTestAsDefault(nonLinOp, stepSize, loadSteps));
  t.subTest(simple2DOperatorLoadControlAsSubsidiaryFunctionTest(nonLinOp, stepSize, loadSteps));
  t.subTest(simple2DOperatorLoadControlTest(nonLinOp, stepSize, loadSteps));
  t.subTest(simple2DOperatorDisplacementControlTest(nonLinOp, stepSize, loadSteps));

  return t.exit();
}
