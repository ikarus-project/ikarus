// SPDX-FileCopyrightText: 2021-2025 The Ikarus Developers ikarus@ibb.uni-stuttgart.de
// SPDX-License-Identifier: LGPL-3.0-or-later

#include <config.h>

#include "testhelpers.hh"

#include <dune/common/test/testsuite.hh>

#include <Eigen/Core>

#include <ikarus/controlroutines/loadcontrol.hh>
#include <ikarus/controlroutines/pathfollowing.hh>
#include <ikarus/solver/nonlinearsolver/newtonraphson.hh>
#include <ikarus/solver/nonlinearsolver/nonlinearsolverfactory.hh>
#include <ikarus/utils/differentiablefunction.hh>
#include <ikarus/utils/init.hh>
#include <ikarus/utils/listener/controllogger.hh>
#include <ikarus/utils/listener/genericlistener.hh>
#include <ikarus/utils/listener/nonlinearsolverlogger.hh>

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

  auto nrSettings            = Ikarus::NewtonRaphsonWithSubsidiaryFunctionConfig({}, linSolver);
  auto nr                    = Ikarus::createNonlinearSolver(nrSettings, f);
  auto alc                   = Ikarus::PathFollowing(nr, loadSteps, stepSize, pft);
  auto nonLinearSolverLogger = Ikarus::NonLinearSolverLogger().subscribeTo(*nr);
  auto pathFollowingLogger   = Ikarus::ControlLogger().subscribeTo(alc);

  const auto controlInfo              = alc.run(req);
  std::vector<int> expectedIterations = {0, 2, 2, 2, 2, 2};
  Eigen::Vector2d expectedDisplacement;
  expectedDisplacement << 0.0883524681713636, 0.3486891415306083;
  double expectedLambda = 0.4877655071723067;

  TestSuite t("Arc Length with Subsidiary function");
  t.check(controlInfo.success, "No convergence");
  for (auto i = 0; i < 2; ++i)
    checkScalars(t, req.globalSolution()[i], expectedDisplacement[i], " --> " + alc.name());
  checkScalars(t, req.parameter(), expectedLambda, " --> " + alc.name());
  checkSolverInfos(t, expectedIterations, controlInfo, loadSteps + 1);
  return t;
}

template <typename DifferentiableFunction>
static auto simple2DOperatorArcLengthTestWithIDBC(DifferentiableFunction& f,
                                                  typename DifferentiableFunction::Domain& req, double stepSize,
                                                  int loadSteps) {
  req.globalSolution().setZero();
  req.parameter() = 0.0;
  auto linSolver  = Ikarus::LinearSolver(Ikarus::SolverTypeTag::d_LDLT);
  auto pft        = Ikarus::ArcLength{}; // Type of path following technique

  // We only solve for w dof and the u is an inhomoegenous boundary condition
  auto updateFunction = []<typename D, typename C>(D& x, const C& dx) {
    if constexpr (not std::is_same_v<C, Ikarus::utils::SyncFERequirements>)
      x[1] += dx[0]; // here x[1] = w
    else {
      if constexpr (requires { x.parameter(); })
        x.globalSolution()[0] = 0.1 * x.parameter(); // update u = lambda * 0.1
    }
  };

  auto idbcForceFunction = [&]<typename D>(const D& x) {
    const auto K = stiffnessMatrix(x.globalSolution(), x.parameter());
    return Eigen::Vector<double, 1>{(K.col(0) * 0.1)[1]};
  };

  auto nrSettings = Ikarus::NewtonRaphsonWithSubsidiaryFunctionConfig({}, linSolver, updateFunction, idbcForceFunction);
  auto nr         = Ikarus::createNonlinearSolver(nrSettings, f);
  auto alc        = Ikarus::PathFollowing(nr, loadSteps, stepSize, pft);
  auto nonLinearSolverObserver = Ikarus::NonLinearSolverLogger().subscribeTo(*nr);
  auto pathFollowingObserver   = Ikarus::ControlLogger().subscribeTo(alc);

  const auto controlInfo              = alc.run(req);
  std::vector<int> expectedIterations = {0, 1, 2, 2, 2, 2};
  Eigen::Vector2d expectedDisplacement;
  expectedDisplacement << 0.05440020617136589, 0.351852377380214;
  double expectedLambda = 0.5440020617136589;

  TestSuite t("Arc Length with Inhomogeneous Dirichlet BCs");
  t.check(controlInfo.success, "No convergence");
  for (auto i = 0; i < 2; ++i)
    checkScalars(t, req.globalSolution()[i], expectedDisplacement[i], " --> " + alc.name());
  checkScalars(t, req.parameter(), expectedLambda, " --> " + alc.name());
  checkSolverInfos(t, expectedIterations, controlInfo, loadSteps + 1);
  return t;
}

template <typename DifferentiableFunction>
static auto simple2DOperatorLoadControlLCWithIDBC(DifferentiableFunction& f,
                                                  typename DifferentiableFunction::Domain& req, double stepSize,
                                                  int loadSteps) {
  req.globalSolution().setZero();
  req.parameter() = 0.0;
  auto linSolver  = Ikarus::LinearSolver(Ikarus::SolverTypeTag::d_LDLT);

  // We only solve for w dof and the u is an inhomoegenous boundary condition
  auto updateFunction = []<typename D, typename C>(D& x, const C& dx) {
    if constexpr (not std::is_same_v<C, Ikarus::utils::SyncFERequirements>)
      x[1] += dx[0]; // here x[1] = w
    else {
      if constexpr (requires { x.parameter(); })
        x.globalSolution()[0] = 0.1 * x.parameter(); // update u = lambda * 0.1
    }
  };

  auto idbcForceFunction = [&]<typename D>(const D& x) {
    const auto K = stiffnessMatrix(x.globalSolution(), x.parameter());
    return Eigen::Vector<double, 1>{(K.col(0) * 0.1)[1]};
  };

  auto nrSettings              = Ikarus::NewtonRaphsonConfig({}, linSolver, updateFunction, idbcForceFunction);
  auto nr                      = Ikarus::createNonlinearSolver(nrSettings, f);
  auto lc                      = Ikarus::LoadControl(nr, loadSteps, {0.0, 0.5});
  auto nonLinearSolverObserver = Ikarus::NonLinearSolverLogger().subscribeTo(*nr);
  auto pathFollowingObserver   = Ikarus::ControlLogger().subscribeTo(lc);

  const auto controlInfo              = lc.run(req);
  std::vector<int> expectedIterations = {0, 2, 2, 2, 3, 3};
  Eigen::Vector2d expectedDisplacement;
  expectedDisplacement << 0.05, 0.322653959150839100;
  double expectedLambda = 0.5;

  TestSuite t("Load Control with Inhomogeneous Dirichlet BCs");
  t.check(controlInfo.success, "No convergence");
  for (auto i = 0; i < 2; ++i)
    checkScalars(t, req.globalSolution()[i], expectedDisplacement[i], " --> " + lc.name());
  checkScalars(t, req.parameter(), expectedLambda, " --> " + lc.name());
  checkSolverInfos(t, expectedIterations, controlInfo, loadSteps + 1);
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
  auto nonLinearSolverLogger = Ikarus::NonLinearSolverLogger().subscribeTo(*nr);
  auto pathFollowingLogger   = Ikarus::ControlLogger().subscribeTo(alc);

  const auto controlInfo              = alc.run(req);
  std::vector<int> expectedIterations = {0, 2, 2, 2, 2, 2};
  Eigen::Vector2d expectedDisplacement;
  expectedDisplacement << 0.0883524681713636, 0.3486891415306083;
  double expectedLambda = 0.4877655071723067;

  TestSuite t("Arc Length as Default Test");
  t.check(controlInfo.success, "No convergence");
  for (auto i = 0; i < 2; ++i)
    checkScalars(t, req.globalSolution()[i], expectedDisplacement[i]);
  checkScalars(t, req.parameter(), expectedLambda);
  checkSolverInfos(t, expectedIterations, controlInfo, loadSteps + 1);
  return t;
}

template <typename DifferentiableFunction>
static auto simple2DOperatorLoadControlTestPF(DifferentiableFunction& f, typename DifferentiableFunction::Domain& req,
                                              double stepSize, int loadSteps) {
  req.globalSolution().setZero();
  req.parameter() = 0.0;
  auto linSolver  = Ikarus::LinearSolver(Ikarus::SolverTypeTag::d_LDLT);
  auto pft        = Ikarus::LoadControlSubsidiaryFunction{}; // Type of path following technique

  auto nrSettings            = Ikarus::NewtonRaphsonWithSubsidiaryFunctionConfig({}, linSolver);
  auto nr                    = Ikarus::createNonlinearSolver(nrSettings, f);
  auto lc                    = Ikarus::PathFollowing(nr, loadSteps, stepSize, pft);
  auto nonLinearSolverLogger = Ikarus::NonLinearSolverLogger().subscribeTo(*nr);
  auto pathFollowingLogger   = Ikarus::ControlLogger().subscribeTo(lc);

  /// Create GenericListener which executes when control routines messages to check displacements at every step
  Eigen::Matrix2Xd dispMat;
  dispMat.setZero(Eigen::NoChange, loadSteps + 1);
  auto genericListenerFunctor = [&](const auto& state) {
    const auto& d        = state.domain.globalSolution();
    int step             = state.loadStep;
    dispMat(0, step + 1) = d[0];
    dispMat(1, step + 1) = d[1];
  };
  auto dispObserver = Ikarus::GenericListener(lc, Ikarus::ControlMessages::SOLUTION_CHANGED, genericListenerFunctor);

  const auto controlInfo              = lc.run(req);
  std::vector<int> expectedIterations = {0, 2, 3, 3, 3, 3};
  Eigen::Matrix2Xd expectedDisplacement;
  expectedDisplacement.setZero(Eigen::NoChange, loadSteps + 1);
  expectedDisplacement << 0.0, 0.01715872957844366, 0.0345464428730192, 0.0524126112865617, 0.0710534689402604,
      0.0908533884835060, 0.0, 0.0691806374841585, 0.1389097864303651, 0.2097895325120464, 0.2825443193976919,
      0.3581294588381901;
  double expectedLambda = 0.5;

  TestSuite t("Load Control with Subsidiary function");
  t.check(controlInfo.success, "No convergence");
  for (auto i = 0; i < 2; ++i)
    for (auto j = 0; j < loadSteps + 1; ++j)
      checkScalars(t, dispMat(i, j), expectedDisplacement(i, j), " --> " + lc.name());
  checkScalars(t, req.parameter(), expectedLambda, " --> " + lc.name());
  checkSolverInfos(t, expectedIterations, controlInfo, loadSteps + 1);
  t.checkThrow<Dune::InvalidStateException>(
      [&]() { auto lc1 = Ikarus::PathFollowing(nr, 0, stepSize, pft); },
      "An object of PathFollowing should not have been constructed with steps = 0.");
  return t;
}

template <typename DifferentiableFunction>
static auto simple2DOperatorLoadControlTestLC(DifferentiableFunction& f, typename DifferentiableFunction::Domain& req,
                                              double stepSize, int loadSteps) {
  req.globalSolution().setZero();
  req.parameter() = 0.0;
  TestSuite t("Load Control Method");
  std::vector<int> expectedIterations = {0, 2, 3, 3, 3, 3};
  Eigen::Matrix2Xd expectedDisplacement;
  expectedDisplacement.setZero(Eigen::NoChange, loadSteps + 1);
  expectedDisplacement << 0.0, 0.01715872957844366, 0.0345464428730192, 0.0524126112865617, 0.0710534689402604,
      0.0908533884835060, 0.0, 0.0691806374841585, 0.1389097864303651, 0.2097895325120464, 0.2825443193976919,
      0.3581294588381901;
  double expectedLambda = 0.5;

  auto linSolver             = Ikarus::LinearSolver(Ikarus::SolverTypeTag::d_LDLT);
  auto nrSettings            = Ikarus::NewtonRaphsonConfig({}, linSolver);
  auto nr                    = Ikarus::createNonlinearSolver(nrSettings, f);
  auto lc                    = Ikarus::LoadControl(nr, loadSteps, {0.0, 0.5});
  auto nonLinearSolverLogger = Ikarus::NonLinearSolverLogger().subscribeTo(*nr);

  auto controlLogger = Ikarus::ControlLogger().subscribeTo(lc);

  /// Create GenericListener which executes when control routines messages to check displacements at every step
  Eigen::Matrix2Xd dispMat;
  dispMat.setZero(Eigen::NoChange, loadSteps + 1);

  auto dispObserver = Ikarus::GenericListener(lc, Ikarus::ControlMessages::SOLUTION_CHANGED, [&](const auto& state) {
    const auto& d        = state.domain.globalSolution();
    int step             = state.loadStep;
    dispMat(0, step + 1) = d[0];
    dispMat(1, step + 1) = d[1];
  });

  const auto controlInfo = lc.run(req);

  t.check(controlInfo.success, "No convergence");
  for (auto i = 0; i < 2; ++i)
    for (auto j = 0; j < loadSteps; ++j)
      checkScalars(t, dispMat(i, j), expectedDisplacement(i, j), " --> " + lc.name());
  checkScalars(t, req.parameter(), expectedLambda, " --> " + lc.name());
  checkSolverInfos(t, expectedIterations, controlInfo, loadSteps + 1);
  t.checkThrow<Dune::InvalidStateException>(
      [&]() { auto lc1 = Ikarus::LoadControl(nr, 0, {0, 1}); },
      "An object of LoadControl should not have been constructed with loadSteps = 0.");
  return t;
}

template <typename DifferentiableFunction>
static auto simple2DOperatorLoadControlTestLCWithDifferentListenerOrder(DifferentiableFunction& f,
                                                                        typename DifferentiableFunction::Domain& req,
                                                                        double stepSize, int loadSteps) {
  req.globalSolution().setZero();
  req.parameter() = 0.0;
  TestSuite t("Load Control Method");
  std::vector<int> expectedIterations = {0, 2, 3, 3, 3, 3};
  Eigen::Matrix2Xd expectedDisplacement;
  expectedDisplacement.setZero(Eigen::NoChange, loadSteps + 1);
  expectedDisplacement << 0.0, 0.01715872957844366, 0.0345464428730192, 0.0524126112865617, 0.0710534689402604,
      0.0908533884835060, 0.0, 0.0691806374841585, 0.1389097864303651, 0.2097895325120464, 0.2825443193976919,
      0.3581294588381901;
  double expectedLambda = 0.5;

  auto linSolver  = Ikarus::LinearSolver(Ikarus::SolverTypeTag::d_LDLT);
  auto nrSettings = Ikarus::NewtonRaphsonConfig({}, linSolver);
  auto nr         = Ikarus::createNonlinearSolver(nrSettings, f);
  auto lc         = Ikarus::LoadControl(nr, loadSteps, {0.0, 0.5});

  /// Create GenericListener which executes when control routines messages to check displacements at every step
  Eigen::Matrix2Xd dispMat;
  dispMat.setZero(Eigen::NoChange, loadSteps + 1);

  auto dispObserver = Ikarus::GenericListener(lc, Ikarus::ControlMessages::SOLUTION_CHANGED, [&](const auto& state) {
    const auto& d        = state.domain.globalSolution();
    int step             = state.loadStep;
    dispMat(0, step + 1) = d[0];
    dispMat(1, step + 1) = d[1];
  });

  auto controlLogger         = Ikarus::ControlLogger().subscribeTo(lc);
  auto nonLinearSolverLogger = Ikarus::NonLinearSolverLogger().subscribeTo(*nr);

  const auto controlInfo = lc.run(req);

  t.check(controlInfo.success, "No convergence");
  for (auto i = 0; i < 2; ++i)
    for (auto j = 0; j < loadSteps; ++j)
      checkScalars(t, dispMat(i, j), expectedDisplacement(i, j), " --> " + lc.name());
  checkScalars(t, req.parameter(), expectedLambda, " --> " + lc.name());
  checkSolverInfos(t, expectedIterations, controlInfo, loadSteps + 1);
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

  auto nrSettings            = Ikarus::NewtonRaphsonWithSubsidiaryFunctionConfig({}, linSolver);
  auto nr                    = Ikarus::createNonlinearSolver(nrSettings, f);
  auto dc                    = Ikarus::PathFollowing(nr, loadSteps, stepSize, pft);
  auto nonLinearSolverLogger = Ikarus::NonLinearSolverLogger().subscribeTo(*nr);
  auto pathFollowingLogger   = Ikarus::ControlLogger().subscribeTo(dc);

  const auto controlInfo              = dc.run(req);
  std::vector<int> expectedIterations = {0, 3, 3, 3, 3, 3};
  Eigen::Vector2d expectedDisplacement;
  expectedDisplacement << 0.5, 1.4781013410920430;
  double expectedLambda = 0.5045466678049050;

  TestSuite t("Displacement Control with Subsidiary function");
  t.check(controlInfo.success, "No convergence");
  for (auto i = 0; i < 2; ++i)
    checkScalars(t, req.globalSolution()[i], expectedDisplacement[i], " --> " + dc.name());
  checkScalars(t, req.parameter(), expectedLambda, " --> " + dc.name());
  checkSolverInfos(t, expectedIterations, controlInfo, loadSteps + 1);
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

  constexpr double stepSize = 0.1;
  constexpr int loadSteps   = 5;

  t.subTest(simple2DOperatorArcLengthTest(f, req, stepSize, loadSteps));
  t.subTest(simple2DOperatorArcLengthTestAsDefault(f, req, stepSize, loadSteps));
  t.subTest(simple2DOperatorLoadControlTestPF(f, req, stepSize, loadSteps));
  t.subTest(simple2DOperatorLoadControlTestLC(f, req, stepSize, loadSteps));
  t.subTest(simple2DOperatorLoadControlTestLCWithDifferentListenerOrder(f, req, stepSize, loadSteps));
  t.subTest(simple2DOperatorDisplacementControlTest(f, req, stepSize, loadSteps));

  auto fvLambdaRed = [&](auto&& req_) {
    return residual(req_.globalSolution(), req_.parameter()).segment(1, 1).eval();
  };
  auto dfvLambdaRed = [&](auto&& req_) {
    return stiffnessMatrix(req_.globalSolution(), req_.parameter()).block(1, 1, 1, 1).eval();
  };

  auto fred = Ikarus::makeDifferentiableFunction(Ikarus::functions(fvLambdaRed, dfvLambdaRed), req);

  t.subTest(simple2DOperatorArcLengthTestWithIDBC(fred, req, stepSize, loadSteps));
  t.subTest(simple2DOperatorLoadControlLCWithIDBC(fred, req, stepSize, loadSteps));

  return t.exit();
}