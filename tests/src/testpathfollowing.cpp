// SPDX-FileCopyrightText: 2021-2022 The Ikarus Developers mueller@ibb.uni-stuttgart.de
// SPDX-License-Identifier: LGPL-3.0-or-later

//
#include <config.h>

#include "testhelpers.hh"

#include <dune/common/test/testsuite.hh>

#include <Eigen/Core>
#include <Eigen/Dense>

#include <ikarus/assembler/simpleassemblers.hh>
#include <ikarus/controlroutines/pathfollowingtechnique.hh>
#include <ikarus/linearalgebra/nonlinearoperator.hh>
#include <ikarus/utils/drawing/griddrawer.hh>
#include <ikarus/utils/functionsanitychecks.hh>
#include <ikarus/utils/init.hh>
#include <ikarus/utils/observer/nonlinearsolverlogger.hh>

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

static auto simple2DOperatorArcLengthTest() {
  double lambda = 0;
  Eigen::VectorXd D;
  D.setZero(2);

  auto fvLambda  = [&](auto&& D_, auto&& lambda_) { return residual(D_, lambda_); };
  auto dfvLambda = [&](auto&& D_, auto&& lambda_) { return stiffnessMatrix(D_, lambda_); };

  auto nonLinOp = Ikarus::NonLinearOperator(Ikarus::functions(fvLambda, dfvLambda), Ikarus::parameter(D, lambda));

  auto linSolver = Ikarus::LinearSolver(Ikarus::SolverTypeTag::d_LDLT);

  double stepSize = 0.1;
  int load_steps  = 50;

  auto pft = Ikarus::StandardArcLength{};  // Path following type

  static_assert(PathFollowingStrategy<decltype(pft), decltype(nonLinOp)>,
                "StandardArcLength is a PathFollowingStrategy");

  auto nr                      = Ikarus::makeNewtonRaphsonWithSubsidiaryFunction(nonLinOp, std::move(linSolver));
  auto alc                     = Ikarus::PathFollowing(nr, load_steps, stepSize, pft);
  auto nonLinearSolverObserver = std::make_shared<Ikarus::NonLinearSolverLogger>();
  nr->subscribeAll(nonLinearSolverObserver);
  const auto controlInfo = alc.run();

  TestSuite t("Arc Length with Subsidiary function");
  t.check(controlInfo.success, "Successful result");
  return t;
}

static auto simple2DOperatorArcLengthTestAsDefault() {
  double lambda = 0;
  Eigen::VectorXd D;
  D.setZero(2);

  auto fvLambda  = [&](auto&& D_, auto&& lambda_) { return residual(D_, lambda_); };
  auto dfvLambda = [&](auto&& D_, auto&& lambda_) { return stiffnessMatrix(D_, lambda_); };

  auto nonLinOp = Ikarus::NonLinearOperator(Ikarus::functions(fvLambda, dfvLambda), Ikarus::parameter(D, lambda));

  auto linSolver = Ikarus::LinearSolver(Ikarus::SolverTypeTag::d_LDLT);

  double stepSize = 0.1;
  int load_steps  = 50;

  auto nr                      = Ikarus::makeNewtonRaphsonWithSubsidiaryFunction(nonLinOp, std::move(linSolver));
  auto alc                     = Ikarus::PathFollowing(nr, load_steps, stepSize);
  auto nonLinearSolverObserver = std::make_shared<Ikarus::NonLinearSolverLogger>();
  nr->subscribeAll(nonLinearSolverObserver);
  const auto controlInfo = alc.run();

  TestSuite t("Arc Length as Default Test");
  t.check(controlInfo.success, "Successful result");
  return t;
}

static auto simple2DOperatorLoadControlTest() {
  double lambda = 0;
  Eigen::VectorXd D;
  D.setZero(2);

  auto fvLambda  = [&](auto&& D_, auto&& lambda_) { return residual(D_, lambda_); };
  auto dfvLambda = [&](auto&& D_, auto&& lambda_) { return stiffnessMatrix(D_, lambda_); };

  auto nonLinOp = Ikarus::NonLinearOperator(Ikarus::functions(fvLambda, dfvLambda), Ikarus::parameter(D, lambda));

  auto linSolver = Ikarus::LinearSolver(Ikarus::SolverTypeTag::d_LDLT);

  double stepSize = 0.1;
  int load_steps  = 50;

  auto pft = Ikarus::LoadControlWithSubsidiaryFunction{};  // Path following type

  static_assert(PathFollowingStrategy<decltype(pft), decltype(nonLinOp)>, "LoadControl is a PathFollowingStrategy");

  auto nr                      = Ikarus::makeNewtonRaphsonWithSubsidiaryFunction(nonLinOp, std::move(linSolver));
  auto lc                      = Ikarus::PathFollowing(nr, load_steps, stepSize, pft);
  auto nonLinearSolverObserver = std::make_shared<Ikarus::NonLinearSolverLogger>();
  nr->subscribeAll(nonLinearSolverObserver);
  const auto controlInfo = lc.run();

  TestSuite t("Load Control with Subsidiary function");
  t.check(controlInfo.success, "Successful result");
  return t;
}

static auto simple2DOperatorDisplacementControlTest() {
  double lambda = 0;
  Eigen::VectorXd D;
  D.setZero(2);

  auto fvLambda  = [&](auto&& D_, auto&& lambda_) { return residual(D_, lambda_); };
  auto dfvLambda = [&](auto&& D_, auto&& lambda_) { return stiffnessMatrix(D_, lambda_); };

  auto nonLinOp = Ikarus::NonLinearOperator(Ikarus::functions(fvLambda, dfvLambda), Ikarus::parameter(D, lambda));

  auto linSolver = Ikarus::LinearSolver(Ikarus::SolverTypeTag::d_LDLT);

  double stepSize                    = 0.05;
  int load_steps                     = 30;
  std::vector<int> controlledIndices = {0};

  auto pft = Ikarus::DisplacementControl{controlledIndices};  // Path following type

  static_assert(PathFollowingStrategy<decltype(pft), decltype(nonLinOp)>,
                "DisplacementControl is a PathFollowingStrategy");

  auto nr                      = Ikarus::makeNewtonRaphsonWithSubsidiaryFunction(nonLinOp, std::move(linSolver));
  auto dc                      = Ikarus::PathFollowing(nr, load_steps, stepSize, pft);
  auto nonLinearSolverObserver = std::make_shared<Ikarus::NonLinearSolverLogger>();
  nr->subscribeAll(nonLinearSolverObserver);
  const auto controlInfo = dc.run();

  TestSuite t("Displacement Control with Subsidiary function");
  t.check(controlInfo.success, "Successful result");
  return t;
}

int main(int argc, char** argv) {
  Ikarus::init(argc, argv);
  TestSuite t;

  t.subTest(simple2DOperatorArcLengthTest());
  t.subTest(simple2DOperatorArcLengthTestAsDefault());
  t.subTest(simple2DOperatorLoadControlTest());
  t.subTest(simple2DOperatorDisplacementControlTest());

  return t.exit();
}
