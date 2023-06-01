// SPDX-FileCopyrightText: 2022 The Ikarus Developers mueller@ibb.uni-stuttgart.de
// SPDX-License-Identifier: LGPL-3.0-or-later

//
#include <config.h>

#include "testHelpers.hh"

#include <dune/common/test/testsuite.hh>

#include <Eigen/Core>

#include <ikarus/linearAlgebra/nonLinearOperator.hh>
#include <ikarus/solver/nonLinearSolver/newtonRaphson.hh>
#include <ikarus/utils/functionSanityChecks.hh>
#include <ikarus/utils/init.hh>
#include <ikarus/utils/observer/nonLinearSolverLogger.hh>

using namespace Ikarus;
using Dune::TestSuite;

template <typename SolutionType, typename SolutionTypeExpected, typename NewtonRhapson>
auto checkNewtonRhapson(NewtonRhapson& nr, SolutionType& x, double tolerance, int maxIter, int iterExpected,
                        const SolutionTypeExpected& xExpected, const auto& x_Predictor) {
  TestSuite t("checkNewtonRhapson");
  nr.setup({tolerance, maxIter});
  const auto solverInfo = nr.solve(x_Predictor);

  if constexpr (std::is_same_v<SolutionType, double>)
    t.check(Dune::FloatCmp::eq(xExpected, x));
  else
    t.check(isApproxSame(x, xExpected, 1e-15));

  t.check(true == solverInfo.success) << "NewtonRhapson wasn't successful.";
  t.check(tolerance >= solverInfo.residualnorm)
      << "Residual norm is not below tolerance " << tolerance << " Actual: " << solverInfo.residualnorm;
  t.check(iterExpected == solverInfo.iterations)
      << "The iteration count does not match the expected number. Expected: " << iterExpected
      << " Actual: " << solverInfo.iterations;
  return t;
}

auto f(double x) { return 0.5 * x * x + x - 2; }
auto df(double x) { return x + 1; }

auto simple1DOperatorNewtonRhapsonTest() {
  TestSuite t("simple1DOperatorNewtonRhapsonTest");

  double x = 13;

  auto fvLambda  = [&](auto&& x_) { return f(x_); };
  auto dfvLambda = [&](auto&& x_) { return df(x_); };
  Ikarus::NonLinearOperator nonLinOp(functions(fvLambda, dfvLambda), parameter(x));

  // Newton method test
  const double eps       = 1e-14;
  const int maxIter      = 20;
  const double xExpected = std::sqrt(5.0) - 1.0;

  Ikarus::NewtonRaphson nr(nonLinOp);
  t.subTest(checkNewtonRhapson(nr, x, eps, maxIter, 7, xExpected, 0.0));
  return t;
}

auto simple1DOperatorNewtonRhapsonCheckThatThePerfectPredictorWorksTest() {
  TestSuite t("simple1DOperatorNewtonRhapsonCheckThatThePerfectPredictorWorksTest");
  double x = 0;

  auto fvLambda  = [](auto&& x_) { return f(x_); };
  auto dfvLambda = [](auto&& x_) { return df(x_); };
  Ikarus::NonLinearOperator nonLinOp(functions(fvLambda, dfvLambda), parameter(x));

  const double eps       = 1e-14;
  const int maxIter      = 20;
  const double xExpected = std::sqrt(5.0) - 1.0;

  Ikarus::NewtonRaphson nr(nonLinOp);

  t.subTest(checkNewtonRhapson(nr, x, eps, maxIter, 0, xExpected, xExpected));
  return t;
}

auto dfFail(double x) { return x + 1000000; }

auto simple1DOperatorNewtonRhapsonWithWrongDerivativeTest() {
  double x = 13;

  auto fvLambda  = [](auto&& x_) { return f(x_); };
  auto dfvLambda = [](auto&& x_) { return dfFail(x_); };
  Ikarus::NonLinearOperator nonLinOp(functions(fvLambda, dfvLambda), parameter(x));

  // Newton method test
  const double eps  = 1e-14;
  const int maxIter = 20;

  TestSuite t("checkNewtonRhapsonFailing");
  Ikarus::NewtonRaphson nr(nonLinOp);
  nr.setup({eps, maxIter});
  const auto solverInfo = nr.solve(1000.0);

  t.check(false == solverInfo.success);
  t.check(maxIter == solverInfo.iterations);

  return t;
}

Eigen::Vector3d fv(const Eigen::Vector3d& x, const Eigen::Matrix3d& A, const Eigen::Vector3d& b) { return b + A * x; }
Eigen::Matrix3d dfv(Eigen::Vector3d&, const Eigen::Matrix3d& A, Eigen::Vector3d&) { return A; }

auto fp(double x, int i) { return 0.5 * x * x + x * i - 2; }
auto dfp(double x, int i) { return x + i; }

auto simple1DOperatorNewtonRhapsonTestWithParamter() {
  TestSuite t("simple1DOperatorNewtonRhapsonTestWithParamter");
  double x = 13;

  for (int i = 0; i < 3; ++i) {
    auto fvLambda  = [](auto&& x_, const int& i_) { return fp(x_, i_); };
    auto dfvLambda = [](auto&& x_, const int& i_) { return dfp(x_, i_); };
    Ikarus::NonLinearOperator nonLinOp(functions(fvLambda, dfvLambda), parameter(x, i));

    // Newton method test
    const double eps       = 1e-14;
    const int maxIter      = 20;
    const double xExpected = std::sqrt(4 + i * i) - i;

    Ikarus::NewtonRaphson nr(nonLinOp);
    const int iterExpected = i == 0 ? 7 : i == 1 ? 5 : 4;
    t.subTest(checkNewtonRhapson(nr, x, eps, maxIter, iterExpected, xExpected, 0.0));
  }
  return t;
}
auto vectorValuedOperatorNewtonRhapsonTest() {
  Eigen::Vector3d x;
  x << 1, 2, 3;
  Eigen::Vector3d b;
  b << 5, 7, 8;
  Eigen::Matrix3d A;
  A = Eigen::Matrix3d::Identity() * 13;

  auto fvLambda  = [&](auto&& x_) { return fv(x_, A, b); };
  auto dfvLambda = [&](auto&& x_) { return dfv(x_, A, b); };
  auto nonLinOp  = Ikarus::NonLinearOperator(functions(fvLambda, dfvLambda), parameter(x));

  // Newton method test
  const double eps  = 1e-14;
  const int maxIter = 20;
  Ikarus::NewtonRaphson nr(nonLinOp, [&](auto& r, auto& A_) { return A_.inverse() * r; });  // special linear solver
  return checkNewtonRhapson(nr, x, eps, maxIter, 1, (-A.ldlt().solve(b)).eval(), Eigen::Vector3d::Zero().eval());
}

double f2v(Eigen::VectorXd& x, Eigen::MatrixXd& A, Eigen::VectorXd& b) { return x.dot(b + A * x); }
Eigen::VectorXd df2v([[maybe_unused]] Eigen::VectorXd& x, Eigen::MatrixXd& A, [[maybe_unused]] Eigen::VectorXd& b) {
  return 2 * A * x + b;
}
Eigen::MatrixXd ddf2v([[maybe_unused]] Eigen::VectorXd& x, Eigen::MatrixXd& A, [[maybe_unused]] Eigen::VectorXd& b) {
  return 2 * A;
}

auto secondOrderVectorValuedOperatorTest() {
  TestSuite t("SecondOrderVectorValuedOperatorTest");
  Eigen::VectorXd x(3);

  x << 1, 2, 3;
  Eigen::VectorXd b(3);
  b << 5, 7, 8;
  Eigen::MatrixXd A(3, 3);
  A = Eigen::MatrixXd::Identity(3, 3) * 13;

  auto fvLambda   = [&](auto&& x_) { return f2v(x_, A, b); };
  auto dfvLambda  = [&](auto&& x_) { return df2v(x_, A, b); };
  auto ddfvLambda = [&](auto&& x_) { return ddf2v(x_, A, b); };
  auto nonLinOp   = Ikarus::NonLinearOperator(functions(fvLambda, dfvLambda, ddfvLambda), parameter(x));

  t.check(checkGradient(nonLinOp, {.draw = false, .writeSlopeStatementIfFailed = false}));

  auto subOperator = nonLinOp.subOperator<1, 2>();
  // Newton method test find root of first derivative
  const double eps  = 1e-14;
  const int maxIter = 20;
  Ikarus::NewtonRaphson nr(subOperator, Ikarus::ILinearSolver<double>(Ikarus::SolverTypeTag::d_LDLT));
  checkNewtonRhapson(nr, x, eps, maxIter, 1, (-0.5 * A.ldlt().solve(b)).eval(), Eigen::VectorXd::Zero(3).eval());
  nonLinOp.update<0>();
  t.check(Dune::FloatCmp::eq(-2.6538461538461533, nonLinOp.value()));
  x << 1, 2, 3;  // Restart and check with predictor
  t.subTest(checkNewtonRhapson(nr, x, eps, maxIter, 1, (-0.5 * A.ldlt().solve(b)).eval(), x));
  nonLinOp.update<0>();
  t.check(Dune::FloatCmp::eq(-2.6538461538461533, nonLinOp.value()));
  return t;
}

#include <autodiff/forward/dual.hpp>
#include <autodiff/forward/dual/eigen.hpp>

using namespace autodiff;
template <typename ScalarType>
ScalarType f2vNL(const Eigen::VectorX<ScalarType>& x, Eigen::MatrixXd&, Eigen::VectorXd&) {
  return x.array().sin().matrix().dot(x);
}

Eigen::VectorXd df2vNL(Eigen::VectorX<autodiff::dual>& x, Eigen::MatrixXd& A, [[maybe_unused]] Eigen::VectorXd& b) {
  return autodiff::gradient(f2vNL<autodiff::dual>, autodiff::wrt(x), autodiff::at(x, A, b));
}

Eigen::MatrixXd ddf2vNL(Eigen::VectorX<autodiff::dual2nd>& x, Eigen::MatrixXd& A, [[maybe_unused]] Eigen::VectorXd& b) {
  return autodiff::hessian(f2vNL<autodiff::dual2nd>, autodiff::wrt(x), autodiff::at(x, A, b));
}

auto secondOrderVectorValuedOperatorNonlinearAutodiff() {
  TestSuite t("SecondOrderVectorValuedOperatorNonlinearAutodiff");
  Eigen::VectorXd x(3);

  x << 1, 2, 3;
  Eigen::VectorXd b(3);
  b << 5, 7, 8;
  Eigen::MatrixXd A(3, 3);
  A = Eigen::MatrixXd::Identity(3, 3) * 13;

  auto fvLambda  = [&](auto&& x_) { return f2vNL<double>(x_, A, b); };
  auto dfvLambda = [&](auto&& x_) {
    auto xR = x_.template cast<autodiff::dual>().eval();
    return df2vNL(xR, A, b);
  };
  auto ddfvLambda = [&](auto&& x_) {
    auto xR = x_.template cast<autodiff::dual2nd>().eval();
    return ddf2vNL(xR, A, b);
  };

  auto nonLinOp = Ikarus::NonLinearOperator(functions(fvLambda, dfvLambda, ddfvLambda), parameter(x));

  t.check(checkGradient(nonLinOp, {.draw = false, .writeSlopeStatementIfFailed = false}));
  t.check(checkHessian(nonLinOp, {.draw = false, .writeSlopeStatementIfFailed = false}));

  auto subOperator = nonLinOp.subOperator<1, 2>();

  // Newton method test find root of first derivative
  const double eps  = 1e-14;
  const int maxIter = 20;
  Ikarus::NewtonRaphson nr(subOperator, Ikarus::ILinearSolver<double>(Ikarus::SolverTypeTag::d_LDLT));

  const Eigen::Vector3d xSol(-4.9131804394348836888, 2.0287578381104342236, 2.0287578381104342236);
  auto nonLinearSolverObserver = std::make_shared<NonLinearSolverLogger>();
  nr.subscribeAll(nonLinearSolverObserver);

  t.subTest(checkNewtonRhapson(nr, x, eps, maxIter, 5, xSol, Eigen::VectorXd::Zero(3).eval()));

  nonLinOp.update<0>();
  t.check(Dune::FloatCmp::eq(-1.1750584073929625716, nonLinOp.value()));
  return t;
}

int main(int argc, char** argv) {
  Ikarus::init(argc, argv);
  TestSuite t;

  t.subTest(simple1DOperatorNewtonRhapsonTest());
  t.subTest(simple1DOperatorNewtonRhapsonCheckThatThePerfectPredictorWorksTest());
  t.subTest(simple1DOperatorNewtonRhapsonWithWrongDerivativeTest());
  t.subTest(simple1DOperatorNewtonRhapsonTestWithParamter());
  t.subTest(vectorValuedOperatorNewtonRhapsonTest());
  t.subTest(secondOrderVectorValuedOperatorTest());
  t.subTest(secondOrderVectorValuedOperatorNonlinearAutodiff());

  return t.exit();
}
