// SPDX-FileCopyrightText: 2021-2025 The Ikarus Developers ikarus@ibb.uni-stuttgart.de
// SPDX-License-Identifier: LGPL-3.0-or-later

//
#include <config.h>

#include "testhelpers.hh"

#include <dune/common/test/testsuite.hh>

#include <Eigen/Core>

#include <ikarus/solver/nonlinearsolver/newtonraphson.hh>
#include <ikarus/utils/derivativetraits.hh>
#include <ikarus/utils/differentiablefunction.hh>
#include <ikarus/utils/functionsanitychecks.hh>
#include <ikarus/utils/init.hh>
#include <ikarus/utils/listener/nonlinearsolverlogger.hh>

using namespace Ikarus;
using Dune::TestSuite;

template <typename SolutionType, typename SolutionTypeExpected, typename NewtonRaphson>
auto checkNewtonRaphson(NewtonRaphson& nr, SolutionType& x, double tolerance, int maxIter, int iterExpected,
                        const SolutionTypeExpected& xExpected) {
  TestSuite t("checkNewtonRaphson");
  t.checkThrow<Dune::InvalidStateException>([&]() { nr.setup({tolerance, maxIter, maxIter + 1}); },
                                            "NewtonRaphson setup should fail if minIter > maxIter");
  nr.setup({tolerance, maxIter});
  const auto solverInfo = nr.solve(x);

  if constexpr (std::is_same_v<SolutionType, double>)
    t.check(Dune::FloatCmp::eq(xExpected, x));
  else
    t.check(isApproxSame(x, xExpected, 1e-15));

  t.check(true == solverInfo.success) << "NewtonRaphson wasn't successful.";
  t.check(tolerance >= solverInfo.residualNorm)
      << "Residual norm is not below tolerance " << tolerance << " Actual: " << solverInfo.residualNorm;
  t.check(iterExpected == solverInfo.iterations)
      << "The iteration count does not match the expected number. Expected: " << iterExpected
      << " Actual: " << solverInfo.iterations;
  return t;
}

static auto f(double x) { return 0.5 * x * x + x - 2; }
static auto df(double x) { return x + 1; }

static auto simple1DOperatorNewtonRaphsonTest() {
  TestSuite t("simple1DOperatorNewtonRaphsonTest");

  double x = 13;

  auto fvLambda  = [&](auto&& x_) { return f(x_); };
  auto dfvLambda = [&](auto&& x_) { return df(x_); };
  auto f         = makeDifferentiableFunction(functions(fvLambda, dfvLambda), x);

  // Newton method test
  const double eps       = 1e-14;
  const int maxIter      = 20;
  const double xExpected = std::sqrt(5.0) - 1.0;

  Ikarus::NewtonRaphson nr(f);
  t.subTest(checkNewtonRaphson(nr, x, eps, maxIter, 7, xExpected));
  return t;
}

static auto simple1DOperatorNewtonRaphsonCheckThatThePerfectPredictorWorksTest() {
  TestSuite t("simple1DOperatorNewtonRaphsonCheckThatThePerfectPredictorWorksTest");

  auto fvLambda  = [](auto&& x_) { return f(x_); };
  auto dfvLambda = [](auto&& x_) { return df(x_); };

  const double eps       = 1e-14;
  const int maxIter      = 20;
  const double xExpected = std::sqrt(5.0) - 1.0;
  double x               = xExpected;
  auto f                 = makeDifferentiableFunction(functions(fvLambda, dfvLambda), x);

  Ikarus::NewtonRaphson nr(f);

  t.subTest(checkNewtonRaphson(nr, x, eps, maxIter, 0, xExpected));
  return t;
}

static auto dfFail(double x) { return x + 1000000; }

static auto simple1DOperatorNewtonRaphsonWithWrongDerivativeTest() {
  double x = 13;

  auto fvLambda  = [](auto&& x_) { return f(x_); };
  auto dfvLambda = [](auto&& x_) { return dfFail(x_); };
  auto f         = makeDifferentiableFunction(functions(fvLambda, dfvLambda), x);

  // Newton method test
  const double eps  = 1e-14;
  const int maxIter = 20;

  TestSuite t("checkNewtonRaphsonFailing");
  Ikarus::NewtonRaphson nr(f);
  nr.setup({eps, maxIter});
  const auto solverInfo = nr.solve(x);

  t.check(false == solverInfo.success);
  t.check(maxIter == solverInfo.iterations);

  return t;
}

static Eigen::Vector3d fv(const Eigen::Vector3d& x, const Eigen::Matrix3d& A, const Eigen::Vector3d& b) {
  return b + A * x;
}
static Eigen::Matrix3d dfv(const Eigen::Vector3d&, const Eigen::Matrix3d& A, Eigen::Vector3d&) { return A; }

static auto fp(double x, int i) { return 0.5 * x * x + x * i - 2; }
static auto dfp(double x, int i) { return x + i; }

static auto simple1DOperatorNewtonRaphsonTestWithParameter() {
  TestSuite t("simple1DOperatorNewtonRaphsonTestWithParameter");
  double x = 13;

  for (int i = 0; i < 3; ++i) {
    auto fvLambda  = [&](auto&& x_) { return fp(x_, i); };
    auto dfvLambda = [&](auto&& x_) { return dfp(x_, i); };
    auto f         = makeDifferentiableFunction(functions(fvLambda, dfvLambda), x);

    // Newton method test
    const double eps       = 1e-14;
    const int maxIter      = 20;
    const double xExpected = std::sqrt(4 + i * i) - i;

    Ikarus::NewtonRaphson nr(f);
    const int iterExpected = i == 0 ? 7 : i == 1 ? 5 : 4;
    t.subTest(checkNewtonRaphson(nr, x, eps, maxIter, iterExpected, xExpected));
  }
  return t;
}

static auto vectorValuedOperatorNewtonRaphsonTest() {
  Eigen::Vector3d x;
  x << 1, 2, 3;
  Eigen::Vector3d b;
  b << 5, 7, 8;
  Eigen::Matrix3d A;
  A = Eigen::Matrix3d::Identity() * 13;

  auto fvLambda  = [&](auto&& x_) { return fv(x_, A, b); };
  auto dfvLambda = [&](auto&& x_) { return dfv(x_, A, b); };
  auto f         = makeDifferentiableFunction(functions(fvLambda, dfvLambda), x);

  // Newton method test
  const double eps  = 1e-14;
  const int maxIter = 20;
  Ikarus::NewtonRaphson nr(f, [&](auto& r, auto& A_) { return A_.inverse() * r; }); // special linear solver
  return checkNewtonRaphson(nr, x, eps, maxIter, 1, (-A.ldlt().solve(b)).eval());
}

static double f2v(const Eigen::VectorXd& x, const Eigen::MatrixXd& A, const Eigen::VectorXd& b) {
  return x.dot(b + A * x);
}
static Eigen::VectorXd df2v([[maybe_unused]] const Eigen::VectorXd& x, const Eigen::MatrixXd& A,
                            [[maybe_unused]] const Eigen::VectorXd& b) {
  return 2 * A * x + b;
}

static Eigen::MatrixXd ddf2v([[maybe_unused]] const Eigen::VectorXd& x, const Eigen::MatrixXd& A,
                             [[maybe_unused]] const Eigen::VectorXd& b) {
  return 2 * A;
}

static auto secondOrderVectorValuedOperatorTest() {
  TestSuite t("SecondOrderVectorValuedOperatorTest");
  Eigen::VectorXd x(3);

  x << 1, 2, 3;
  Eigen::VectorXd b(3);
  b << 5, 7, 8;
  Eigen::MatrixXd A(3, 3);
  A = Eigen::MatrixXd::Identity(3, 3) * 13;

  auto fvLambda            = [&](auto&& x_) { return f2v(x_, A, b); };
  auto dfvLambda           = [&](auto&& x_) { return df2v(x_, A, b); };
  auto ddfvLambda          = [&](auto&& x_) { return ddf2v(x_, A, b); };
  DifferentiableFunction f = makeDifferentiableFunction(functions(fvLambda, dfvLambda, ddfvLambda), x);

  t.check(Ikarus::utils::checkGradient(f, x, {.draw = false, .writeSlopeStatementIfFailed = false}));

  auto subOperator = derivative(f);
  // Newton method test find root of first derivative
  const double eps  = 1e-14;
  const int maxIter = 20;
  Ikarus::NewtonRaphson nr(subOperator, Ikarus::LinearSolver(Ikarus::SolverTypeTag::d_LDLT));
  checkNewtonRaphson(nr, x, eps, maxIter, 1, (-0.5 * A.ldlt().solve(b)).eval());
  const double e = f(x);
  t.check(Dune::FloatCmp::eq(-2.6538461538461533, e));
  x << 1, 2, 3; // Restart and check with predictor
  t.subTest(checkNewtonRaphson(nr, x, eps, maxIter, 1, (-0.5 * A.ldlt().solve(b)).eval()));
  const double e2 = f(x);
  t.check(Dune::FloatCmp::eq(-2.6538461538461533, e2));
  return t;
}

#include <autodiff/forward/dual.hpp>
#include <autodiff/forward/dual/eigen.hpp>

using namespace autodiff;
template <typename ScalarType>
ScalarType f2vNL(const Eigen::Vector3<ScalarType>& x, const Eigen::Matrix3d&, const Eigen::Vector3d&) {
  return x.array().sin().matrix().dot(x);
}

static Eigen::Vector3d df2vNL(const Eigen::Vector3d& x, const Eigen::Matrix3d& A,
                              [[maybe_unused]] const Eigen::Vector3d& b) {
  Eigen::Vector3<autodiff::dual> xD = x;
  return autodiff::gradient(f2vNL<autodiff::dual>, autodiff::wrt(xD), autodiff::at(xD, A, b));
}

static Eigen::Matrix3d ddf2vNL(const Eigen::Vector3d& x, const Eigen::Matrix3d& A,
                               [[maybe_unused]] const Eigen::Vector3d& b) {
  Eigen::Vector3<autodiff::dual2nd> xD = x;
  return autodiff::hessian(f2vNL<autodiff::dual2nd>, autodiff::wrt(xD), autodiff::at(xD, A, b));
}

static auto secondOrderVectorValuedOperatorNonlinearAutodiff() {
  TestSuite t("SecondOrderVectorValuedOperatorNonlinearAutodiff");
  Eigen::Vector3d x(3);

  x << 1, 2, 3;
  Eigen::Vector3d b(3);
  b << 5, 7, 8;
  Eigen::Matrix3d A(3, 3);
  A = Eigen::Matrix3d::Identity() * 13;

  auto fvLambda   = [&](auto&& x_) { return f2vNL(x_, A, b); };
  auto dfvLambda  = [&](auto&& x_) { return df2vNL(x_, A, b); };
  auto ddfvLambda = [&](auto&& x_) { return ddf2vNL(x_, A, b); };

  auto f = makeDifferentiableFunction(functions(fvLambda, dfvLambda, ddfvLambda), x);

  t.check(Ikarus::utils::checkGradient(f, x, {.draw = false, .writeSlopeStatementIfFailed = false}));
  t.check(Ikarus::utils::checkHessian(f, x, {.draw = false, .writeSlopeStatementIfFailed = false}));

  auto subOperator = derivative(f);

  // Newton method test find root of first derivative
  const double eps  = 1e-14;
  const int maxIter = 20;
  Ikarus::NewtonRaphson nr(subOperator, [&](auto& r, auto& A_) { return A_.inverse() * r; });
  auto nonLinearSolverObserver = NonLinearSolverLogger().subscribeTo(nr);

  const Eigen::Vector3d xSol(-4.9131804394348836888, 2.0287578381104342236, 2.0287578381104342236);

  t.subTest(checkNewtonRaphson(nr, x, eps, maxIter, 5, xSol));

  const double e = f(x);
  t.check(Dune::FloatCmp::eq(-1.1750584073929625716, e));
  return t;
}

int main(int argc, char** argv) {
  Ikarus::init(argc, argv);
  TestSuite t;

  t.subTest(simple1DOperatorNewtonRaphsonTest());
  t.subTest(simple1DOperatorNewtonRaphsonCheckThatThePerfectPredictorWorksTest());
  t.subTest(simple1DOperatorNewtonRaphsonWithWrongDerivativeTest());
  t.subTest(simple1DOperatorNewtonRaphsonTestWithParameter());
  t.subTest(vectorValuedOperatorNewtonRaphsonTest());
  t.subTest(secondOrderVectorValuedOperatorTest());
  t.subTest(secondOrderVectorValuedOperatorNonlinearAutodiff());

  return t.exit();
}
