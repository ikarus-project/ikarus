//
// Created by Alex on 21.07.2021.
//
#include <gmock/gmock.h>
#include <gtest/gtest.h>

#include "testHelpers.h"
#define EIGEN_SPARSEMATRIX_PLUGIN "eigenSparseAddon.h"
#include <Eigen/Core>
#include <Eigen/Dense>

#include "ikarus/Assembler/SimpleAssemblers.h"
#include "ikarus/utils/drawing/griddrawer.h"
#include "ikarus/Solver/NonLinearSolver/NewtonRaphson.hpp"
#include "ikarus/utils/Observer/nonLinearSolverLogger.h"
#include <ikarus/LinearAlgebra/NonLinearOperator.h>
#include <ikarus/utils/functionSanityChecks.h>

template <typename SolutionType, typename SolutionTypeExpected, typename NewtonRhapson>
void checkNewtonRhapson(NewtonRhapson& nr, SolutionType& x, double tolerance, int maxIter, int iterExpected,
                        const SolutionTypeExpected& xExpected, const auto& x_Predictor) {
  nr.setup({tolerance, maxIter});
  const auto solverInfo = nr.solve(x_Predictor);

  if constexpr (std::is_same_v<SolutionType, double>)
    EXPECT_DOUBLE_EQ(x, xExpected);
  else
    EXPECT_THAT(x, EigenApproxEqual(xExpected, 1e-15));

  EXPECT_EQ(solverInfo.sucess, true);
  EXPECT_LE(solverInfo.residualnorm, tolerance);
  EXPECT_EQ(solverInfo.iterations, iterExpected);
}

auto f(double& x) { return 0.5 * x * x + x - 2; }
auto df(double& x) { return x + 1; }

TEST(NonLinearOperator, SimpleOperatorNewtonRhapsonTest) {
  double x = 13;

  auto fvLambda  = [&](auto&& x) { return f(x); };
  auto dfvLambda = [&](auto&& x) { return df(x); };
  Ikarus::NonLinearOperator nonLinOp(linearAlgebraFunctions(fvLambda, dfvLambda), parameter(x));

  // Newton method test
  const double eps       = 1e-14;
  const int maxIter      = 20;
  const double xExpected = std::sqrt(5.0) - 1.0;

  Ikarus::NewtonRaphson nr(nonLinOp);

  checkNewtonRhapson(nr, x, eps, maxIter, 7, xExpected, 0.0);
}

Eigen::Vector3d fv(Eigen::Vector3d& x, Eigen::Matrix3d& A, Eigen::Vector3d& b) { return b + A * x; }
Eigen::Matrix3d dfv([[maybe_unused]] Eigen::Vector3d& x, Eigen::Matrix3d& A, [[maybe_unused]] Eigen::Vector3d& b) {
  return A;
}

TEST(NonLinearOperator, VectorValuedOperatorNewtonMethod) {
  Eigen::Vector3d x;
  x << 1, 2, 3;
  Eigen::Vector3d b;
  b << 5, 7, 8;
  Eigen::Matrix3d A;
  A = Eigen::Matrix3d::Identity() * 13;

  auto fvLambda  = [&](auto&& x) { return fv(x, A, b); };
  auto dfvLambda = [&](auto&& x) { return dfv(x, A, b); };
  auto nonLinOp  = Ikarus::NonLinearOperator(linearAlgebraFunctions(fvLambda, dfvLambda), parameter(x));

  // Newton method test
  const double eps  = 1e-14;
  const int maxIter = 20;
  Ikarus::NewtonRaphson nr(nonLinOp, [&](auto& r, auto& A_) { return A.inverse() * r; });  // special linear solver
  checkNewtonRhapson(nr, x, eps, maxIter, 1, (-A.ldlt().solve(b)).eval(), Eigen::Vector3d::Zero().eval());
}

double f2v(Eigen::VectorXd& x, Eigen::MatrixXd& A, Eigen::VectorXd& b) { return x.dot(b + A * x); }
Eigen::VectorXd df2v([[maybe_unused]] Eigen::VectorXd& x, Eigen::MatrixXd& A, [[maybe_unused]] Eigen::VectorXd& b) {
  return 2 * A * x + b;
}
Eigen::MatrixXd ddf2v([[maybe_unused]] Eigen::VectorXd& x, Eigen::MatrixXd& A, [[maybe_unused]] Eigen::VectorXd& b) {
  return 2 * A;
}

TEST(NonLinearOperator, SecondOrderVectorValuedOperator) {
  Eigen::VectorXd x(3);

  x << 1, 2, 3;
  Eigen::VectorXd b(3);
  b << 5, 7, 8;
  Eigen::MatrixXd A(3, 3);
  A = Eigen::MatrixXd::Identity(3, 3) * 13;

  auto fvLambda   = [&](auto&& xL) { return f2v(xL, A, b); };
  auto dfvLambda  = [&](auto&& xL) { return df2v(xL, A, b); };
  auto ddfvLambda = [&](auto&& xL) { return ddf2v(xL, A, b); };
  auto nonLinOp   = Ikarus::NonLinearOperator(linearAlgebraFunctions(fvLambda, dfvLambda, ddfvLambda), parameter(x));

  EXPECT_TRUE(checkGradient(nonLinOp, false));

  auto subOperator = nonLinOp.subOperator<1, 2>();
  // Newton method test find root of first derivative
  const double eps  = 1e-14;
  const int maxIter = 20;
  Ikarus::NewtonRaphson nr(subOperator, Ikarus::ILinearSolver<double>(Ikarus::SolverTypeTag::d_LDLT));
  checkNewtonRhapson(nr, x, eps, maxIter, 1, (-0.5 * A.ldlt().solve(b)).eval(), Eigen::VectorXd::Zero(3).eval());
  nonLinOp.update<0>();
  EXPECT_DOUBLE_EQ(nonLinOp.value(), -2.6538461538461533);
  x << 1, 2, 3;  // Restart and check with predictor
  checkNewtonRhapson(nr, x, eps, maxIter, 2, (-0.5 * A.ldlt().solve(b)).eval(), x);
  nonLinOp.update<0>();
  EXPECT_DOUBLE_EQ(nonLinOp.value(), -2.6538461538461533);
}

#include <autodiff/forward/dual.hpp>
#include <autodiff/forward/dual/eigen.hpp>

using namespace autodiff;
template <typename ScalarType>
ScalarType f2vNL(const Eigen::VectorX<ScalarType>& x, Eigen::MatrixXd&, Eigen::VectorXd&) {
  return x.array().sin().matrix().dot(x);
}

Eigen::VectorXd df2vNL(Eigen::VectorX<autodiff::dual>& x, Eigen::MatrixXd& A, [[maybe_unused]] Eigen::VectorXd& b) {
  return autodiff::gradient(f2vNL<autodiff::dual>, wrt(x), at(x, A, b));
}

Eigen::MatrixXd ddf2vNL(Eigen::VectorX<autodiff::dual2nd>& x, Eigen::MatrixXd& A, [[maybe_unused]] Eigen::VectorXd& b) {
  return autodiff::hessian(f2vNL<autodiff::dual2nd>, wrt(x), at(x, A, b));
}

TEST(NonLinearOperator, SecondOrderVectorValuedOperatorNonlinearAutodiff) {
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

  auto nonLinOp = Ikarus::NonLinearOperator(linearAlgebraFunctions(fvLambda, dfvLambda, ddfvLambda), parameter(x));

  EXPECT_TRUE(checkGradient(nonLinOp, false));
  EXPECT_TRUE(checkHessian(nonLinOp, false));

  auto subOperator = nonLinOp.subOperator<1, 2>();

  // Newton method test find root of first derivative
  const double eps  = 1e-14;
  const int maxIter = 20;
  Ikarus::NewtonRaphson nr(subOperator, Ikarus::ILinearSolver<double>(Ikarus::SolverTypeTag::d_LDLT));

  const Eigen::Vector3d xSol(-4.9131804394348836888, 2.0287578381104342236, 2.0287578381104342236);
  auto nonLinearSolverObserver = std::make_shared<NonLinearSolverLogger>();
  nr.subscribeAll(nonLinearSolverObserver);

  checkNewtonRhapson(nr, x, eps, maxIter, 5, xSol, Eigen::VectorXd::Zero(3).eval());

  nonLinOp.update<0>();
  EXPECT_DOUBLE_EQ(nonLinOp.value(), -1.1750584073929625716);
}
