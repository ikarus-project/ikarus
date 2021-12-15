//
// Created by Alex on 21.07.2021.
//
#include <gmock/gmock.h>
#include <gtest/gtest.h>

#include "testHelpers.h"

#include <Eigen/Core>
#include <Eigen/Dense>

#include <ikarus/LinearAlgebra/NonLinearOperator.h>

auto f(double& x) { return 0.5 * x * x + x - 2; }
auto df(double& x) { return x + 1; }

TEST(NonLinearOperator, SimpleOperator) {
  double x      = 13;

  Ikarus::NonLinearOperator nonLinOp(&f, derivatives(&df), parameter(x));

  auto& val      = nonLinOp.value();
  auto& gradient = nonLinOp.derivative();

  // Newton method test
  const double eps  = 1e-14;
  const int maxIter = 20;
  int iter          = 0;
  while (val > eps && iter < maxIter) {
    x -= val / gradient;
    nonLinOp.updateAll();
    ++iter;
  }
  EXPECT_DOUBLE_EQ(val, 0.0);
  double xExpected = std::sqrt(5.0) - 1.0;
  EXPECT_DOUBLE_EQ(gradient, df(xExpected));
  EXPECT_DOUBLE_EQ(x, xExpected);
}

Eigen::VectorXd fv(Eigen::VectorXd& x, Eigen::MatrixXd& A, Eigen::VectorXd& b) { return b + A * x; }
Eigen::MatrixXd dfv([[maybe_unused]] Eigen::VectorXd& x, Eigen::MatrixXd& A, [[maybe_unused]] Eigen::VectorXd& b) {
  return A;
}

TEST(NonLinearOperator, VectorValuedOperator) {
  Eigen::VectorXd x(3);

  x << 1, 2, 3;
  Eigen::VectorXd b(3);
  b << 5, 7, 8;
  Eigen::MatrixXd A(3, 3);
  A = Eigen::MatrixXd::Identity(3, 3) * 13;

  auto fvLambda  = [&](auto&& x) { return fv(x, A, b); };
  auto dfvLambda = [&](auto&& x) { return dfv(x, A, b); };
  auto nonLinOp  = Ikarus::NonLinearOperator(fvLambda, derivatives(dfvLambda), parameter(x));

  auto& val      = nonLinOp.value();
  auto& jacobian = nonLinOp.derivative();

  // Newton method test
  const double eps  = 1e-14;
  const int maxIter = 20;
  int iter          = 0;
  while (val.norm() > eps && iter < maxIter) {
    x -= jacobian.inverse() * val;
    nonLinOp.updateAll();
    ++iter;
  }
  EXPECT_EQ(iter, 1);  // Linear System should be solved in one step
  EXPECT_THAT(b, EigenApproxEqual(-A * x, 1e-15));
}


double f2v(Eigen::VectorXd& x, Eigen::MatrixXd& A, Eigen::VectorXd& b) { return x.dot(b + A * x); }
Eigen::VectorXd df2v([[maybe_unused]] Eigen::VectorXd& x, Eigen::MatrixXd& A, [[maybe_unused]] Eigen::VectorXd& b) {
  return 2*A*x + b;
}
Eigen::MatrixXd ddf2v([[maybe_unused]] Eigen::VectorXd& x, Eigen::MatrixXd& A, [[maybe_unused]] Eigen::VectorXd& b) {
  return 2*A;
}

TEST(NonLinearOperator, SecondOrderVectorValuedOperator) {
  Eigen::VectorXd x(3);

  x << 1, 2, 3;
  Eigen::VectorXd b(3);
  b << 5, 7, 8;
  Eigen::MatrixXd A(3, 3);
  A = Eigen::MatrixXd::Identity(3, 3) * 13;

  auto fvLambda  = [&](auto&& x) { return f2v(x, A, b); };
  auto dfvLambda = [&](auto&& x) { return df2v(x, A, b); };
  auto ddfvLambda = [&](auto&& x) { return ddf2v(x, A, b); };
  auto nonLinOp  = Ikarus::NonLinearOperator(fvLambda, derivatives(dfvLambda,ddfvLambda), parameter(x));

  auto& val      = nonLinOp.value();
  auto& residual = nonLinOp.derivative();
  auto& hessian = nonLinOp.secondDerivative();

  // Newton method test find root of first derivative
  const double eps  = 1e-14;
  const int maxIter = 20;
  int iter          = 0;
  while (residual.norm() > eps && iter < maxIter) {
    x -= hessian.inverse() * residual;
    nonLinOp.updateAll();
    std::cout<<val<<std::endl;
    ++iter;
  }

  EXPECT_EQ(iter, 1);  // Linear System should be solved in one step
  EXPECT_EQ(val, -2.6538461538461533);  // Linear System should be solved in one step
  EXPECT_THAT(b, EigenApproxEqual(-2*A * x, 1e-15));
}