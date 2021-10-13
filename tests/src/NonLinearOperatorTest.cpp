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
  auto nonLinOp = NonLinearOperator(&f, &df, x);

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
  auto nonLinOp  = NonLinearOperator(fvLambda, dfvLambda, x);

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