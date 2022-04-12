//
// Created by Alex on 21.04.2021.
//
#include <gmock/gmock.h>
#include <gtest/gtest.h>

#include "testHelpers.hh"

#include <matplot/matplot.h>

#include <Eigen/Core>

#include <ikarus/utils/polyfit.hh>

TEST(PolyFitTest, PolyFitTest1) {
  Eigen::VectorXd x = Eigen::VectorXd::LinSpaced(10, 0, 10);
  Eigen::VectorXd y = Eigen::VectorXd::LinSpaced(10, 2, 20);

  auto [poly, normE] = Ikarus::polyfit(x, y, 1);
  EXPECT_DOUBLE_EQ(poly.coefficients()[0], 2.0);
  EXPECT_DOUBLE_EQ(poly.coefficients()[1], 1.8);
  EXPECT_LT(normE, 1e-14);
}

TEST(PolyFitTest, PolyFitTest2) {
  const double factor = 7.6;
  Eigen::VectorXd x   = Eigen::VectorXd::LinSpaced(10, 0, 10);
  Eigen::VectorXd y   = 7 * x.array().cwiseProduct(x.array()).matrix();
  for (int i = 0; i < y.size(); ++i) {
    y[i] += (1 - i / 10.0) * factor - (1 - i * i / 10.0) * factor + std::sin(i / 10.0);
  }

  auto [poly, normE] = Ikarus::polyfit(x, y, 2);

  EXPECT_DOUBLE_EQ(poly.coefficients()[0], -0.0038062785674569739);
  EXPECT_DOUBLE_EQ(poly.coefficients()[1], -0.58760441700969401);
  EXPECT_DOUBLE_EQ(poly.coefficients()[2], 7.6138682871655829);
  EXPECT_DOUBLE_EQ(normE, 0.0082367593944499204);
}
