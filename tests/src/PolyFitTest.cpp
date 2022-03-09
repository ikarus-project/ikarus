//
// Created by Alex on 21.04.2021.
//
#include <gmock/gmock.h>
#include <gtest/gtest.h>

#include "testHelpers.h"

#include <Eigen/Core>
#include <ikarus/utils/utils/polyfit.h>
#include <matplot/matplot.h>


TEST(PolyFitTest, PolyFitTest1) {

  Eigen::VectorXd x = Eigen::VectorXd::LinSpaced(10,0,10);
  Eigen::VectorXd y = Eigen::VectorXd::LinSpaced(10,2,20);

  auto [poly,normE]= Ikarus::polyfit(x,y,1);
  EXPECT_DOUBLE_EQ(poly.coefficients()[0],2.0);
  EXPECT_DOUBLE_EQ(poly.coefficients()[1],1.8);
  EXPECT_LT(normE,1e-14);
}


TEST(PolyFitTest, PolyFitTest2) {

  Eigen::VectorXd x = Eigen::VectorXd::LinSpaced(10,0,10);
  Eigen::VectorXd y = 7*x.array().cwiseProduct(x.array()).matrix()+x+Eigen::VectorXd::Random(10)*10;

  auto [poly,normE]= Ikarus::polyfit(x,y,2);

  EXPECT_DOUBLE_EQ(poly.coefficients()[0],-1.6516491022370672);
  EXPECT_DOUBLE_EQ(poly.coefficients()[1],0.1378228751314157);
  EXPECT_DOUBLE_EQ(poly.coefficients()[2],7.0960748007968437);
  EXPECT_DOUBLE_EQ(normE,12.069293754483319);
}
