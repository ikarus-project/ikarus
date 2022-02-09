//
// Created by Alex on 21.04.2021.
//

#include <gmock/gmock.h>
#include <gtest/gtest.h>

#include "testHelpers.h"

#include <array>
#include <complex>
#include <vector>

#include <Eigen/Core>

#include "ikarus/Geometries/GeometryWithExternalInput.h"
#include "ikarus/LocalBasis/localBasis.h"
#include <dune/localfunctions/lagrange/lagrangecube.hh>
#include <dune/common/fvector.hh>
#include <dune/common/fmatrix.hh>

const double tol = 1e-15;

TEST(GeometryTest, CreateJacobianDeterminantAndJacobian2DFlat) {
  using namespace Ikarus::Geometry;
  const Dune::FieldVector<double, 2> xieta({0.0, 0.0});

  Dune::Impl::LagrangeCubeLocalBasis<double,double, 2, 1> localBasis_d;

  Ikarus::LocalBasis localBasis(localBasis_d);
  Eigen::Matrix<double,4,2> dN;
  Eigen::Vector4d N;
  localBasis.evaluateFunctionAndJacobian(xieta,N,dN);

  Eigen::Matrix<double, 2, 4> x;
  x.col(0) << 0, 0;
  x.col(1) << 4, 0;
  x.col(2) << 0, 4;
  x.col(3) << 4, 4;
  Eigen::Matrix<double, 2, 2> JT = GeometryWithExternalInput<double,2,2>::jacobianTransposed(dN, x);

  Eigen::Matrix2d JTexpected;
  JTexpected << 4, 0, 0, 4;
  EXPECT_EQ(JTexpected, JT);

  Eigen::Matrix2d JTinv = GeometryWithExternalInput<double,2,2>::jacobianInverseTransposed(dN, x);

  Eigen::Matrix2d JTinvexpected;
  JTinvexpected << 0.25, 0.0, 0.0, 0.25;

  EXPECT_EQ(JTinvexpected, JTinv);

  const double detJ = GeometryWithExternalInput<double,2,2>::determinantJacobian(dN, x);

  EXPECT_EQ(16.0, detJ);
}


TEST(GeometryTest, CreateJacobianDeterminantAndJacobian2DSurfIn3D) {
  const Dune::FieldVector<double, 2> xieta({-1.0, -1.0});
  Dune::Impl::LagrangeCubeLocalBasis<double,double, 2, 1> localBasis_d;

  Ikarus::LocalBasis localBasis(localBasis_d);
  Eigen::Matrix<double,4,2> dN;
  Eigen::Vector4d N;
  localBasis.evaluateFunctionAndJacobian(xieta,N,dN);

  Eigen::Matrix<double, 2, 4> dNCorrect;
  Eigen::Matrix<double, 3, 4> x;
  x.col(0) << 0, 0, 0;
  x.col(1) << 4, 0, 0;
  x.col(2) << 0, 4, 0;
  x.col(3) << 4, 4, 0;
  Eigen::Matrix<double, 2, 3> JT = Ikarus::Geometry::GeometryWithExternalInput<double,2,2>::jacobianTransposed(dN, x);
  Eigen::Matrix<double, 2, 3> JTexpected;
  JTexpected << 4, 0, 0, 0, 4, 0;

  ASSERT_EQ(JTexpected, JT);
}
