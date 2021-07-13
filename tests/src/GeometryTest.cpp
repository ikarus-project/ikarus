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

#include "ikarus/AnsatzFunctions/Lagrange.h"
#include "ikarus/Geometries/GeometryWithExternalInput.h"

const double tol = 1e-15;

TEST(GeometryTest, CreateJacobianDeterminantAndJacobian2DFlat) {
  using namespace Ikarus::Geometry;
  const Eigen::Matrix<double, 2, 1> xieta({-1.0, -1.0});

  Eigen::Matrix<double, 4, 2> dN = Ikarus::LagrangeCube<double, 2, 1>::evaluateJacobian(xieta);

  Eigen::Matrix<double, 2, 4> x;
  x.col(0) << 0, 0;
  x.col(1) << 4, 0;
  x.col(2) << 0, 4;
  x.col(3) << 4, 4;
  Eigen::Matrix<double, 2, 2> JT = ExternalPlaneGeometry<double>::jacobianTransposed(dN, x);

  Eigen::Matrix2d JTexpected;
  JTexpected << 2, 0, 0, 2;
  EXPECT_EQ(JTexpected, JT);

  Eigen::Matrix2d JTinv = ExternalPlaneGeometry<double>::jacobianInverseTransposed(dN, x);

  Eigen::Matrix2d JTinvexpected;
  JTinvexpected << 0.5, 0.0, 0.0, 0.5;

  EXPECT_EQ(JTinvexpected, JTinv);

  double detJ = ExternalPlaneGeometry<double>::determinantJacobian(dN, x);

  EXPECT_EQ(4.0, detJ);
}

#include <ikarus/Geometries/GeometricElementDefinitions.h>

TEST(GeometryTest, WithInteralAnsatzandVertices) {
  using namespace Ikarus::Geometry;

  const Eigen::Matrix<double, 3, 1> xieta({-1.0, -1.0, -1.0});

  Eigen::Matrix<double, 3, 8> x;
  x.col(0) << 0, 0, 0;
  x.col(1) << 4, 0, 0;
  x.col(2) << 0, 4, 0;
  x.col(3) << 4, 4, 0;
  x.col(4) << 0, 0, 1;
  x.col(5) << 4, 0, 1;
  x.col(6) << 0, 4, 1;
  x.col(7) << 4, 4, 1;
  LinearBrickGeometry<double> geoEle(x);

  Eigen::Matrix<double, 3, 3> JT = geoEle.jacobianTransposed(xieta);

  Eigen::Matrix3d JTexpected;
  JTexpected << 2, 0, 0, 0, 2, 0, 0, 0, 0.5;

  EXPECT_THAT(JT, EigenApproxEqual(JTexpected, tol));

  Eigen::Matrix3d JTinv = geoEle.jacobianInverseTransposed(xieta);

  Eigen::Matrix3d JTinvexpected;
  JTinvexpected << 0.5, 0.0, 0.0, 0.0, 0.5, 0.0, 0.0, 0.0, 2.0;

  EXPECT_THAT(JTinv, EigenApproxEqual(JTinvexpected, tol));

  double detJ = geoEle.determinantJacobian(xieta);

  EXPECT_EQ(2.0, detJ);
}

TEST(GeometryTest, WithInteralAnsatzandVerticesQuadratic) {
  using namespace Ikarus::Geometry;

  const Eigen::Matrix<double, 2, 1> xieta({
      -0.35,
      0.8,
  });

  Eigen::Matrix<double, 2, 9> x;
  x.col(0) << 0, 0;
  x.col(1) << 3, 0;
  x.col(2) << 7, 0;
  x.col(3) << 0, 1;
  x.col(4) << 2, 2;
  x.col(5) << 4, 3;
  x.col(6) << -1, 1;
  x.col(7) << 2, 7;
  x.col(8) << 10, 10;
  QuadraticPlaneGeometry<double> geoEle(x);

  Eigen::Matrix<double, 2, 2> JT = geoEle.jacobianTransposed(xieta);

  Eigen::Matrix2d JTexpected;
  JTexpected << 3.168000000, 4.356000000, -1.033500000, 4.173625000;

  EXPECT_THAT(JT, EigenApproxEqual(JTexpected, tol));

  Eigen::Matrix2d JTinv = geoEle.jacobianInverseTransposed(xieta);

  Eigen::Matrix2d JTinvexpected;
  JTinvexpected << .23547912798317758380, -.24576886555325923030, 0.58310863762464052920e-1, .17874099312964307658;

  EXPECT_THAT(JTinv, EigenApproxEqual(JTinvexpected, tol));

  double detJ = geoEle.determinantJacobian(xieta);

  EXPECT_DOUBLE_EQ(17.72397000, detJ);

  auto dNCart = geoEle.transformCurvLinearDerivativesToCartesian(xieta);

  Eigen::Matrix<double, 9, 2> dNCartExpected;
  dNCartExpected << 0.12624897551105789618e-1, 0.10920614574095834274e-1, -0.10396974453851826743e-1,
      0.88743757452822190506e-1, -0.22279230972539628740e-2, -0.84966090236083447561e-2, -0.56812038979976053284e-1,
      -0.67091418924707847485e-1, 0.46786385042333220342e-1, -.46601333523387006079, 0.10025653937642832933e-1,
      0.46876684807592948154e-1, -.11362407795995210656, .18889292229356298362, 0.93572770084666440706e-1,
      .26796901006332354172, 0.20051307875285665863e-1, -0.61801626009211245234e-1;

  EXPECT_THAT(dNCart, EigenApproxEqual(dNCartExpected, tol));
}

TEST(GeometryTest, CreateJacobianDeterminantAndJacobian2DSurfIn3D) {
  const Eigen::Matrix<double, 2, 1> xieta({-1.0, -1.0});
  Eigen::Matrix<double, 4, 2> dN = Ikarus::LagrangeCube<double, 2, 1>::evaluateJacobian(xieta);

  Eigen::Matrix<double, 2, 4> dNCorrect;
  Eigen::Matrix<double, 3, 4> x;
  x.col(0) << 0, 0, 0;
  x.col(1) << 4, 0, 0;
  x.col(2) << 0, 4, 0;
  x.col(3) << 4, 4, 0;
  Eigen::Matrix<double, 2, 3> JT = Ikarus::Geometry::ExternalSurfaceGeometry<double>::jacobianTransposed(dN, x);
  Eigen::Matrix<double, 2, 3> JTexpected;
  JTexpected << 2, 0, 0, 0, 2, 0;

  ASSERT_EQ(JTexpected, JT);
}
