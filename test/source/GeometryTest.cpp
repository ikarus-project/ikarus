////
//// Created by Alex on 21.04.2021.
////
//
//#include <vector>
//#include <array>
//#include <complex>
//
//#include "spdlog/spdlog.h"
//#define EIGEN_MATRIXBASE_PLUGIN "IBB_Eigen_MatrixBaseAddon.h"
//#include <Eigen/Core>
//
//#include "gtest/gtest.h"
//#include <gmock/gmock.h>
//
//#include "ikarus/Geometries/GeometryWithExternalInput.h"
//#include "ikarus/AnsatzFunctions/Lagrange.h"
//
//#include "testHelpers.h"
//
//const double tol = 1e-15;
//
//TEST(GeometryTest, CreateJacobianDeterminantAndJacobian2DFlat) {
//using namespace Ikarus::Geometry;
//    const Eigen::Matrix< double, 2, 1 > xieta({-1.0,-1.0});
//    auto N = Ikarus::LagrangeCube<double,2,1>::getAnsatzFunctions(xieta);
//    Eigen::Matrix< double, 4, 2 > dN = Ikarus::LagrangeCube<double,2,1>::getAnsatzFunctionJacobian(xieta);
//
//    Eigen::Matrix< double, 2, 4 > x;
//    x.col(0) << 0,0;
//    x.col(1) << 4,0;
//    x.col(2) << 0,4;
//    x.col(3) << 4,4;
//    Eigen::Matrix< double, 2, 2 > JT = PlaneGeometry <double>::getJacobianTransposed(dN,x );
//
//    Eigen::Matrix2d JTexpected;
//    JTexpected<<2,0,
//                0,2;
//    EXPECT_EQ(JTexpected,JT);
//
//    Eigen::Matrix2d JTinv = PlaneGeometry<double>::JacobianInverseTransposed(dN,x );
//
//    Eigen::Matrix2d JTinvexpected;
//    JTinvexpected<<0.5,0.0,
//                   0.0,0.5;
//
//    EXPECT_EQ(JTinvexpected,JTinv);
//
//    double detJ = PlaneGeometry<double>::determinantJacobian(dN,x );
//
//    EXPECT_EQ(4.0,detJ);
//
//}
//
//
//#include <ikarus/Geometries/GeometricElementDefinitions.h>
//
//TEST(GeometryTest,WithInteralAnsatzandVertices) {
//  using namespace Ikarus::Geometry;
//
//
//  const Eigen::Matrix< double, 3, 1 > xieta({-1.0,-1.0,-1.0});
//
//  Eigen::Matrix< double, 3, 8 > x;
//  x.col(0) << 0,0,0;
//  x.col(1) << 4,0,0;
//  x.col(2) << 0,4,0;
//  x.col(3) << 4,4,0;
//  x.col(4) << 0,0,1;
//  x.col(5) << 4,0,1;
//  x.col(6) << 0,4,1;
//  x.col(7) << 4,4,1;
//  BrickGeometry1<double> geoEle(x);
//
//  Eigen::Matrix< double, 3, 3 > JT = geoEle.jacobianTransposed(xieta);
//
//  Eigen::Matrix3d JTexpected;
//  JTexpected<<2,0,0,
//              0,2,0,
//              0,0,0.5;
//
//  EXPECT_THAT(JT,EigenApproxEqual( JTexpected,tol));
//
//  Eigen::Matrix3d JTinv = geoEle.jacobianInverseTransposed(xieta );
//
//  Eigen::Matrix3d JTinvexpected;
//  JTinvexpected<<0.5,0.0,0.0,
//                 0.0,0.5,0.0,
//                 0.0,0.0,2.0;
//
//  EXPECT_THAT(JTinv,EigenApproxEqual(JTinvexpected,tol));
//
//  double detJ = geoEle.determinantJacobian(xieta);
//
//  EXPECT_EQ(2.0,detJ);
//
//}
//
////TEST(GeometryTest, CreateJacobianDeterminantAndJacobian2DSurfIn3D) {
////
////  const Eigen::Matrix< double, 2, 1 > xieta({-1.0,-1.0});
////  auto N = Ikarus::LagrangeCube<double,2,1>::getAnsatzFunction(xieta);
////  Eigen::Matrix< double, 4, 2 > dN = Ikarus::LagrangeCube<double,2,1>::getAnsatzFunctionJacobian(xieta);
////  std::cout<<"N: "<<N<<std::endl;
////  std::cout<<"dN : "<<dN<<std::endl;
//////    Eigen::Matrix< double, 2, 4 > dNCorrect;
////  Eigen::Matrix< double, 3, 4 > x;
////  x.col(0) << 0,0,0;
////  x.col(1) << 4,0,0;
////  x.col(2) << 0,4,0;
////  x.col(3) << 4,4,0;
////  Eigen::Matrix< double, 2, 3 > JT = Ikarus::SurfaceGeometry<double>::getJacobianTransposed(dN,x );
////  std::cout<<JT<<std::endl;
////  std::cout<<"=================="<<std::endl;
////  std::cout<<dN<<std::endl;
////  Eigen::Matrix< double, 2, 3 > JTexpected;
////  JTexpected<<2,0,0,
////      0,2,0;
////
////  ASSERT_EQ(JTexpected,JT);
////}