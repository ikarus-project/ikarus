//
// Created by Alex on 21.04.2021.
//

#include <gmock/gmock.h>

#include "testHelpers.hh"

#include <ikarus/linearAlgebra/nonLinearOperator.hh>
#include <ikarus/manifolds/realTuple.hh>
#include <ikarus/manifolds/unitVector.hh>
#include <ikarus/utils/linearAlgebraHelper.hh>


#include <autodiff/forward/dual.hpp>
#include <autodiff/forward/dual/eigen.hpp>

TEST(UnitQuaternion, RotationMatrixDerivatives) {
  Eigen::Vector<autodiff::dual , 4> q;
  q.setRandom();
  q.normalize();

   const auto D = Ikarus::rotationMatrixColumnDerivatives(q.cast<double>().eval());
     std::array<Eigen::Matrix<double,3,3>,4> RA;
    for (int i = 0; i < 4; ++i) {
        q[i].grad=1;
        const Eigen::Matrix<autodiff::dual,3,3> R =Ikarus::toRotationMatrix(q);
        for (int j = 0; j < 3; ++j) {
            for (int k = 0; k < 3; ++k) {
                RA[i](j,k) = R(j,k).grad;
            }
        }
        q[i].grad=0;
    }

    EXPECT_THAT(D[0].col(0), EigenApproxEqual(RA[0].col(0), 1e-15));
    EXPECT_THAT(D[0].col(1), EigenApproxEqual(RA[1].col(0), 1e-15));
    EXPECT_THAT(D[0].col(2), EigenApproxEqual(RA[2].col(0), 1e-15));
    EXPECT_THAT(D[0].col(3), EigenApproxEqual(RA[3].col(0), 1e-15));

    EXPECT_THAT(D[1].col(0), EigenApproxEqual(RA[0].col(1), 1e-15));
    EXPECT_THAT(D[1].col(1), EigenApproxEqual(RA[1].col(1), 1e-15));
    EXPECT_THAT(D[1].col(2), EigenApproxEqual(RA[2].col(1), 1e-15));
    EXPECT_THAT(D[1].col(3), EigenApproxEqual(RA[3].col(1), 1e-15));

    EXPECT_THAT(D[2].col(0), EigenApproxEqual(RA[0].col(2), 1e-15));
    EXPECT_THAT(D[2].col(1), EigenApproxEqual(RA[1].col(2), 1e-15));
    EXPECT_THAT(D[2].col(2), EigenApproxEqual(RA[2].col(2), 1e-15));
    EXPECT_THAT(D[2].col(3), EigenApproxEqual(RA[3].col(2), 1e-15));
}

