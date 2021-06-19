//
// Created by Alex on 21.04.2021.
//

#include <vector>
#include <array>
//#define EIGEN_MATRIXBASE_PLUGIN "IBB_Eigen_MatrixBaseAddon.h"
#include <Eigen/Core>
#include <ikarus/utils/LinearAlgebraTypedefs.h>


#include "gtest/gtest.h"

using FunctionSignature = double(const Ikarus::FixedVector2d&);
TEST(SimpleShapeFunction, 1DEvaluation) {

    Ikarus::FixedMatrixd<2,4> corners;
    corners.col(0)<<-1.0,-1.0;
    corners.col(1)<<1.0,-1.0;
    corners.col(2)<<1.0,1.0;
    corners.col(3)<<-1.0,1.0;

    const std::function<FunctionSignature > N1 = [](const Ikarus::FixedVector2d& paraPoint )->double { return 0.25*(1-paraPoint(0))*(1-paraPoint(1));};
    const std::function<FunctionSignature > N2 = [](const Ikarus::FixedVector2d& paraPoint )->double { return 0.25*(1+paraPoint(0))*(1-paraPoint(1));};
    const std::function<FunctionSignature > N3 = [](const Ikarus::FixedVector2d& paraPoint )->double { return 0.25*(1-paraPoint(0))*(1+paraPoint(1));};
    const std::function<FunctionSignature > N4 = [](const Ikarus::FixedVector2d& paraPoint )->double { return 0.25*(1+paraPoint(0))*(1+paraPoint(1));};
    using AnsatzFunctionQ1 = std::array<const std::function<double(const Ikarus::FixedVector2d&) >,4>;
    AnsatzFunctionQ1 N{N1,N2,N3,N4};
    EXPECT_EQ(N[0](corners.col(0)), 1.0);
    EXPECT_EQ(N[0](corners.col(1)), 0.0);
    EXPECT_EQ(N[0](corners.col(2)), 0.0);
    EXPECT_EQ(N[0](corners.col(3)), 0.0);
}
