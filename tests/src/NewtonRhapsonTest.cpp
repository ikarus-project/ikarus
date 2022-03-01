//
// Created by Alex on 21.04.2021.
//
#include <gmock/gmock.h>

#include "../../config.h"
#include "ikarus/Solver/NonLinearSolver/NewtonRaphson.hpp"

#include <ikarus/LinearAlgebra/NonLinearOperator.h>
auto f1(double& x) { return 0.5 * x * x + x - 2; }
auto df1(double& x) { return x + 1; }


TEST(NR, NRTest) {

    double x = 13;

    auto fvLambda  = [](auto&& xL) { return f1(xL); };
    auto dfvLambda = [](auto&& xL) { return df1(xL); };
    Ikarus::NonLinearOperator nonLinOp(linearAlgebraFunctions(fvLambda, dfvLambda), parameter(x));

    const double eps       = 1e-10;
    const int maxIter      = 20;
    const double xExpected = std::sqrt(5.0) - 1.0;

    Ikarus::NewtonRaphson nr(nonLinOp);
    nr.setup({eps, maxIter});
    const auto solverInfo = nr.solve();


    EXPECT_EQ(solverInfo.sucess,true);
    EXPECT_EQ(solverInfo.iterations,8);
    EXPECT_LT(solverInfo.iterations,1);
    EXPECT_LT(solverInfo.residualnorm,eps);
    EXPECT_DOUBLE_EQ(x,xExpected);

}
