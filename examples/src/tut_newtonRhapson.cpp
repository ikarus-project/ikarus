//
// Created by ac120950 on 23.02.2022.
//
#include "ikarus/Solver/NonLinearSolver/NewtonRaphson.hpp"
#include "ikarus/utils/Observer/nonLinearSolverLogger.h"
#include <ikarus/LinearAlgebra/NonLinearOperator.h>

auto f(double& x) { return 0.5 * x * x + x - 2; }
auto df(double& x) { return x + 1; }

void newtonRhapsonVeryBasicExample() {
  double x = 13;

  auto fvLambda  = [&](auto&& x) { return f(x); };
  auto dfvLambda = [&](auto&& x) { return df(x); };
  Ikarus::NonLinearOperator nonLinOp(linearAlgebraFunctions(fvLambda, dfvLambda), parameter(x));

  const double eps       = 1e-10;
  const int maxIter      = 20;
  const double xExpected = std::sqrt(5.0) - 1.0;
  for (int i = 0; i < maxIter; ++i) {
    x -= nonLinOp.value() / nonLinOp.derivative();
    nonLinOp.updateAll();
    if (std::abs(x) < eps) break;
  }

  Ikarus::NewtonRaphson nr(nonLinOp);
  nr.setup({eps, maxIter});
  const auto solverInfo = nr.solve(x);

  std::cout << "success: " << solverInfo.sucess << "\n";
  std::cout << "iterations: " << solverInfo.iterations << "\n";
  std::cout << "residuum: " << solverInfo.residualnorm << "\n";
  std::cout << "solution: " << x << "\n";
  std::cout << "expected solution: " << xExpected << "\n";
}

void newtonRhapsonBasicExampleWithLogger() {
  double x = 13;

  auto fvLambda  = [&](auto&& x) { return f(x); };
  auto dfvLambda = [&](auto&& x) { return df(x); };
  Ikarus::NonLinearOperator nonLinOp(linearAlgebraFunctions(fvLambda, dfvLambda), parameter(x));

  const double eps       = 1e-10;
  const int maxIter      = 20;
  const double xExpected = std::sqrt(5.0) - 1.0;

  Ikarus::NewtonRaphson nr(nonLinOp);
  nr.setup({eps, maxIter});

  // create observer and subscribe to Newton-Rhapson
  auto nonLinearSolverObserver = std::make_shared<NonLinearSolverLogger>();
  nr.subscribe(NonLinearSolverMessages::FINISHED_SUCESSFULLY, nonLinearSolverObserver);
  // nr.subscribeAll(nonLinearSolverObserver);

  const auto solverInfo = nr.solve(x);

  std::cout << "solution: " << x << "\n";
}

int main() {
  // newtonRhapsonVeryBasicExample();
  newtonRhapsonBasicExampleWithLogger();
}