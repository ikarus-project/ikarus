//
// Created by ac120950 on 23.02.2022.
//
#include "ikarus/Solver/NonLinearSolver/NewtonRaphson.hpp"
#include "ikarus/utils/Observer/nonLinearSolverLogger.h"
#include <ikarus/LinearAlgebra/NonLinearOperator.h>

auto f(double& x) { return 0.5 * x * x + x - 2; }
auto df(double& x) { return x + 1; }

void newtonRaphsonVeryBasicExample() {
  double x               = 13;
  const double eps       = 1e-10;
  const int maxIter      = 20;
  const double xExpected = std::sqrt(5.0) - 1.0;

  auto fvLambda  = [&](auto&& x) { return f(x); };
  auto dfvLambda = [&](auto&& x) { return df(x); };
  Ikarus::NonLinearOperator nonLinOp(linearAlgebraFunctions(fvLambda, dfvLambda), parameter(x));

  int iterCount = 1;
  while (abs(nonLinOp.value()) > eps and iterCount <= maxIter) {
    x -= nonLinOp.value() / nonLinOp.derivative();
    nonLinOp.updateAll();
    iterCount++;

    std::cout << "nonlinearOperator, value(): " << nonLinOp.value() << "\n";
    std::cout << "nonlinearOperator, x: " << nonLinOp.firstParameter() << "\n";
  }

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

class OurFirstObserver : public IObserver<NonLinearSolverMessages> {
public:
  void updateImpl(NonLinearSolverMessages message) override {
    if (message == NonLinearSolverMessages::ITERATION_STARTED) std::cout << "Yeah, the iteration started. Let's go!\n";
  }
};

void newtonRaphsonBasicExampleWithLogger() {
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
  auto ourSimpleObserver = std::make_shared<OurFirstObserver>();
  nr.subscribe(NonLinearSolverMessages::ITERATION_STARTED, ourSimpleObserver);
  // nr.subscribeAll(ourSimpleObserver);
  // auto nonLinearSolverObserver = std::make_shared<NonLinearSolverLogger>();
  // nr.subscribe(NonLinearSolverMessages::FINISHED_SUCESSFULLY, nonLinearSolverObserver);
  // nr.subscribeAll(nonLinearSolverObserver);

  const auto solverInfo = nr.solve(x);

  std::cout << "solution: " << x << "\n";
}

int main() {
  newtonRaphsonVeryBasicExample();
  newtonRaphsonBasicExampleWithLogger();
}