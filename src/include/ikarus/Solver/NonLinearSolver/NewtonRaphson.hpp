//
// Created by lex on 30/12/2021.
//

#pragma once
#include <iostream>

#include "ikarus/LinearAlgebra/NonLinearOperator.h"
#include "ikarus/Solver/LinearSolver/LinearSolver.h"
#include "ikarus/utils/LinearAlgebraHelper.h"
#include "ikarus/utils/Observer/observerMessages.h"
#include <ikarus/utils/Observer/observer.h>

namespace Ikarus {

  struct NonlinearSolverSettings {
    double tol{1e-8};
    int maxIter{20};
  };

  struct SolverInformation {
    bool sucess{false};
    double residualnorm{0.0};
    int iterations{0};
  };
  template <typename NonLinearOperatorImpl,
            typename LinearSolver   = std::function<typename NonLinearOperatorImpl::ValueType(
                  const typename NonLinearOperatorImpl::ValueType&, const typename NonLinearOperatorImpl::ValueType&)>,
            typename UpdateFunction = std::function<void(typename NonLinearOperatorImpl::ValueType&,
                                                         const typename NonLinearOperatorImpl::ValueType&)> >
  class NewtonRaphson : public IObservable<NonLinearSolverMessages> {
  public:
    using LinearSolverStdFunctionType = std::function<typename NonLinearOperatorImpl::ValueType(
        const typename NonLinearOperatorImpl::ValueType&, const typename NonLinearOperatorImpl::ValueType&)>;
    //    explicit NewtonRaphson(
    //        const NonLinearOperatorImpl& p_nonLinearOperator, LinearSolver p_linearSolver=[](const typename
    //        NonLinearOperatorImpl::ValueType& a,
    //                                                                                            const typename
    //                                                                                            NonLinearOperatorImpl::ValueType&
    //                                                                                            b) { return a / b; },
    //        UpdateFunction p_updateFunction = [](typename NonLinearOperatorImpl::ValueType& a,
    //                                             const typename NonLinearOperatorImpl::ValueType& b) { a += b; })
    //        : nonLinearOperator{p_nonLinearOperator}, linearSolver{std::move(p_linearSolver)},
    //        updateFunction{p_updateFunction} {}

    explicit NewtonRaphson(
        const NonLinearOperatorImpl& p_nonLinearOperator,
        LinearSolver&& p_linearSolver   = [](const typename NonLinearOperatorImpl::ValueType& a,
                                           const typename NonLinearOperatorImpl::ValueType& b) { return a / b; },
        UpdateFunction p_updateFunction = [](typename NonLinearOperatorImpl::ValueType& a,
                                             const typename NonLinearOperatorImpl::ValueType& b) { a += b; })
        : nonLinearOperator{p_nonLinearOperator},
          linearSolver{std::move(p_linearSolver)},
          updateFunction{p_updateFunction} {}

    void setup(const NonlinearSolverSettings& p_settings) { settings = p_settings; }
    template <typename SolutionType>
    SolverInformation solve(SolutionType& x) {
      std::cout<<"Solve:"<<std::endl;
      nonLinearOperator.updateAll();
      const auto& rx = nonLinearOperator.value();
      const auto& Ax = nonLinearOperator.derivative();
      auto rNorm     = norm(rx);
      int iter{0};
      if constexpr (!std::is_same_v<LinearSolver, LinearSolverStdFunctionType>)
        linearSolver.analyzePattern(Ax);
      while (rNorm > settings.tol && iter <= settings.maxIter) {
        this->notify(NonLinearSolverMessages::ITERATION_STARTED);
        if constexpr (!std::is_same_v<LinearSolver, LinearSolverStdFunctionType>) {
          linearSolver.compute(Ax);
          const Eigen::VectorXd D = -linearSolver.solve(rx);
          updateFunction(x, D);
        } else {
          const auto D = -linearSolver(rx, Ax);
          updateFunction(x, D);
        }
        this->notify(NonLinearSolverMessages::SOLUTION_CHANGED);
        nonLinearOperator.updateAll();
        rNorm = norm(rx);
        this->notify(NonLinearSolverMessages::RESIDUALNORM_UPDATED, rNorm);
        this->notify(NonLinearSolverMessages::ITERATION_ENDED);
        ++iter;
      }
      SolverInformation solverInformation;
      solverInformation.iterations   = iter;
      solverInformation.residualnorm = rNorm;
      solverInformation.sucess       = true;
      return solverInformation;
    }

  private:
    NonLinearOperatorImpl nonLinearOperator;
    LinearSolver linearSolver;
    UpdateFunction updateFunction;
    NonlinearSolverSettings settings;
  };

}  // namespace Ikarus