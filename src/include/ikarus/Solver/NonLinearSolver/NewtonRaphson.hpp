//
// Created by lex on 30/12/2021.
//

#pragma once
#include <iostream>

#include "ikarus/LinearAlgebra/NonLinearOperator.h"
#include "ikarus/utils/Observer/observerMessages.h"
#include "ikarus/utils/LinearAlgebraHelper.h"
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
            typename UpdateFunction = std::function<void(typename NonLinearOperatorImpl::ValueType&,
                                                  const typename NonLinearOperatorImpl::ValueType& b)> >
  class NewtonRaphson : public IObservable<NonLinearSolverMessages> {
  public:
    explicit NewtonRaphson(
        const NonLinearOperatorImpl& p_nonLinearOperator,
        UpdateFunction p_updateFunction = [](typename NonLinearOperatorImpl::ValueType& a,
                                             const typename NonLinearOperatorImpl::ValueType& b) { a += b; })
        : nonLinearOperator{p_nonLinearOperator}, updateFunction{p_updateFunction} {}

    void setup(const NonlinearSolverSettings& p_settings) { settings = p_settings; }
    template <typename SolutionType>
    SolverInformation solve(SolutionType& x) {
      const auto& rx = nonLinearOperator.value();
      const auto& Ax = nonLinearOperator.derivative();
      auto rNorm     = norm(rx);
      int iter{0};
      while (rNorm > settings.tol && iter <= settings.maxIter) {
        this->notify(NonLinearSolverMessages::ITERATION_STARTED);
        if constexpr (std::is_same_v<typename NonLinearOperatorImpl::ValueType, Eigen::VectorXd>) {
          const Eigen::VectorXd D = -Ax.ldlt().solve(rx);
          updateFunction(x, D);
        } else {
          const auto D = -rx / Ax;
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
      solverInformation.iterations=iter;
      solverInformation.residualnorm=rNorm;
      solverInformation.sucess=true;
      return solverInformation;
    }

  private:
    NonLinearOperatorImpl nonLinearOperator;
    UpdateFunction updateFunction;
    NonlinearSolverSettings settings;
  };

}  // namespace Ikarus