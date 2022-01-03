//
// Created by lex on 30/12/2021.
//

#pragma once
#include "ikarus/LinearAlgebra/NonLinearOperator.h"
#include <ikarus/utils/Observer/observer.h>
#include "ikarus/utils/Observer/observerMessages.h"

namespace Ikarus {

  struct NonlinearSolverSettings
  {
    double tol{1e-8};
    int maxIter{20};
  };
  template <typename... LinearAlgebraFunctionArgs>
  class NewtonRaphson : public IObservable<NonLinearSolverMessages> {

    NewtonRaphson(const nonLinearOperator&)
        : linearAlgebraFunctions_{linearAlgebraFunctions} {}

        void setup(const NonlinearSolverSettings& settings) { settings_= settings;}
        template<typename SolutionType>
        void solve(SolutionType& x)
        {
          auto rNorm = nonLinearOperator.value().norm();
          int iter   = 0;
          while (rNorm > settings_.tol && iter <= settings_.maxIter) {
            const auto& residual = nonLinearOperator.value();
            const auto& K        = nonLinearOperator.derivative();

            const Eigen::VectorXd D = -K.ldlt().solve(residual);
            x += D;
            nonLinearOperator.updateAll();
            rNorm = residual.norm();
            this->notify(ControlMessages::RESIDUALNORM_UPDATED, rNorm);
            this->notify(ControlMessages::SOLUTION_CHANGED);
            this->notify(ControlMessages::ITERATION_ENDED);
            ++iter;
          }
        }
        LinearAlgebraFunctions<LinearAlgebraFunctionArgs...> linearAlgebraFunctions_;
        NonlinearSolverSettings settings_;
  };

}  // namespace Ikarus