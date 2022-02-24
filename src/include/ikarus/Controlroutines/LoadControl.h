//
// Created by lex on 17/12/2021.
//

#pragma once
#include <memory>

#include "ikarus/LinearAlgebra/NonLinearOperator.h"
#include "ikarus/Solver/NonLinearSolver/NewtonRaphson.hpp"
#include "ikarus/Variables/ParameterFactory.h"
#include <ikarus/utils/Observer/observer.h>
#include <ikarus/utils/Observer/observerMessages.h>

namespace Ikarus {

  template <typename NonLinearSolver>
  class LoadControl : public IObservable<ControlMessages> {
  public:
    LoadControl(NonLinearSolver&& p_nonLinearSolver, int loadSteps, const std::pair<double, double>& tbeginEnd)
        : nonLinearSolver{std::move(p_nonLinearSolver)},
          loadSteps_{loadSteps},
          parameterBegin_{tbeginEnd.first},
          parameterEnd_{tbeginEnd.second},
          stepSize_{(parameterEnd_ - parameterBegin_) / loadSteps_} {
      static_assert(
          requires {
            nonLinearSolver.nonLinearOperator().lastParameter() = 0.0;
            nonLinearSolver.nonLinearOperator().lastParameter() += 0.0;
          }, "The last parameter (load factor) must be assignable and incremenntable with a double!");
    }

    void run() {
      auto& nonOp = nonLinearSolver.nonLinearOperator();
      this->notify(ControlMessages::CONTROL_STARTED);
      auto& loadParameter = nonOp.lastParameter();
      //      auto& x             = nonOp.template nthParameter<1>();
      // assemble u
      //      auto& PredictRHS    = nonOp.value()(u=(0,0,0,0,uhat));

      loadParameter = 0.0;
      //      auto& x                = feManager_->getVariables();
      //      x+=u=(0,0,0,0,uhat);
      //      nonOp.update<0>();
      //      x += nonOp.deriv().solve(-nonOp.value()),
      //      nonLinearSolver.solve(x);
      //      x += nonOp.value().solve(PredictRHS);
      for (int ls = 0; ls < loadSteps_ + 1; ++ls) {
        this->notify(ControlMessages::STEP_STARTED);
        nonLinearSolver.solve();
        this->notify(ControlMessages::SOLUTION_CHANGED);
        loadParameter += stepSize_;
        this->notify(ControlMessages::STEP_ENDED);
      }
      this->notify(ControlMessages::CONTROL_ENDED);
    }

  private:
    NonLinearSolver nonLinearSolver;
    int loadSteps_;
    double parameterBegin_;
    double parameterEnd_;
    double stepSize_;
  };
}  // namespace Ikarus
