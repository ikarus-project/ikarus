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

  struct LoadControlInformation {
    bool sucess{false};
  };

  template <typename NonLinearSolver>
  class LoadControl : public IObservable<ControlMessages> {
  public:
    LoadControl(const std::shared_ptr<NonLinearSolver>& p_nonLinearSolver, int loadSteps, const std::pair<double, double>& tbeginEnd)
        : nonLinearSolver{p_nonLinearSolver},
          loadSteps_{loadSteps},
          parameterBegin_{tbeginEnd.first},
          parameterEnd_{tbeginEnd.second},
          stepSize_{(parameterEnd_ - parameterBegin_) / loadSteps_} {
      static_assert(
          requires {
            nonLinearSolver->nonLinearOperator().lastParameter() = 0.0;
            nonLinearSolver->nonLinearOperator().lastParameter() += 0.0;
          },
          "The last parameter (load factor) must be assignable and incremenntable with a double!");
    }

    LoadControlInformation run() {
      LoadControlInformation info({false});
      auto& nonOp = nonLinearSolver->nonLinearOperator();
      this->notify(ControlMessages::CONTROL_STARTED);
      auto& loadParameter = nonOp.lastParameter();

      loadParameter = 0.0;
      this->notify(ControlMessages::STEP_STARTED);
       auto solverInfo = nonLinearSolver->solve();
      if(not solverInfo.sucess)
        return info;
      this->notify(ControlMessages::SOLUTION_CHANGED);
      this->notify(ControlMessages::STEP_ENDED);

      for (int ls = 0; ls < loadSteps_ ; ++ls) {
        this->notify(ControlMessages::STEP_STARTED);
        loadParameter += stepSize_;
        solverInfo = nonLinearSolver->solve();
        if(not solverInfo.sucess)
          return info;
        this->notify(ControlMessages::SOLUTION_CHANGED);
        this->notify(ControlMessages::STEP_ENDED);
      }
      this->notify(ControlMessages::CONTROL_ENDED);
      info.sucess = true;
      return info;
    }

  private:
    std::shared_ptr<NonLinearSolver> nonLinearSolver;
    int loadSteps_;
    double parameterBegin_;
    double parameterEnd_;
    double stepSize_;
  };
}  // namespace Ikarus
