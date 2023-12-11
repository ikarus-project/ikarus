// SPDX-FileCopyrightText: 2021-2023 The Ikarus Developers mueller@ibb.uni-stuttgart.de
// SPDX-License-Identifier: LGPL-3.0-or-later

#pragma once
#include <memory>

#include <ikarus/linearalgebra/nonlinearoperator.hh>
#include <ikarus/solver/nonlinearsolver/newtonraphson.hh>
#include <ikarus/utils/observer/observer.hh>
#include <ikarus/utils/observer/observermessages.hh>

namespace Ikarus {

  struct LoadControlInformation {
    bool success{false};
  };

  /**  The loadControl control routine simply increases the last parameter of a nonlinear operator and then calls
   * a nonlinear solver, e.g. Newton's method */
  template <typename NonLinearSolver>
  class LoadControl : public IObservable<ControlMessages> {
  public:
    LoadControl(const std::shared_ptr<NonLinearSolver>& p_nonLinearSolver, int loadSteps,
                const std::array<double, 2>& tbeginEnd)
        : nonLinearSolver{p_nonLinearSolver},
          loadSteps_{loadSteps},
          parameterBegin_{tbeginEnd[0]},
          parameterEnd_{tbeginEnd[1]},
          stepSize_{(parameterEnd_ - parameterBegin_) / loadSteps_} {
      static_assert(
          requires {
            nonLinearSolver->nonLinearOperator().lastParameter() = 0.0;
            nonLinearSolver->nonLinearOperator().lastParameter() += 0.0;
          },
          "The last parameter (load factor) must be assignable and incrementable with a double!");
    }

    LoadControlInformation run() {
      LoadControlInformation info({false});
      auto& nonOp = nonLinearSolver->nonLinearOperator();
      this->notify(ControlMessages::CONTROL_STARTED);
      auto& loadParameter = nonOp.lastParameter();

      loadParameter = 0.0;
      this->notify(ControlMessages::STEP_STARTED);
      auto solverInfo = nonLinearSolver->solve();
      if (not solverInfo.success) return info;
      this->notify(ControlMessages::SOLUTION_CHANGED);
      this->notify(ControlMessages::STEP_ENDED);

      for (int ls = 0; ls < loadSteps_; ++ls) {
        this->notify(ControlMessages::STEP_STARTED);
        loadParameter += stepSize_;
        solverInfo = nonLinearSolver->solve();
        if (not solverInfo.success) return info;
        this->notify(ControlMessages::SOLUTION_CHANGED);
        this->notify(ControlMessages::STEP_ENDED);
      }
      this->notify(ControlMessages::CONTROL_ENDED);
      info.success = true;
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
