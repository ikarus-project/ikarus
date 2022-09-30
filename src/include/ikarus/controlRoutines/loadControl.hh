/*
 * This file is part of the Ikarus distribution (https://github.com/IkarusRepo/Ikarus).
 * Copyright (c) 2022. The Ikarus developers.
 *
 * This library is free software; you can redistribute it and/or
 *  modify it under the terms of the GNU Lesser General Public
 *  License as published by the Free Software Foundation; either
 *  version 2.1 of the License, or (at your option) any later version.
 *
 * This library is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
 * Lesser General Public License for more details.
 *
 * You should have received a copy of the GNU Lesser General Public
 * License along with this library; if not, write to the Free Software
 * Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301 USA
 */

#pragma once
#include <memory>

#include <ikarus/linearAlgebra/nonLinearOperator.hh>
#include <ikarus/solver/nonLinearSolver/newtonRaphson.hh>
#include <ikarus/utils/observer/observer.hh>
#include <ikarus/utils/observer/observerMessages.hh>

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
