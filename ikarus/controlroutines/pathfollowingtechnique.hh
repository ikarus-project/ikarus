// SPDX-FileCopyrightText: 2021-2023 The Ikarus Developers mueller@ibb.uni-stuttgart.de
// SPDX-License-Identifier: LGPL-3.0-or-later

#pragma once
#include <memory>

#include <spdlog/spdlog.h>

#include <Eigen/Core>

#include <ikarus/linearalgebra/nonlinearoperator.hh>
#include <ikarus/solver/nonlinearsolver/newtonraphsonwithscalarsubsidiaryfunction.hh>
#include <ikarus/utils/observer/observer.hh>
#include <ikarus/utils/observer/observermessages.hh>
#include <ikarus/utils/pathfollowingfunctions.hh>

namespace Ikarus {

  struct PathFollowingInformation {
    bool success{false};
  };

  template <typename NonLinearSolver, typename PathFollowingType = Ikarus::StandardArcLength>
  requires Concepts::PathFollowingStrategy<PathFollowingType, typename NonLinearSolver::NonLinearOperator>
  class PathFollowing : public IObservable<ControlMessages> {
  public:
    PathFollowing(const std::shared_ptr<NonLinearSolver>& p_nonLinearSolver, int loadSteps, double stepSize,
                  PathFollowingType p_pathFollowingType = Ikarus::StandardArcLength{})
        : nonLinearSolver{p_nonLinearSolver},
          loadSteps_{loadSteps},
          stepSize_{stepSize},
          pathFollowingType_{p_pathFollowingType} {}

    PathFollowingInformation run() {
      spdlog::info("Started path following with subsidiary equation: {}", pathFollowingType_.name);
      PathFollowingInformation info({false});
      auto& nonOp = nonLinearSolver->nonLinearOperator();
      this->notify(ControlMessages::CONTROL_STARTED);

      SubsidiaryArgs subsidiaryArgs;

      subsidiaryArgs.stepSize = stepSize_;
      subsidiaryArgs.DD.resizeLike(nonOp.firstParameter());
      subsidiaryArgs.DD.setZero();

      /// Initializing solver
      this->notify(ControlMessages::STEP_STARTED);
      auto subsidiaryFunctionPtr = [&](auto&& args) { return pathFollowingType_.evaluateSubsidiaryFunction(args); };
      pathFollowingType_.initialPrediction(nonOp, subsidiaryArgs);
      auto solverInfo = nonLinearSolver->solve(subsidiaryFunctionPtr, subsidiaryArgs);
      if (not solverInfo.success) return info;
      this->notify(ControlMessages::SOLUTION_CHANGED);
      this->notify(ControlMessages::STEP_ENDED);

      /// Calculate predictor for a particular step
      for (int ls = 0; ls < loadSteps_; ++ls) {
        this->notify(ControlMessages::STEP_STARTED);
        if (ls > 0) pathFollowingType_.intermediatePrediction(nonOp, subsidiaryArgs);

        solverInfo = nonLinearSolver->solve(subsidiaryFunctionPtr, subsidiaryArgs);

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
    double stepSize_;
    PathFollowingType pathFollowingType_;
  };
}  // namespace Ikarus
