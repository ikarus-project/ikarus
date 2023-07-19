// SPDX-FileCopyrightText: 2022 The Ikarus Developers mueller@ibb.uni-stuttgart.de
// SPDX-License-Identifier: LGPL-3.0-or-later

#pragma once
#include <memory>

#include <spdlog/spdlog.h>

#include <Eigen/Core>

#include "ikarus/solver/nonLinearSolver/newtonRaphson.hh"
#include <ikarus/controlRoutines/adaptiveStepSizing.hh>
#include <ikarus/linearAlgebra/nonLinearOperator.hh>
#include <ikarus/solver/nonLinearSolver/newtonRaphsonWithScalarSubsidiaryFunction.hh>
#include <ikarus/utils/observer/observer.hh>
#include <ikarus/utils/observer/observerMessages.hh>
#include <ikarus/utils/pathFollowingFunctions.hh>

namespace Ikarus {

  struct PathFollowingInformation {
    bool success{false};
  };

  template <typename NonLinearSolver, typename PathFollowingType = Ikarus::StandardArcLength,
            typename AdaptiveStepSizing = Ikarus::DefaultAdaptiveStepSizing<>>
  requires Concepts::PathFollowingStrategy<PathFollowingType, typename NonLinearSolver::NonLinearOperator>
  class PathFollowing : public IObservable<ControlMessages> {
  public:
    PathFollowing(const std::shared_ptr<NonLinearSolver>& p_nonLinearSolver, int loadSteps, double stepSize,
                  PathFollowingType p_pathFollowingType   = Ikarus::StandardArcLength{},
                  AdaptiveStepSizing p_adaptiveStepSizing = {})
        : nonLinearSolver{p_nonLinearSolver},
          loadSteps_{loadSteps},
          stepSize_{stepSize},
          pathFollowingType_{p_pathFollowingType},
          adaptiveStepSizing{p_adaptiveStepSizing} {}

    PathFollowingInformation run() {
      spdlog::info("Started path following with subsidiary equation: {}", pathFollowingType_.name);
      PathFollowingInformation info({false});
      auto& nonOp = nonLinearSolver->nonLinearOperator();
      this->notify(ControlMessages::CONTROL_STARTED);

      SubsidiaryArgs subsidiaryArgs;

      subsidiaryArgs.stepSize  = stepSize_;
      subsidiaryArgs.loadSteps = loadSteps_;
      subsidiaryArgs.DD.resizeLike(nonOp.firstParameter());
      subsidiaryArgs.DD.setZero();
      /// adaptive step-sizing
      if (std::is_same_v<Ikarus::FancyAdaptiveStepSizing<>, decltype(adaptiveStepSizing)>) {
        auto numEigen = 3;
        subsidiaryArgs.eigenValues.setZero(numEigen, subsidiaryArgs.loadSteps);
      }

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
        subsidiaryArgs.actStep = ls;
        if (ls > 0) pathFollowingType_.intermediatePrediction(nonOp, subsidiaryArgs);

        solverInfo = nonLinearSolver->solve(subsidiaryFunctionPtr, subsidiaryArgs);

        adaptiveStepSizing(solverInfo, nonOp, subsidiaryArgs);

        //
        if (not solverInfo.success) return info;
        this->notify(ControlMessages::SOLUTION_CHANGED);
        this->notify(ControlMessages::STEP_ENDED, subsidiaryArgs.stepSize);
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
    AdaptiveStepSizing adaptiveStepSizing;
  };
}  // namespace Ikarus
