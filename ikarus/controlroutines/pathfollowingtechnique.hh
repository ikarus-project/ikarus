// SPDX-FileCopyrightText: 2021-2024 The Ikarus Developers mueller@ibb.uni-stuttgart.de
// SPDX-License-Identifier: LGPL-3.0-or-later

#pragma once

#include <memory>

#include <Eigen/Core>

#include <ikarus/controlroutines/adaptivestepsizing.hh>
#include <ikarus/controlroutines/controlinfos.hh>
#include <ikarus/controlroutines/pathfollowingfunctions.hh>
#include <ikarus/linearalgebra/nonlinearoperator.hh>
#include <ikarus/solver/nonlinearsolver/newtonraphsonwithscalarsubsidiaryfunction.hh>
#include <ikarus/utils/observer/observer.hh>
#include <ikarus/utils/observer/observermessages.hh>

namespace Ikarus {

  template <typename NonLinearSolver, typename PathFollowingType = Ikarus::StandardArcLength,
            typename AdaptiveStepSizing = Ikarus::AdaptiveStepSizing::NoOp>
  requires((Concepts::PathFollowingStrategy<
               PathFollowingType, typename NonLinearSolver::NonLinearOperator,
               Ikarus::SubsidiaryArgs>)and(Concepts::
                                               AdaptiveStepSizingStrategy<
                                                   AdaptiveStepSizing, Ikarus::NonLinearSolverInformation,
                                                   Ikarus::SubsidiaryArgs,
                                                   std::remove_cvref_t<typename NonLinearSolver::NonLinearOperator>>)
           and (Concepts::NonLinearSolverCheckForPathFollowing<NonLinearSolver>)) class PathFollowing
      : public IObservable<Ikarus::ControlMessages> {
  public:
    PathFollowing(const std::shared_ptr<NonLinearSolver>& p_nonLinearSolver, int loadSteps, double stepSize,
                  PathFollowingType p_pathFollowingType   = Ikarus::StandardArcLength{},
                  AdaptiveStepSizing p_adaptiveStepSizing = {})
        : nonLinearSolver{p_nonLinearSolver},
          loadSteps_{loadSteps},
          stepSize_{stepSize},
          pathFollowingType_{p_pathFollowingType},
          adaptiveStepSizing{p_adaptiveStepSizing} {}

    ControlInformation run() {
      ControlInformation info;
      auto& nonOp = nonLinearSolver->nonLinearOperator();
      this->notify(ControlMessages::CONTROL_STARTED, pathFollowingType_.name);

      SubsidiaryArgs subsidiaryArgs;

      info.totalIterations    = 0;
      subsidiaryArgs.stepSize = stepSize_;
      subsidiaryArgs.DD.resizeLike(nonOp.firstParameter());
      subsidiaryArgs.DD.setZero();

      /// Initializing solver
      this->notify(ControlMessages::STEP_STARTED, 0, subsidiaryArgs.stepSize);
      auto subsidiaryFunctionPtr = [&](auto&& args) { return pathFollowingType_.evaluateSubsidiaryFunction(args); };
      pathFollowingType_.initialPrediction(nonOp, subsidiaryArgs);
      auto solverInfo = nonLinearSolver->solve(subsidiaryFunctionPtr, subsidiaryArgs);
      info.solverInfos.push_back(solverInfo);
      info.totalIterations += solverInfo.iterations;
      if (not solverInfo.success) return info;
      this->notify(ControlMessages::SOLUTION_CHANGED);
      this->notify(ControlMessages::STEP_ENDED);

      /// Calculate predictor for a particular step
      for (int ls = 1; ls < loadSteps_; ++ls) {
        subsidiaryArgs.currentStep = ls;

        adaptiveStepSizing(solverInfo, subsidiaryArgs, nonOp);

        this->notify(ControlMessages::STEP_STARTED, subsidiaryArgs.currentStep, subsidiaryArgs.stepSize);

        pathFollowingType_.intermediatePrediction(nonOp, subsidiaryArgs);

        solverInfo = nonLinearSolver->solve(subsidiaryFunctionPtr, subsidiaryArgs);

        info.solverInfos.push_back(solverInfo);
        info.totalIterations += solverInfo.iterations;
        if (not solverInfo.success) return info;
        this->notify(ControlMessages::SOLUTION_CHANGED);
        this->notify(ControlMessages::STEP_ENDED);
      }

      this->notify(ControlMessages::CONTROL_ENDED, info.totalIterations, pathFollowingType_.name);
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
