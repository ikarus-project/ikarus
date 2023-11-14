// SPDX-FileCopyrightText: 2022 The Ikarus Developers mueller@ibb.uni-stuttgart.de
// SPDX-License-Identifier: LGPL-3.0-or-later

#pragma once

#include "adaptivestepsizing.hh"

#include <chrono>
#include <memory>

#include <Eigen/Core>

#include <ikarus/linearalgebra/nonlinearoperator.hh>
#include <ikarus/solver/nonlinearsolver/newtonraphsonwithscalarsubsidiaryfunction.hh>
#include <ikarus/utils/observer/observer.hh>
#include <ikarus/utils/observer/observermessages.hh>
#include <ikarus/utils/pathfollowingfunctions.hh>

namespace Ikarus {

  struct ControlInformation {
    bool success{false};
    std::vector<Ikarus::NonLinearSolverInformation> solverInfos{};
  };

  template <typename NonLinearSolver, typename PathFollowingType = Ikarus::StandardArcLength,
            typename AdaptiveStepSizing = Ikarus::DefaultAdaptiveStepSizing<>,
            typename ActivationFunction = Ikarus::AdaptiveStepSizing::LinearPiecewiseFunction>
    requires Concepts::PathFollowingStrategy<PathFollowingType, typename NonLinearSolver::NonLinearOperator>
  class PathFollowing : public IObservable<Ikarus::ControlMessages> {
  public:
    PathFollowing(const std::shared_ptr<NonLinearSolver>& p_nonLinearSolver, int loadSteps, double stepSize,
                  PathFollowingType p_pathFollowingType   = Ikarus::StandardArcLength{},
                  AdaptiveStepSizing p_adaptiveStepSizing = {},
                  const ActivationFunction& p_activationFunction
                  = Ikarus::AdaptiveStepSizing::LinearPiecewiseFunction{})
        : nonLinearSolver{p_nonLinearSolver},
          loadSteps_{loadSteps},
          stepSize_{stepSize},
          pathFollowingType_{p_pathFollowingType},
          adaptiveStepSizing{p_adaptiveStepSizing},
          activationFunction{p_activationFunction} {}

    ControlInformation run() {
      ControlInformation info({false});
      auto& nonOp = nonLinearSolver->nonLinearOperator();
      this->notify(ControlMessages::CONTROL_STARTED, pathFollowingType_.name);

      auto start = std::chrono::high_resolution_clock::now();

      SubsidiaryArgs subsidiaryArgs;

      subsidiaryArgs.totalIte  = 0;
      subsidiaryArgs.stepSize  = stepSize_;
      subsidiaryArgs.loadSteps = loadSteps_;
      subsidiaryArgs.DD.resizeLike(nonOp.firstParameter());
      subsidiaryArgs.DD.setZero();
      /// adaptive step-sizing
      if (std::is_same_v<Ikarus::FancyAdaptiveStepSizing<ActivationFunction>, decltype(adaptiveStepSizing)>) {
        auto numEigen = 3;
        subsidiaryArgs.eigenValues.setZero(numEigen, subsidiaryArgs.loadSteps);
      }

      /// Initializing solver
      this->notify(ControlMessages::STEP_STARTED, 0, subsidiaryArgs.stepSize);
      auto subsidiaryFunctionPtr = [&](auto&& args) { return pathFollowingType_.evaluateSubsidiaryFunction(args); };
      pathFollowingType_.initialPrediction(nonOp, subsidiaryArgs);
      auto solverInfo = nonLinearSolver->solve(subsidiaryFunctionPtr, subsidiaryArgs);
      info.solverInfos.push_back(solverInfo);
      if (not solverInfo.success) return info;
      this->notify(ControlMessages::SOLUTION_CHANGED);
      this->notify(ControlMessages::STEP_ENDED);

      /// Calculate predictor for a particular step
      for (int ls = 1; ls < loadSteps_; ++ls) {
        subsidiaryArgs.actStep = ls;
        this->notify(ControlMessages::STEP_STARTED, subsidiaryArgs.actStep, subsidiaryArgs.stepSize);

        pathFollowingType_.intermediatePrediction(nonOp, subsidiaryArgs);

        solverInfo = nonLinearSolver->solve(subsidiaryFunctionPtr, subsidiaryArgs);

        adaptiveStepSizing(solverInfo, nonOp, subsidiaryArgs, activationFunction);

        info.solverInfos.push_back(solverInfo);
        if (not solverInfo.success) return info;
        this->notify(ControlMessages::SOLUTION_CHANGED);
        this->notify(ControlMessages::STEP_ENDED);
      }

      auto stop     = std::chrono::high_resolution_clock::now();
      auto duration = duration_cast<std::chrono::milliseconds>(stop - start);

      this->notify(ControlMessages::CONTROL_ENDED, subsidiaryArgs.totalIte, duration.count(), pathFollowingType_.name);
      info.success = true;
      return info;
    }

  private:
    std::shared_ptr<NonLinearSolver> nonLinearSolver;
    int loadSteps_;
    double stepSize_;
    PathFollowingType pathFollowingType_;
    AdaptiveStepSizing adaptiveStepSizing;
    ActivationFunction activationFunction;
  };
}  // namespace Ikarus
