// SPDX-FileCopyrightText: 2021-2024 The Ikarus Developers mueller@ibb.uni-stuttgart.de
// SPDX-License-Identifier: LGPL-3.0-or-later

/*!
 * \file pathfollowing.inl
 * \brief Implementation of the run function
 * @ingroup  controlroutines
 */

#pragma once

#include <Eigen/Core>

#include <ikarus/controlroutines/adaptivestepsizing.hh>
#include <ikarus/controlroutines/controlinfos.hh>
#include <ikarus/controlroutines/pathfollowingfunctions.hh>
#include <ikarus/solver/nonlinearsolver/newtonraphsonwithscalarsubsidiaryfunction.hh>
#include <ikarus/utils/nonlinearoperator.hh>
#include <ikarus/utils/observer/observer.hh>
#include <ikarus/utils/observer/observermessages.hh>

namespace Ikarus {

template <typename NLS, typename PF, typename ASS>
requires(Impl::checkPathFollowingTemplates<NLS, PF, ASS>())
ControlInformation PathFollowing<NLS, PF, ASS>::run() {
  ControlInformation info{.success = false, .totalIterations = 0};
  ControlLoggerInformation logInfo{.currentStep = 0, .totalIterations = 0, .stepSize = stepSize_, .name = this->name()};
  auto& nonOp = nonLinearSolver_->nonLinearOperator();
  this->notify(ControlMessages::CONTROL_STARTED, logInfo);

  SubsidiaryArgs subsidiaryArgs;

  subsidiaryArgs.stepSize = stepSize_;
  subsidiaryArgs.DD.resizeLike(nonOp.firstParameter());
  subsidiaryArgs.DD.setZero();

  /// Initializing solver
  this->notify(ControlMessages::STEP_STARTED, logInfo);
  auto subsidiaryFunctionPtr = [&](auto&& args) { return pathFollowingType_(args); };
  pathFollowingType_.initialPrediction(nonOp, subsidiaryArgs);
  auto solverInfo = nonLinearSolver_->solve(pathFollowingType_, subsidiaryArgs);
  info.solverInfos.push_back(solverInfo);
  info.totalIterations += solverInfo.iterations;
  logInfo.totalIterations = info.totalIterations;
  if (not solverInfo.success)
    return info;
  logInfo.lambda = nonOp.lastParameter();
  this->notify(ControlMessages::SOLUTION_CHANGED, logInfo);
  this->notify(ControlMessages::STEP_ENDED, logInfo);

  /// Calculate predictor for a particular step
  for (int ls = 1; ls < steps_; ++ls) {
    subsidiaryArgs.currentStep = ls;

    adaptiveStepSizing_(solverInfo, subsidiaryArgs, nonOp);

    logInfo.currentStep = subsidiaryArgs.currentStep;
    logInfo.stepSize    = subsidiaryArgs.stepSize;

    this->notify(ControlMessages::STEP_STARTED, logInfo);

    pathFollowingType_.intermediatePrediction(nonOp, subsidiaryArgs);

    solverInfo = nonLinearSolver_->solve(pathFollowingType_, subsidiaryArgs);

    info.solverInfos.push_back(solverInfo);
    info.totalIterations += solverInfo.iterations;
    logInfo.totalIterations = info.totalIterations;
    if (not solverInfo.success)
      return info;
    logInfo.lambda = nonOp.lastParameter();
    this->notify(ControlMessages::SOLUTION_CHANGED, logInfo);
    this->notify(ControlMessages::STEP_ENDED, logInfo);
  }

  this->notify(ControlMessages::CONTROL_ENDED, logInfo);
  info.success = true;
  return info;
}
} // namespace Ikarus
