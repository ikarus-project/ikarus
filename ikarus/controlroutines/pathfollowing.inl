// SPDX-FileCopyrightText: 2021-2024 The Ikarus Developers mueller@ibb.uni-stuttgart.de
// SPDX-License-Identifier: LGPL-3.0-or-later

/*!
 * \file pathfollowing.inl
 * \brief Implementation of the run function
 * @ingroup  controlroutines
 */

#pragma once

#include <Eigen/Core>

#include <ikarus/controlroutines/controlstate.hh>
#include <ikarus/controlroutines/pathfollowingfunctions.hh>
#include <ikarus/solver/nonlinearsolver/newtonraphsonwithscalarsubsidiaryfunction.hh>
#include <ikarus/utils/observer/observermessages.hh>

namespace Ikarus {

template <typename NLS, typename PF, typename ASS>
requires(Impl::checkPathFollowingTemplates<NLS, PF, ASS>())
[[nodiscard(
    "The run method returns information of the control routine. You should store this information and check if "
    "it was successful")]] ControlState
PathFollowing<NLS, PF, ASS>::run() {
  ControlState controlState{.success = false, .currentStep = 0, .stepSize = stepSize_, .name = this->name()};
  auto& nonOp         = nonLinearSolver_->nonLinearOperator();
  controlState.sol    = &nonOp.firstParameter();
  controlState.lambda = nonOp.lastParameter();
  this->notify(ControlMessages::CONTROL_STARTED, controlState);

  SubsidiaryArgs subsidiaryArgs;

  subsidiaryArgs.stepSize = stepSize_;
  subsidiaryArgs.DD.resizeLike(nonOp.firstParameter());
  subsidiaryArgs.DD.setZero();
  subsidiaryArgs.Dlambda = 0.0;

  this->notify(ControlMessages::STEP_STARTED, controlState);
  pathFollowingType_.initialPrediction(nonOp, subsidiaryArgs);
  auto solverState = nonLinearSolver_->solve(pathFollowingType_, subsidiaryArgs);
  updateAndNotifyControlState(controlState, nonOp, solverState);
  if (not solverState.success)
    return controlState;

  for (int ls = 1; ls < steps_; ++ls) {
    subsidiaryArgs.currentStep = ls;
    controlState.currentStep   = subsidiaryArgs.currentStep;
    adaptiveStepSizing_(solverState, subsidiaryArgs, nonOp);
    controlState.stepSize = subsidiaryArgs.stepSize;

    this->notify(ControlMessages::STEP_STARTED, controlState);
    pathFollowingType_.intermediatePrediction(nonOp, subsidiaryArgs);
    solverState = nonLinearSolver_->solve(pathFollowingType_, subsidiaryArgs);
    updateAndNotifyControlState(controlState, nonOp, solverState);
    if (not solverState.success)
      return controlState;
  }

  controlState.success = true;
  this->notify(ControlMessages::CONTROL_ENDED, controlState);
  return controlState;
}
} // namespace Ikarus
