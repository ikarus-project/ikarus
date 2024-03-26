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
#include <ikarus/controlroutines/controlstate.hh>
#include <ikarus/controlroutines/pathfollowingfunctions.hh>
#include <ikarus/solver/nonlinearsolver/newtonraphsonwithscalarsubsidiaryfunction.hh>
#include <ikarus/utils/nonlinearoperator.hh>
#include <ikarus/utils/observer/observer.hh>
#include <ikarus/utils/observer/observermessages.hh>

namespace Ikarus {

template <typename NLS, typename PF, typename ASS>
requires(Impl::checkPathFollowingTemplates<NLS, PF, ASS>())
ControlState PathFollowing<NLS, PF, ASS>::run() {
  ControlState controlState{
      .success = false, .currentStep = 0, .totalIterations = 0, .stepSize = stepSize_, .name = this->name()};
  auto& nonOp = nonLinearSolver_->nonLinearOperator();
  this->notify(ControlMessages::CONTROL_STARTED, controlState);

  SubsidiaryArgs subsidiaryArgs;

  subsidiaryArgs.stepSize = stepSize_;
  subsidiaryArgs.DD.resizeLike(nonOp.firstParameter());
  subsidiaryArgs.DD.setZero();
  subsidiaryArgs.Dlambda = 0.0;

  /// Initializing solver
  /// Dummy execution of the code inorder to save information of the undeformed (or initial)
  /// configuration while using, for example, the class ControlSubsamplingVertexVTKWriter
  this->notify(ControlMessages::STEP_STARTED, controlState);
  auto solverState = nonLinearSolver_->solve(pathFollowingType_, subsidiaryArgs);
  if (not solverState.success)
    return controlState;
  updateAndNotifyControlState(controlState, nonOp, solverState);

  for (int ls = 0; ls < steps_; ++ls) {
    subsidiaryArgs.currentStep = ls;
    controlState.currentStep   = subsidiaryArgs.currentStep;
    if (ls == 0) {
      this->notify(ControlMessages::STEP_STARTED, controlState);
      pathFollowingType_.initialPrediction(nonOp, subsidiaryArgs);
      solverState = nonLinearSolver_->solve(pathFollowingType_, subsidiaryArgs);
      if (not solverState.success)
        return controlState;
      updateAndNotifyControlState(controlState, nonOp, solverState);
    } else {
      adaptiveStepSizing_(solverState, subsidiaryArgs, nonOp);
      controlState.stepSize = subsidiaryArgs.stepSize;

      this->notify(ControlMessages::STEP_STARTED, controlState);
      pathFollowingType_.intermediatePrediction(nonOp, subsidiaryArgs);
      solverState = nonLinearSolver_->solve(pathFollowingType_, subsidiaryArgs);
      if (not solverState.success)
        return controlState;
      updateAndNotifyControlState(controlState, nonOp, solverState);
    }
  }

  this->notify(ControlMessages::CONTROL_ENDED, controlState);
  controlState.success = true;
  return controlState;
}
} // namespace Ikarus
