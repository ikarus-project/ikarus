// SPDX-FileCopyrightText: 2021-2024 The Ikarus Developers mueller@ibb.uni-stuttgart.de
// SPDX-License-Identifier: LGPL-3.0-or-later

/*!
 * \file loadcontrol.inl
 * \brief Implementation of the run function
 * @ingroup  controlroutines
 */

#pragma once

namespace Ikarus {
template <typename NLS>
[[nodiscard(
    "The run method returns information of the control routine. You should store this information and check if "
    "it was successful")]] ControlState
LoadControl<NLS>::run() {
  ControlState controlState{.success = false, .currentStep = 0, .stepSize = stepSize_, .name = this->name()};
  auto& nonOp         = nonLinearSolver_->nonLinearOperator();
  controlState.sol    = &nonOp.firstParameter();
  controlState.lambda = nonOp.lastParameter();
  this->notify(ControlMessages::CONTROL_STARTED, controlState);
  auto& loadParameter = nonOp.lastParameter();

  for (int ls = 0; ls < loadSteps_; ++ls) {
    this->notify(ControlMessages::STEP_STARTED, controlState);
    controlState.currentStep = ls;
    loadParameter += stepSize_;
    auto solverState = nonLinearSolver_->solve();
    updateAndNotifyControlState(controlState, nonOp, solverState);
    if (not solverState.success)
      return controlState;
  }
  controlState.success = true;
  this->notify(ControlMessages::CONTROL_ENDED, controlState);
  return controlState;
}
} // namespace Ikarus
