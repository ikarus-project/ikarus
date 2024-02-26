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
ControlState LoadControl<NLS>::run() {
  ControlState controlState{
      .success = false, .currentStep = 0, .totalIterations = 0, .stepSize = stepSize_, .name = this->name()};
  auto& nonOp = nonLinearSolver_->nonLinearOperator();
  this->notify(ControlMessages::CONTROL_STARTED, controlState);
  auto& loadParameter = nonOp.lastParameter();

  loadParameter += stepSize_;
  this->notify(ControlMessages::STEP_STARTED, controlState);
  auto solverInfo = nonLinearSolver_->solve();
  controlState.solverInfos.push_back(solverInfo);
  controlState.totalIterations += solverInfo.iterations;
  if (not solverInfo.success)
    return controlState;
  controlState.lambda = nonOp.lastParameter();
  this->notify(ControlMessages::SOLUTION_CHANGED, controlState);
  this->notify(ControlMessages::STEP_ENDED, controlState);

  for (int ls = 1; ls < loadSteps_; ++ls) {
    controlState.currentStep = ls;
    this->notify(ControlMessages::STEP_STARTED, controlState);
    loadParameter += stepSize_;
    solverInfo = nonLinearSolver_->solve();
    controlState.solverInfos.push_back(solverInfo);
    controlState.totalIterations += solverInfo.iterations;
    if (not solverInfo.success)
      return controlState;
    controlState.lambda = nonOp.lastParameter();
    this->notify(ControlMessages::SOLUTION_CHANGED, controlState);
    this->notify(ControlMessages::STEP_ENDED, controlState);
  }
  this->notify(ControlMessages::CONTROL_ENDED, controlState);
  controlState.success = true;
  return controlState;
}
} // namespace Ikarus
