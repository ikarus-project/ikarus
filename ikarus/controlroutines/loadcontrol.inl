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

  /// Dummy execution of the code with loadParameter as 0 in order to save information of the undeformed (or initial)
  /// configuration while using, for example, the class ControlSubsamplingVertexVTKWriter
  NonLinearSolverState solverState{.success = true, .iterations = 0, .currentIter = 0};
  updateAndNotifyControlState(controlState, nonOp, solverState);
  controlState.initialConfig = false;

  for (int ls = 0; ls < loadSteps_; ++ls) {
    this->notify(ControlMessages::STEP_STARTED, controlState);
    controlState.currentStep = ls;
    loadParameter += stepSize_;
    solverState = nonLinearSolver_->solve();
    updateAndNotifyControlState(controlState, nonOp, solverState);
    if (not solverState.success)
      return controlState;
  }
  controlState.success = true;
  this->notify(ControlMessages::CONTROL_ENDED, controlState);
  return controlState;
}
} // namespace Ikarus
