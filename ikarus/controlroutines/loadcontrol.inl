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
ControlInformation LoadControl<NLS>::run() {
  ControlInformation info{.success = false, .totalIterations = 0};
  ControlLoggerInformation logInfo{.currentStep = 0, .totalIterations = 0, .stepSize = stepSize_, .name = this->name()};
  auto& nonOp = nonLinearSolver_->nonLinearOperator();
  this->notify(ControlMessages::CONTROL_STARTED, logInfo);
  auto& loadParameter = nonOp.lastParameter();

  loadParameter += stepSize_;
  this->notify(ControlMessages::STEP_STARTED, logInfo);
  auto solverInfo = nonLinearSolver_->solve();
  info.solverInfos.push_back(solverInfo);
  info.totalIterations += solverInfo.iterations;
  logInfo.totalIterations = info.totalIterations;
  if (not solverInfo.success)
    return info;
  logInfo.lambda = nonOp.lastParameter();
  this->notify(ControlMessages::SOLUTION_CHANGED, logInfo);
  this->notify(ControlMessages::STEP_ENDED, logInfo);

  for (int ls = 1; ls < loadSteps_; ++ls) {
    logInfo.currentStep = ls;
    this->notify(ControlMessages::STEP_STARTED, logInfo);
    loadParameter += stepSize_;
    solverInfo = nonLinearSolver_->solve();
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
