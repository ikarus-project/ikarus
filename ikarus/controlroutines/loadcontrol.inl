// SPDX-FileCopyrightText: 2021-2025 The Ikarus Developers ikarus@ibb.uni-stuttgart.de
// SPDX-License-Identifier: LGPL-3.0-or-later

/*!
 * \file loadcontrol.inl
 * \brief Implementation of the run function
 * \ingroup  controlroutines
 */

#pragma once

#include <ikarus/controlroutines/common.hh>
#include <ikarus/controlroutines/loadcontrol.hh>
#include <ikarus/utils/defaultfunctions.hh>

namespace Ikarus {
template <typename NLS>
[[nodiscard(
    "The run method returns information of the control routine. You should store this information and check if "
    "it was successful")]] ControlInformation
LoadControl<NLS>::run(typename NLS::Domain& x) {
  using enum ControlMessages;
  decltype(auto) nonOp = nonLinearSolver_->residual();

  ControlInformation info(this->name());
  auto state = typename LoadControl::State(x, info);
  this->notify(CONTROL_STARTED, state);

  // Initial step to check if the undeformed (or initial) state is in equilibrium
  this->notify(ControlMessages::STEP_STARTED, state);
  auto solverInfo = nonLinearSolver_->solve(x);
  updateAndNotifyControlInfo(info, solverInfo, state);
  if (not solverInfo.success)
    return info;

  state.stepSize = stepSize_;

  for (int ls = 0; ls < loadSteps_; ++ls) {
    state.loadStep = ls;
    this->notify(STEP_STARTED, state);
    this->predictor(x);
    solverInfo = nonLinearSolver_->solve(x, stepSize_);
    if (solverInfo.iterations == 0)
      x.syncParameterAndGlobalSolution(nonLinearSolver_->updateFunction());
    updateAndNotifyControlInfo(info, solverInfo, state);
    if (not solverInfo.success)
      return info;
  }
  this->notify(CONTROL_ENDED, state);
  info.success = true;
  return info;
}

template <typename NLS>
void LoadControl<NLS>::predictor(typename NLS::Domain& x) const {
  x.parameter() += stepSize_;
}
} // namespace Ikarus
