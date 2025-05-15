// SPDX-FileCopyrightText: 2021-2025 The Ikarus Developers ikarus@ibb.uni-stuttgart.de
// SPDX-License-Identifier: LGPL-3.0-or-later

/*!
 * \file loadcontrol.inl
 * \brief Implementation of the run function
 * \ingroup  controlroutines
 */

#pragma once

#include <ikarus/controlroutines/loadcontrol.hh>

namespace Ikarus {
template <typename NLS>
ControlInformation LoadControl<NLS>::run(typename NLS::Domain& x) {
  using enum ControlMessages;
  decltype(auto) nonOp = nonLinearSolver_->residual();

  ControlInformation info(this->name());
  auto state = typename LoadControl::State(x, info);
  this->notify(CONTROL_STARTED, state);

  auto& loadParameter = x.parameter();

  loadParameter = 0.0;
  this->notify(ControlMessages::STEP_STARTED, state);
  auto solverInfo = nonLinearSolver_->solve(x);
  info.solverInfos.push_back(solverInfo);
  info.totalIterations += solverInfo.iterations;
  if (not solverInfo.success)
    return info;

  this->notify(SOLUTION_CHANGED, state);
  this->notify(STEP_ENDED, state);

  state.stepSize = stepSize_;

  for (int ls = 0; ls < loadSteps_; ++ls) {
    state.loadStep = ls;
    this->notify(STEP_STARTED, state);
    loadParameter += stepSize_;
    solverInfo = nonLinearSolver_->solve(x);
    info.solverInfos.push_back(solverInfo);
    info.totalIterations += solverInfo.iterations;
    if (not solverInfo.success)
      return info;

    this->notify(SOLUTION_CHANGED, state);
    this->notify(STEP_ENDED, state);
  }
  this->notify(CONTROL_ENDED, state);
  info.success = true;
  return info;
}
} // namespace Ikarus
