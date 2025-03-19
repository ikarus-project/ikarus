// SPDX-FileCopyrightText: 2021-2025 The Ikarus Developers ikarus@ibb.uni-stuttgart.de
// SPDX-License-Identifier: LGPL-3.0-or-later

/*!
 * \file loadcontrol.inl
 * \brief Implementation of the run function
 * \ingroup  controlroutines
 */

#pragma once

#include <type_traits>

#include <ikarus/controlroutines/common.hh>
#include <ikarus/utils/defaultfunctions.hh>

namespace Ikarus {
template <typename NLS>
ControlInformation LoadControl<NLS>::run(typename NLS::Domain& x) {
  using enum ControlMessages;
  ControlInformation info({false});
  this->notify(CONTROL_STARTED, static_cast<std::string>(this->name()));
  auto& loadParameter = x.parameter();

  loadParameter = 0.0;
  this->notify(ControlMessages::STEP_STARTED, 0, stepSize_);
  auto solverInfo = nonLinearSolver_->solve(x);
  info.solverInfos.push_back(solverInfo);
  info.totalIterations += solverInfo.iterations;
  if (not solverInfo.success)
    return info;

  auto state = typename LoadControl::State{.domain = x};

  this->notify(SOLUTION_CHANGED, state);
  this->notify(STEP_ENDED, state);

  state.stepSize = stepSize_;

  this->initialPrediction(x);

  for (int ls = 0; ls < loadSteps_; ++ls) {
    this->notify(STEP_STARTED, ls, stepSize_);
    loadParameter += stepSize_;
    solverInfo = nonLinearSolver_->solve(x);
    info.solverInfos.push_back(solverInfo);
    info.totalIterations += solverInfo.iterations;
    if (not solverInfo.success)
      return info;

    state.loadStep = ls;
    this->notify(SOLUTION_CHANGED, state);
    this->notify(STEP_ENDED, state);
  }
  this->notify(CONTROL_ENDED, info.totalIterations, static_cast<std::string>(this->name()));
  info.success = true;
  return info;
}

template <typename NLS>
void LoadControl<NLS>::initialPrediction(typename NLS::Domain& x) const {
  auto y = x;
  y.parameter() += stepSize_;
  x += predictorForNewLoadLevel(*nonLinearSolver_, x, y);
}
} // namespace Ikarus
