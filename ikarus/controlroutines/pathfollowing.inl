// SPDX-FileCopyrightText: 2021-2025 The Ikarus Developers ikarus@ibb.uni-stuttgart.de
// SPDX-License-Identifier: LGPL-3.0-or-later

/*!
 * \file pathfollowing.inl
 * \brief Implementation of the run function
 * \ingroup  controlroutines
 */

#pragma once

#include <Eigen/Core>

#include <ikarus/controlroutines/adaptivestepsizing.hh>
#include <ikarus/controlroutines/controlinfos.hh>
#include <ikarus/controlroutines/pathfollowing.hh>
#include <ikarus/solver/nonlinearsolver/newtonraphsonwithscalarsubsidiaryfunction.hh>
#include <ikarus/utils/broadcaster/broadcastermessages.hh>
#include <ikarus/utils/differentiablefunction.hh>

namespace Ikarus {

template <typename NLS, typename PF, typename ASS>
requires(Impl::checkPathFollowingTemplates<NLS, PF, ASS>())
[[nodiscard(
    "The run method returns information of the control routine. You should store this information and check if "
    "it was successful")]] ControlInformation
PathFollowing<NLS, PF, ASS>::run(typename NLS::Domain& req) {
  using enum ControlMessages;
  auto& residual = nonLinearSolver_->residual();

  ControlInformation info(this->name());
  info.totalIterations = 0;
  subsidiaryArgs_.setZero(req.globalSolution());
  auto state = typename PathFollowing::State(req, info, subsidiaryArgs_);

  this->notify(CONTROL_STARTED, state);
  subsidiaryArgs_.stepSize = stepSize_;

  // Initial step to check if the undeformed (or initial) state is in equilibrium
  state.loadStep = -1;
  this->notify(STEP_STARTED, state);
  auto solverInfo = nonLinearSolver_->solve(req, pathFollowingType_, subsidiaryArgs_);
  updateAndNotifyControlInfo(info, solverInfo, state);
  if (not solverInfo.success)
    return info;

  state.loadStep = 0;
  state.stepSize = stepSize_;
  this->notify(STEP_STARTED, state);
  pathFollowingType_.initialPrediction(req, *nonLinearSolver_, subsidiaryArgs_);
  solverInfo = nonLinearSolver_->solve(req, pathFollowingType_, subsidiaryArgs_);
  updateAndNotifyControlInfo(info, solverInfo, state);
  if (not solverInfo.success)
    return info;

  /// Calculate predictor for a particular step
  for (int ls = 1; ls < steps_; ++ls) {
    subsidiaryArgs_.currentStep = ls;
    state.loadStep              = ls;

    adaptiveStepSizing_(solverInfo, subsidiaryArgs_, residual);
    pathFollowingType_.intermediatePrediction(req, *nonLinearSolver_, subsidiaryArgs_);

    state.stepSize = subsidiaryArgs_.stepSize;
    this->notify(STEP_STARTED, state);

    solverInfo = nonLinearSolver_->solve(req, pathFollowingType_, subsidiaryArgs_);
    updateAndNotifyControlInfo(info, solverInfo, state);
    if (not solverInfo.success)
      return info;
  }

  this->notify(CONTROL_ENDED, state);
  info.success = true;
  return info;
}
} // namespace Ikarus
