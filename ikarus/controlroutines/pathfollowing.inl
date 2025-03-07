// SPDX-FileCopyrightText: 2021-2025 The Ikarus Developers mueller@ibb.uni-stuttgart.de
// SPDX-License-Identifier: LGPL-3.0-or-later

/*!
 * \file pathfollowing.inl
 * \brief Implementation of the run function
 * @ingroup  controlroutines
 */

#pragma once

#include <Eigen/Core>

#include <ikarus/controlroutines/adaptivestepsizing.hh>
#include <ikarus/controlroutines/controlinfos.hh>
#include <ikarus/controlroutines/pathfollowingfunctions.hh>
#include <ikarus/solver/nonlinearsolver/newtonraphsonwithscalarsubsidiaryfunction.hh>
#include <ikarus/utils/differentiablefunction.hh>
#include <ikarus/utils/observer/observer.hh>
#include <ikarus/utils/observer/observermessages.hh>

namespace Ikarus {

template <typename NLS, typename PF, typename ASS>
requires(Impl::checkPathFollowingTemplates<NLS, PF, ASS>())
ControlInformation PathFollowing<NLS, PF, ASS>::run(typename NLS::Domain& req) {
  ControlInformation info;
  auto& residual = nonLinearSolver_->residual();
  this->notify(ControlMessages::CONTROL_STARTED, pathFollowingType_.name());

  SubsidiaryArgs subsidiaryArgs;

  info.totalIterations    = 0;
  subsidiaryArgs.stepSize = stepSize_;
  subsidiaryArgs.DD.resizeLike(req.globalSolution());
  subsidiaryArgs.DD.setZero();

  /// Initializing solver
  this->notify(ControlMessages::STEP_STARTED, 0, subsidiaryArgs.stepSize);
  pathFollowingType_.initialPrediction(req, residual, subsidiaryArgs);
  auto solverInfo = nonLinearSolver_->solve(req, pathFollowingType_, subsidiaryArgs);
  info.solverInfos.push_back(solverInfo);
  info.totalIterations += solverInfo.iterations;
  if (not solverInfo.success)
    return info;
  this->notify(ControlMessages::SOLUTION_CHANGED);
  this->notify(ControlMessages::STEP_ENDED);

  /// Calculate predictor for a particular step
  for (int ls = 1; ls < steps_; ++ls) {
    subsidiaryArgs.currentStep = ls;

    adaptiveStepSizing_(solverInfo, subsidiaryArgs, residual);

    this->notify(ControlMessages::STEP_STARTED, subsidiaryArgs.currentStep, subsidiaryArgs.stepSize);

    pathFollowingType_.intermediatePrediction(req, residual, subsidiaryArgs);

    solverInfo = nonLinearSolver_->solve(req, pathFollowingType_, subsidiaryArgs);

    info.solverInfos.push_back(solverInfo);
    info.totalIterations += solverInfo.iterations;
    if (not solverInfo.success)
      return info;
    this->notify(ControlMessages::SOLUTION_CHANGED);
    this->notify(ControlMessages::STEP_ENDED);
  }

  this->notify(ControlMessages::CONTROL_ENDED, info.totalIterations, pathFollowingType_.name());
  info.success = true;
  return info;
}
} // namespace Ikarus
