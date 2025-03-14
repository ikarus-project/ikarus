// SPDX-FileCopyrightText: 2021-2025 The Ikarus Developers mueller@ibb.uni-stuttgart.de
// SPDX-License-Identifier: LGPL-3.0-or-later

/*!
 * \file loadcontrol.inl
 * \brief Implementation of the run function
 * \ingroup  controlroutines
 */

#pragma once

#include <type_traits>

#include "ikarus/utils/defaultfunctions.hh"

namespace Ikarus {
template <typename NLS>
ControlInformation LoadControl<NLS>::run(typename NLS::Domain& x) {
  ControlInformation info({false});
  this->notify(ControlMessages::CONTROL_STARTED, static_cast<std::string>(this->name()));
  auto& loadParameter = x.parameter();

  loadParameter = 0.0;
  this->notify(ControlMessages::STEP_STARTED, 0, stepSize_);
  auto solverInfo = nonLinearSolver_->solve(x);
  info.solverInfos.push_back(solverInfo);
  info.totalIterations += solverInfo.iterations;
  if (not solverInfo.success)
    return info;
  this->notify(ControlMessages::SOLUTION_CHANGED);
  this->notify(ControlMessages::STEP_ENDED);

  this->initialPrediction(x);

  for (int ls = 0; ls < loadSteps_; ++ls) {
    this->notify(ControlMessages::STEP_STARTED, ls, stepSize_);
    loadParameter += stepSize_;
    solverInfo = nonLinearSolver_->solve(x);
    info.solverInfos.push_back(solverInfo);
    info.totalIterations += solverInfo.iterations;
    if (not solverInfo.success)
      return info;
    this->notify(ControlMessages::SOLUTION_CHANGED);
    this->notify(ControlMessages::STEP_ENDED);
  }
  this->notify(ControlMessages::CONTROL_ENDED, info.totalIterations, static_cast<std::string>(this->name()));
  info.success = true;
  return info;
}

template <typename NLS>
void LoadControl<NLS>::initialPrediction(typename NLS::Domain& x) const {
  SolverTypeTag solverTag;

  auto&& residual = nonLinearSolver_->residual();

  if constexpr (traits::isSpecializationTypeAndNonTypes<Eigen::Matrix, typename NLS::JacobianType>::value)
    solverTag = SolverTypeTag::d_LDLT;
  else
    solverTag = SolverTypeTag::sd_SimplicialLDLT;

  auto&& R = residual(x);
  auto&& K = derivative(residual)(x);

  auto y = x;
  y.parameter() += stepSize_; // artificially set the lambda from the next step
  nonLinearSolver_->updateFunction()(
      y, utils::zeroIncrementTag); // get the inhomogeneous part of the solution of the next step
  R -= K * (y.globalSolution() -
            x.globalSolution()); // compute the internal forces due to the displcament increment
                                 // F_Int_dir = K* delta_u_dir, delta_u_dir is only non-zero for the inhomogeneous part

  auto linearSolver = LinearSolver(solverTag); // solve for new displacements
  linearSolver.analyzePattern(K);
  linearSolver.factorize(K);
  Eigen::VectorXd dPredictor;
  linearSolver.solve(dPredictor, -R);

  x += dPredictor;
}
} // namespace Ikarus
