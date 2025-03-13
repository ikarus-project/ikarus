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

  SolverTypeTag solverTag;
  // initial prediction
  if constexpr (traits::isSpecializationTypeAndNonTypes<Eigen::Matrix, typename NLS::JacobianType>::value)
  solverTag = SolverTypeTag::d_LDLT;
else
  solverTag = SolverTypeTag::sd_SimplicialLDLT;

auto&& R = nonLinearSolver_->residual()(x);
auto&& K = derivative(nonLinearSolver_->residual())(x);
std::remove_cvref_t<typename NLS::CorrectionType> d_dir;
auto y = x;
nonLinearSolver_->updateFunction()(y,utils::zeroIncrementTag);
R-= K*y.globalSolution();

auto linearSolver = LinearSolver(solverTag); // for the linear predictor step
linearSolver.analyzePattern(K);
linearSolver.factorize(K);
Eigen::VectorXd dispPredictor;
linearSolver.solve(dispPredictor, -R);
nonLinearSolver_->updateFunction()(x,dispPredictor);

  // end initial prediction

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
} // namespace Ikarus
