// SPDX-FileCopyrightText: 2021-2024 The Ikarus Developers mueller@ibb.uni-stuttgart.de
// SPDX-License-Identifier: LGPL-3.0-or-later

/**
 * \file
 * \brief Defines the ControlInformation structure for storing control results.
 * \ingroup  controlroutines
 */

#pragma once

#include <ikarus/solver/nonlinearsolver/solverinfos.hh>

namespace Ikarus {

/**
 * \struct ControlState
 * \brief Structure containing information about the state of the control routine.
 */
struct ControlState
{
  bool success{false};                                      ///< Flag indicating the success of the control.
  int currentStep{-1};                                      ///< Information about the current step number.
  int totalIterations{-1};                                  ///< Information about the total iterations performed.
  double stepSize{std::numeric_limits<double>::infinity()}; ///< Information about the step size.
  double lambda{std::numeric_limits<double>::infinity()};   ///< Value of the load factor.
  std::vector<NonLinearSolverState> solverInfos{};          ///< Vector containing information from nonlinear solvers.
  std::string name{};                                       ///< Information about the name of the control method.
};

} // namespace Ikarus
