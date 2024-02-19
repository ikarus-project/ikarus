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
 * \struct ControlInformation
 * \brief Structure containing information about the control results.
 */
struct ControlInformation
{
  bool success{false};                                   ///< Flag indicating the success of the control.
  std::vector<NonLinearSolverInformation> solverInfos{}; ///< Vector containing information from nonlinear solvers.
  int totalIterations{-1};                               ///< Total number of iterations performed.
};

/**
 * \struct ControlLoggerInformation
 * \brief Structure containing information about the control loggers.
 *
 * \details This structure holds information about the value of the current step number, the step size,
 * the name of the control routine and the total number of iterations performed.
 */
struct ControlLoggerInformation
{
  int currentStep{-1};                                      ///< Information about the current step number.
  int totalIterations{-1};                                  ///< Information about the total iterations performed.
  double stepSize{std::numeric_limits<double>::infinity()}; ///< Information about the step size.
  double lambda{std::numeric_limits<double>::infinity()};   ///< Value of the load factor.
  std::string name{};                                       ///< Information about the name of the control method.
};

} // namespace Ikarus
