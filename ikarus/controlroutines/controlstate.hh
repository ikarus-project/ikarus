// SPDX-FileCopyrightText: 2021-2024 The Ikarus Developers mueller@ibb.uni-stuttgart.de
// SPDX-License-Identifier: LGPL-3.0-or-later

/**
 * \file controlstate.hh
 * \brief Defines the ControlState structure for storing results of a control routine.
 * \ingroup  controlroutines
 */

#pragma once

#include <ikarus/solver/nonlinearsolver/solverstate.hh>

namespace Ikarus {

/**
 * \struct ControlState
 * \brief Structure containing information about the state of the control routine.
 *
 * \details This structure holds information about the success of a control routine,
 * including current load step number, total number of iterations performed by the non-linear solver,
 * load step size, value of the load factor, state of the non-linear solver and the name of the control routine.
 */
struct ControlState
{
  bool success{false};      ///< Flag indicating the success of the control.
  bool initialConfig{true}; ///< Flag indicating if the state of deformation is in the initial configuration
  int currentStep{-1};      ///< Information about the current step number.
  int totalIterations{-1};  ///< Information about the total iterations performed.
  double stepSize{std::numeric_limits<double>::infinity()}; ///< Information about the step size.
  double lambda{std::numeric_limits<double>::infinity()};   ///< Value of the load factor.
  std::vector<NonLinearSolverState>
      solverState{};   ///< Vector containing information about the state of the nonlinear solver.
  std::string name{};  ///< Information about the name of the control method.
  Eigen::VectorXd sol; ///< The solution vector, for example, the displacement vector.
};

} // namespace Ikarus
