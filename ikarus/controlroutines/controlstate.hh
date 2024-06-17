// SPDX-FileCopyrightText: 2021-2024 The Ikarus Developers mueller@ibb.uni-stuttgart.de
// SPDX-License-Identifier: LGPL-3.0-or-later

/**
 * \file controlstate.hh
 * \brief Defines the ControlState structure for storing results of a control routine.
 * \ingroup  controlroutines
 */

#pragma once

#include <vector>

#include <Eigen/Core>

#include <ikarus/solver/nonlinearsolver/solverstate.hh>

namespace Ikarus {

/**
 * \struct ControlState
 * \brief Structure containing information about the state of the control routine.
 *
 * \details This structure holds information about the success of a control routine,
 * including current load step number, load step size, value of the load factor, the solution vector,
 * state of the non-linear solver and the name of the control routine.
 *
 * \ingroup  controlroutines
 */
struct ControlState
{
  bool success{false};                                      ///< Flag indicating the success of the control.
  int currentStep{-1};                                      ///< Information about the current step number.
  double stepSize{std::numeric_limits<double>::infinity()}; ///< Information about the step size.
  double lambda{std::numeric_limits<double>::infinity()};   ///< Value of the load factor.
  std::vector<NonLinearSolverState>
      solverStates{};                ///< Vector containing information about the state of the nonlinear solver.
  std::string name{};                ///< Information about the name of the control method.
  const Eigen::VectorX<double>* sol; ///< A pointer to the solution vector, for example, the displacement vector.
};

} // namespace Ikarus
