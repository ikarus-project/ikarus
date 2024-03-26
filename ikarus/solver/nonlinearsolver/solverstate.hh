// SPDX-FileCopyrightText: 2021-2024 The Ikarus Developers mueller@ibb.uni-stuttgart.de
// SPDX-License-Identifier: LGPL-3.0-or-later

/**
 * \file solverstate.hh
 * \brief Implementation of the Newton-Raphson method for solving nonlinear equations.
 */

#pragma once

namespace Ikarus {
/**
 * \brief Information about the state of a non-linear solver.
 *
 * This structure holds information about the success of a non-linear solver,
 * including the residual norm, correction norm, and the number of iterations performed.
 */
struct NonLinearSolverState
{
  /**
   * \brief Conversion operator to check the success of the solver.
   *
   * \return `true` if the solver was successful, `false` otherwise.
   */
  explicit operator bool() const { return success; }
  bool success{false};                                            ///< Flag indicating the success of the solver.
  double residualNorm{std::numeric_limits<double>::infinity()};   ///< Value of the residual norm.
  double correctionNorm{std::numeric_limits<double>::infinity()}; ///< Value of the correction norm.
  int iterations{-1};  ///< Total number of iterations performed by the non-linear solver.
  int currentIter{-1}; ///< Current iteration number.
};
} // namespace Ikarus
