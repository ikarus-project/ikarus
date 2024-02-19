// SPDX-FileCopyrightText: 2021-2024 The Ikarus Developers mueller@ibb.uni-stuttgart.de
// SPDX-License-Identifier: LGPL-3.0-or-later

/**
 * \file solverinfos.hh
 * \brief Implementation of the Newton-Raphson method for solving nonlinear equations.
 */

#pragma once

namespace Ikarus {
/**
 * \brief Information about the result of a non-linear solver.
 *
 * This structure holds information about the success of a non-linear solver,
 * including the residual norm, correction norm, and the number of iterations performed.
 */
struct NonLinearSolverInformation
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
  int iterations{-1}; ///< Total number of iterations performed by the non-linear solver.
};

/**
 * \brief Information about the loggers of a non-linear solver.
 *
 * \details This structure holds information about the value of the load factor, the residual norm,
 * the correction norm and the number of iterations performed.
 */
struct NonLinearSolverLoggingInformation
{
  int currentIter{-1}; ///< Current iteration number.
  int iterations{-1};  ///< Total number of iterations performed by the non-linear solver.
  double residualNorm{std::numeric_limits<double>::infinity()};   ///< Value of the residual norm.
  double correctionNorm{std::numeric_limits<double>::infinity()}; ///< Value of the correction norm.
};
} // namespace Ikarus
