// SPDX-FileCopyrightText: 2021-2025 The Ikarus Developers ikarus@ibb.uni-stuttgart.de
// SPDX-License-Identifier: LGPL-3.0-or-later

/**
 * \file nonlinearsolverlogger.hh
 * \brief Observer implementation for logging non-linear solvers
 */

#pragma once
#include <ikarus/utils/broadcaster/broadcastermessages.hh>
#include <ikarus/utils/listener/listener.hh>

namespace Ikarus {
/**
 * \brief Implementation of an observer for logging non-linear solvers.
 * \ingroup observer
 * This class inherits from the IObserver class and provides specific
 * implementations for updating based on NonLinearSolverMessages.
 */
class NonLinearSolverLogger : public Listener
{
public:
  template <typename BC>
  NonLinearSolverLogger& subscribeTo(BC& bc) {
    this->subscribe(bc,
                    [&](NonLinearSolverMessages message, const BC::State& state) { this->updateImpl(message, state); });

    return *this;
  }

  /**
   * \brief Handles the update when a NonLinearSolverMessages is received.
   *
   * \param message The NonLinearSolverMessages received.
   */
  void updateImpl(NonLinearSolverMessages message);

  /**
   * \brief Handles the update when a NonLinearSolverMessages with a double value is received.
   *
   * \param message The NonLinearSolverMessages received.
   * \param val The double value associated with the message.
   */
  void updateImpl(NonLinearSolverMessages message, double val);

  /**
   * \brief Handles the update when a NonLinearSolverMessages with an integer value is received.
   *
   * \param message The NonLinearSolverMessages received.
   * \param intVal The integer value associated with numberOfIterations.
   */
  void updateImpl(NonLinearSolverMessages message, int intVal);

  template <typename State>
  void updateImpl(NonLinearSolverMessages message, const State& state) {
    updateImpl(message);
    updateImpl(message, state.information.residualNorm);
    updateImpl(message, state.information.correctionNorm);
    updateImpl(message, state.information.iterations);
  }

private:
  int iters_{0};
  double dNorm_{0};
  double rNorm_{0};
  double lambda_{0};
};
} // namespace Ikarus
