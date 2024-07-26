// SPDX-FileCopyrightText: 2021-2024 The Ikarus Developers mueller@ibb.uni-stuttgart.de
// SPDX-License-Identifier: LGPL-3.0-or-later

/**
 * \file nonlinearsolverlogger.hh
 * \brief Observer implementation for logging non-linear solvers
 */

#pragma once
#include "observable.hh"
#include "observer.hh"
#include "observermessages.hh"

#include <ikarus/solver/nonlinearsolver/solverstate.hh>

namespace Ikarus {
/**
 * \brief Implementation of an observer for logging non-linear solvers.
 * \ingroup observer
 * This class inherits from the IObserver class and provides specific
 * implementations for updating based on NonLinearSolverMessages.
 */
class NonLinearSolverLogger : public IObserver<NonLinearSolverObservable>
{
public:
  /**
   * \brief Handles the update when a NonLinearSolverMessages with a NonLinearSolverState is received.
   *
   * \param message The NonLinearSolverMessages received.
   * \param state The state of the non-linear solver needed for logging.
   */
  void updateImpl(MessageType message, const StateType& state) final;
};
} // namespace Ikarus
