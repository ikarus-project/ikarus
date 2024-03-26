// SPDX-FileCopyrightText: 2021-2024 The Ikarus Developers mueller@ibb.uni-stuttgart.de
// SPDX-License-Identifier: LGPL-3.0-or-later

/**
 * \file nonlinearsolverlogger.hh
 * \brief Observer implementation for logging non-linear solvers
 */

#pragma once
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
class NonLinearSolverLogger : public IObserver<IObservable<NonLinearSolverMessages, NonLinearSolverState>>
{
public:
  /**
   * \brief Handles the update when a NonLinearSolverMessages with a NonLinearSolverLoggingInformation is received.
   *
   * \param message The NonLinearSolverMessages received.
   * \param info The non-linear solver information needed for logging.
   */
  void updateImpl(NonLinearSolverMessages message, const NonLinearSolverState& info) final;
};
} // namespace Ikarus
