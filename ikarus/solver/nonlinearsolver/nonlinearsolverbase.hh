// SPDX-FileCopyrightText: 2021-2025 The Ikarus Developers ikarus@ibb.uni-stuttgart.de
// SPDX-License-Identifier: LGPL-3.0-or-later

/**
 * \file nonlinearsolverbase.hh
 * \brief Base for all nonlinear solvers
 */

#pragma once
#include <ikarus/solver/nonlinearsolver/nonlinearsolverstate.hh>
#include <ikarus/utils/broadcaster/broadcaster.hh>
#include <ikarus/utils/broadcaster/broadcastermessages.hh>

namespace Ikarus {

/**
 * \brief Base for all nonlinear solvers. Defines the message interface that can be broadcasted to listeners.
 *
 * \tparam F Type of the differentiable function to solve.
 * \tparam Args Additional message signatures can be broadcasted
 */
template <typename F, typename... Args>
struct NonlinearSolverBase : Broadcaster<NonLinearSolverMessages, NonlinearSolverStateType<F>>
{
  using State = NonlinearSolverStateType<F>;
};

} // namespace Ikarus
