// SPDX-FileCopyrightText: 2021-2025 The Ikarus Developers ikarus@ibb.uni-stuttgart.de
// SPDX-License-Identifier: LGPL-3.0-or-later

/**
 * \file controlroutinebase.hh
 * \brief Base for all control routines
 */

#pragma once
#include <ikarus/controlroutines/controlroutinestate.hh>
#include <ikarus/utils/broadcaster/broadcaster.hh>
#include <ikarus/utils/broadcaster/broadcastermessages.hh>

namespace Ikarus {

/**
 * \brief Base for all control routines. Defines the message interface that can be broadcasted to listeners.
 *
 * \tparam F Type of the differentiable function to solve.
 * \tparam Args Additional custom message signatures, that can be broadcasted
 */
template <typename F, typename S = ControlRoutineStateType<F>, typename... Args>
struct ControlRoutineBase
    : public Broadcasters<void(ControlMessages), void(ControlMessages, const std::string&),
                          void(ControlMessages, int, const std::string&), void(ControlMessages, int, double),
                          void(ControlMessages, const S&), Args...>
{
  using State = S;
};

} // namespace Ikarus
