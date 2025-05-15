// SPDX-FileCopyrightText: 2021-2025 The Ikarus Developers ikarus@ibb.uni-stuttgart.de
// SPDX-License-Identifier: LGPL-3.0-or-later

/**
 * \file genericlistener.hh
 * \brief Listener implementation for calling a callback function
 * \ingroup observer
 */

#pragma once
#include <spdlog/spdlog.h>

#include <ikarus/utils/broadcaster/broadcastermessages.hh>
#include <ikarus/utils/listener/listener.hh>

namespace Ikarus {

/**
 * \brief GenericListener class for observing specific messages. This class template implements an listener for a
 * specific message type.
 *
 * \tparam M The type of messages to be listend to.
 */
template <typename BC>
class GenericListener : public Listener
{
  using Messages = BC::MessageType;
  using State    = BC::State;

public:
  /**
   * \brief Constructor for GenericListener.
   *
   * Initializes the listener with a specific message and a function to be executed upon listening.
   *
   * \tparam F Type of the function to be executed.
   * \param message The message to be listend to.
   * \param f The function to be executed with the current step.
   */
  template <typename F>
  GenericListener(BC& bc, Messages message, F&& f)
      : message_{message},
        f_{std::forward<F>(f)} {
    this->subscribe(bc, [&](Messages message, const BC::State& state) { this->updateImpl(message, state); });
  }

  void updateImpl(Messages message, const State& state) {
    if (message_ == message)
      f_(state);
  }

private:
  Messages message_;
  std::function<void(const State&)> f_;
};

template <typename BC, typename MT, typename F>
GenericListener(BC&, MT, F&&) -> GenericListener<BC>;

} // namespace Ikarus
