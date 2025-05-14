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
template <typename MT>
class GenericListener : public Listener
{
  using MessageType = MT;
  // using State    = BC::State;

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
  GenericListener(MessageType message)
      : message_(message) {}

  /**
   * \brief Registers a given function to the broadcaster
   *
   * \tparam BC the type of the broadcaster
   * \tparam F the type of the function
   * \param bc the broadcaster
   * \param f the function
   * \return GenericListener&
   */
  template <typename BC, typename F>
  GenericListener& subscribeTo(BC& bc, F&& f) {
    this->subscribe(bc, [&](BC::MessageType message, const BC::State& state) {
      if (message_ == message)
        f(state);
    });
    return *this;
  }

private:
  MessageType message_;
};

template <typename MT>
GenericListener(MT) -> GenericListener<MT>;

} // namespace Ikarus
