// SPDX-FileCopyrightText: 2021-2024 The Ikarus Developers mueller@ibb.uni-stuttgart.de
// SPDX-License-Identifier: LGPL-3.0-or-later

/**
 * \file genericobserver.hh
 * \brief Observer implementation for calling a callback function
 */

#pragma once

#include <ikarus/utils/observer/observer.hh>

namespace Ikarus {

/**
 * \brief GenericObserver class for observing specific messages.
 *
 * This class template implements an observer for a specific message type.
 *
 * \tparam OBS The type of observable.
 */
template <typename OBS>
class GenericObserver : public IObserver<OBS>
{
public:
  using MessageType = typename OBS::MessageType;
  using StateType   = typename OBS::StateType;

  /**
   * \brief Constructor for GenericObserver.
   *
   * Initializes the observer with a specific message and a function to be executed upon observation.
   *
   * \tparam F Type of the function to be executed.
   * \param message The message to be observed.
   * \param f The function to be executed with the current step .
   */
  template <typename F>
  GenericObserver(MessageType message, F&& f)
      : message_{message},
        f_{f} {}
  void updateImpl(MessageType message, const StateType&) override {
    if (message_ == message) {
      f_(step_);
      ++step_;
    }
  }

private:
  MessageType message_;
  std::function<void(int)> f_;
  int step_{0};
};
} // namespace Ikarus
