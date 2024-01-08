// SPDX-FileCopyrightText: 2021-2024 The Ikarus Developers mueller@ibb.uni-stuttgart.de
// SPDX-License-Identifier: LGPL-3.0-or-later

/**
 * \file genericobserver.hh
 * \brief Observer implementation for calling a callback function
 */

#pragma once
#include <concepts>
#include <string>

#include <spdlog/spdlog.h>

#include <ikarus/utils/observer/observer.hh>
#include <ikarus/utils/observer/observermessages.hh>

namespace Ikarus {

  /**
   * @brief GenericObserver class for observing specific messages.
   *
   * This class template implements an observer for a specific message type.
   *
   * @tparam Messages The type of messages to be observed.
   */
  template <typename Messages>
  class GenericObserver : public IObserver<Messages> {
  public:
    /**
     * @brief Constructor for GenericObserver.
     *
     * Initializes the observer with a specific message and a function to be executed upon observation.
     *
     * @tparam F Type of the function to be executed.
     * @param p_message The message to be observed.
     * @param p_f The function to be executed with the current step .
     */
    template <typename F>
    GenericObserver(Messages p_message, F&& p_f) : message{p_message}, f{p_f} {}
    void updateImpl(Messages p_message) override {
      if (p_message == message) {
        f(step);
        ++step;
      }
    }

    Messages message;
    std::function<void(int)> f;
    int step{0};
  };
}  // namespace Ikarus
