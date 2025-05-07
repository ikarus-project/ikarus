// SPDX-FileCopyrightText: 2021-2025 The Ikarus Developers ikarus@ibb.uni-stuttgart.de
// SPDX-License-Identifier: LGPL-3.0-or-later

/**
 * \file broadcaster.hh
 * \brief Implementation of the observer design pattern with broadcasters
 * \ingroup observer
 */

#pragma once
#include <algorithm>
#include <functional>
#include <memory>
#include <vector>

namespace Ikarus {

/**
 * \brief Implements a Broadcaster for a given MessageType and BroadcasterState
 *
 * \tparam MT the message type
 * \tparam S the broadcaster state
 */
template <typename MT, typename S>
class Broadcaster
{
public:
  using MessageType = MT;
  using State       = S;

  using Callback = std::function<void(MT, const State&)>;
  using Token    = std::shared_ptr<Callback>;

  /**
   * \brief This method is used to register a Listener function.
   * \details The function that is passed in is first stored in a shared_ptr. After this, the shared_ptr is added to the
   * vector of listener functions, which leads to a implicit conversion to a weak_ptr. The shared_ptr is then returned
   * to the Listener that has called this function to be stored in a vector of shared_ptr<void> \ref listener.hh.
   *
   * \param callback the callback function
   * \return Token
   */
  Token registerListener(Callback callback) {
    auto sp = std::make_shared<Callback>(std::move(callback));
    listeners.push_back(sp);
    return sp;
  }

  /**
   * \brief deregisters a specific function
   */
  void unregisterListener(Token token) {
    if (token) {
      token.reset();
    }
  }

  /**
   * \brief This calls all the registered functions.
   */
  void notify(MT message, const State& data) {
    trim();
    for (auto& weakCb : listeners) {
      if (auto callback = weakCb.lock()) {
        (*callback)(message, data);
      }
    }
  }

private:
  std::vector<std::weak_ptr<Callback>> listeners;

  void trim() {
    listeners.erase(
        std::remove_if(listeners.begin(), listeners.end(), [](const auto& weakCb) { return weakCb.expired(); }),
        listeners.end());
  }
};

} // namespace Ikarus