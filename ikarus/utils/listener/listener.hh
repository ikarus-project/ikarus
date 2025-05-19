// SPDX-FileCopyrightText: 2021-2025 The Ikarus Developers ikarus@ibb.uni-stuttgart.de
// SPDX-License-Identifier: LGPL-3.0-or-later

/**
 * \file listener.hh
 * \brief Implementation of the observer design pattern with broadcasters
 * \ingroup observer
 */

#pragma once
#include <algorithm>
#include <cassert>
#include <functional>
#include <iostream>
#include <memory>
#include <vector>

namespace Ikarus {

/**
 * \brief Implements a listener
 * \details The functions that the listener is listening to are stored in a vector of shared_ptr<void>. This type
 * erasure has the advantage that we can listen to different function signatures and the Listener not being a template.
 * This works, because the deleter of the stored objects is not bound to the type information and thus there are no
 * memory leaks possible.
 */
class Listener
{
public:
  using Token = std::shared_ptr<void>;

  /**
   * \brief Function to subscribe to a broadcaster with a given function (either a lambda, std::function or function
   * pointer).
   *
   * \tparam Broadcaster the type of the Broadcaster (for example a NonlinearSolver or ControlRoutine)
   * \param broadcaster the broadcaster
   * \param callback the function
   */
  template <typename Broadcaster>
  auto subscribe(Broadcaster& broadcaster,
                 std::function<void(typename Broadcaster::MessageType, const typename Broadcaster::State&)> callback) {
    auto token = broadcaster.registerListener(std::move(callback));
    tokens.push_back(token);
    return token;
  }

  /**
   * \brief Unsubscribe from all listeners. At the moment unsubscribing can't be done more granularly.
   */
  void unSubscribeAll() {
    for (auto& token : tokens) {
      if (token) {
        token.reset();
      }
    }
    tokens.clear();
  }

  /**
   * \brief Unsubscribe from the last subscribed listener.
   */
  void unSubscribeLast() {
    if (!tokens.empty()) {
      tokens.back().reset();
      tokens.pop_back();
    }
  }

  /**
   * \brief Unsubscribe from a specific token.
   */
  void unSubscribe(const Token& token) {
    auto it = std::ranges::find(tokens, token);
    if (it != tokens.end()) {
      (*it).reset();
      tokens.erase(it);
    }
  }

private:
  std::vector<Token> tokens;
};

} // namespace Ikarus