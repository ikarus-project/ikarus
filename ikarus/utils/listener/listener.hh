#pragma once
#include <any>
#include <memory>
#include <vector>

#include <ikarus/utils/functionhelper.hh>
#include <ikarus/utils/traits.hh>

namespace Ikarus {

class Listener
{
public:
  // For the listener we simply use the token provided by the broadcaster.
  using Token = std::shared_ptr<void>;

  // Subscribe to a broadcaster with a callback.
  template <typename Broadcaster>
  void subscribe(Broadcaster& broadcaster,
                 std::function<void(typename Broadcaster::MessageType, const typename Broadcaster::State&)> callback) {
    auto token = broadcaster.registerListener(std::move(callback));
    tokens.push_back(token);
  }

  // Unsubscribe from all registered listeners.
  void unSubscribeAll() {
    for (auto& token : tokens) {
      if (token) {
        token.reset();
      }
    }
    tokens.clear();
  }

  // Unsubscribe from the last registered listener.
  void unSubscribeLast() {
    if (!tokens.empty()) {
      tokens.back().reset();
      tokens.pop_back();
    }
  }

  // Unsubscribe a specific token.
  void unSubscribe(const Token& token) {
    auto it = std::find(tokens.begin(), tokens.end(), token);
    if (it != tokens.end()) {
      assert(it->use_count() == 1 && "The token has external references");
      (*it).reset();
      tokens.erase(it);
    }
  }

private:
  std::vector<Token> tokens;
};

} // namespace Ikarus