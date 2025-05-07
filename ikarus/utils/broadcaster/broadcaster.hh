#pragma once
#include <algorithm>
#include <any>
#include <functional>
#include <memory>
#include <vector>

namespace Ikarus {

template <typename MT, typename S>
class Broadcaster
{
public:
using MessageType = MT;
using State = S;

using Callback = std::function<void(MT, const State&)>;
  using Token    = std::shared_ptr<Callback>;

  Token registerListener(Callback callback) {
    auto sp = std::make_shared<Callback>(std::move(callback));
    listeners.push_back(sp);
    return sp;
  }

  void unregisterListener(Token token) {
    if (token) {
      token.reset();
    }
  }

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