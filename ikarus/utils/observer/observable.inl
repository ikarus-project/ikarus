// SPDX-FileCopyrightText: 2021-2024 The Ikarus Developers mueller@ibb.uni-stuttgart.de
// SPDX-License-Identifier: LGPL-3.0-or-later

/**
 * \file observable.hh
 * \brief Implementation of the member functions of the class IObservable
 */

#pragma once
#include <map>

#include <Eigen/Core>

#include <ikarus/utils/concepts.hh>
#include <ikarus/utils/observer/observer.hh>

namespace Ikarus {

template <typename MT, typename ST>
requires(Concepts::Observable<MT, ST>)
void IObservable<MT, ST>::subscribe(MT message, std::shared_ptr<IObserver<IObservable>> observer) {
  observers_[message];
  auto&& vectorOfObserversOfASpecificMessage = observers_[message];
  vectorOfObserversOfASpecificMessage.push_back(observer);
}

template <typename MT, typename ST>
requires(Concepts::Observable<MT, ST>)
void IObservable<MT, ST>::subscribeAll(std::shared_ptr<IObserver<IObservable>> observer) {
  for (auto& msg : messages_)
    subscribe(msg, observer);
}

template <typename MT, typename ST>
requires(Concepts::Observable<MT, ST>)
void IObservable<MT, ST>::subscribeAll(std::initializer_list<std::shared_ptr<IObserver<IObservable>>> observers) {
  for (auto& observer : observers)
    for (auto& msg : messages_)
      subscribe(msg, observer);
}

template <typename MT, typename ST>
requires(Concepts::Observable<MT, ST>)
void IObservable<MT, ST>::unSubscribe(MT message, std::shared_ptr<IObserver<IObservable>> observer) {
  auto vectorOfObserversOfASpecificMessage = observers_[message];
  std::ranges::remove_if(vectorOfObserversOfASpecificMessage, [&observer](auto&& obs) { return obs == observer; });
}

template <typename MT, typename ST>
requires(Concepts::Observable<MT, ST>)
void IObservable<MT, ST>::unSubscribeAll(std::shared_ptr<IObserver<IObservable>> observer) {
  for (auto& msg : messages_)
    unSubscribe(msg, observer);
}

template <typename MT, typename ST>
requires(Concepts::Observable<MT, ST>)
void IObservable<MT, ST>::notify(MT message, const ST& state) {
  auto vectorOfObserversOfASpecificMessage = observers_[message];
  for (auto&& obs : vectorOfObserversOfASpecificMessage)
    obs->update(message, state);
}
} // namespace Ikarus
