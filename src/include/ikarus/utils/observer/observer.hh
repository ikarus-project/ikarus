/*
 * This file is part of the Ikarus distribution (https://github.com/IkarusRepo/Ikarus).
 * Copyright (c) 2022. The Ikarus developers.
 *
 * This library is free software; you can redistribute it and/or
 *  modify it under the terms of the GNU Lesser General Public
 *  License as published by the Free Software Foundation; either
 *  version 2.1 of the License, or (at your option) any later version.
 *
 * This library is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
 * Lesser General Public License for more details.
 *
 * You should have received a copy of the GNU Lesser General Public
 * License along with this library; if not, write to the Free Software
 * Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301 USA
 */



#pragma once
#include <map>

#include <Eigen/Core>

template <typename MessageType>
class IObserver {
public:
  virtual ~IObserver() = default;
  void update(MessageType message) {
    assert(MessageType::END != message && "The END enum type should not be used");
    assert(MessageType::BEGIN != message && "The BEGIN enum type should not be used");
    updateImpl(message);
  };
  void update(MessageType message, double val) {
    assert(MessageType::END != message && "The END enum type should not be used");
    assert(MessageType::BEGIN != message && "The BEGIN enum type should not be used");
    updateImpl(message, val);
  };

  void update(MessageType message, int intVal, double val1, double val2) {
    assert(MessageType::END != message && "The END enum type should not be used");
    assert(MessageType::BEGIN != message && "The BEGIN enum type should not be used");
    updateImpl(message, intVal, val1, val2);
  };

  void update(MessageType message, const Eigen::VectorXd& vec) {
    assert(MessageType::END != message && "The END enum type should not be used");
    assert(MessageType::BEGIN != message && "The BEGIN enum type should not be used");
    updateImpl(message, vec);
  };

protected:
  virtual void updateImpl([[maybe_unused]] MessageType message){};
  virtual void updateImpl([[maybe_unused]] MessageType message, [[maybe_unused]] double val){};
  virtual void updateImpl([[maybe_unused]] MessageType message, [[maybe_unused]] int intVal,
                          [[maybe_unused]] double val1, [[maybe_unused]] double val2){};
  virtual void updateImpl([[maybe_unused]] MessageType message, [[maybe_unused]] const Eigen::VectorXd& vec){};
};

template <typename MessageType>
MessageType& increment(MessageType& e) {
  if (e == MessageType::END) {
    throw std::out_of_range("for MessageType& operator ++ (MessageType&)");
  }
  e = MessageType(static_cast<typename std::underlying_type<MessageType>::type>(e) + 1);
  return e;
}

template <typename MessageType>
class IObservable {
public:
  IObservable() {
    for (MessageType msg = MessageType::BEGIN; msg != MessageType::END; increment(msg))
      messages_.push_back(msg);
  }
  virtual ~IObservable() = default;
  void subscribe(MessageType message, std::shared_ptr<IObserver<MessageType>> observer);
  void subscribeAll(std::shared_ptr<IObserver<MessageType>> observer);
  void subscribeAll(std::initializer_list<std::shared_ptr<IObserver<MessageType>>> observers);
  void unSubscribe(MessageType message, std::shared_ptr<IObserver<MessageType>> observer);
  void unSubscribeAll(std::shared_ptr<IObserver<MessageType>> observer);
  void notify(MessageType message);
  template <std::floating_point ScalarType>
  void notify(MessageType message, ScalarType val);
  void notify(MessageType message, int intVal, double val1, double val2);
  template <std::floating_point ScalarType>
  void notify(MessageType message, Eigen::VectorX<ScalarType> val);

private:
  using ObserverVector = std::vector<std::shared_ptr<IObserver<MessageType>>>;
  using ObserverMap    = std::map<MessageType, ObserverVector>;
  ObserverMap observers_;
  std::vector<MessageType> messages_;
};

template <typename MessageType>
void IObservable<MessageType>::subscribe(MessageType message, std::shared_ptr<IObserver<MessageType>> observer) {
  observers_[message];
  auto&& vectorOfObserversOfASpecificMessage = observers_[message];
  vectorOfObserversOfASpecificMessage.push_back(observer);
}

template <typename MessageType>
void IObservable<MessageType>::subscribeAll(std::shared_ptr<IObserver<MessageType>> observer) {
  for (auto& msg : messages_)
    subscribe(msg, observer);
}

template <typename MessageType>
void IObservable<MessageType>::subscribeAll(std::initializer_list<std::shared_ptr<IObserver<MessageType>>> observers) {
  for (auto& observer : observers)
    for (auto& msg : messages_)
      subscribe(msg, observer);
}

template <typename MessageType>
void IObservable<MessageType>::unSubscribe(MessageType message, std::shared_ptr<IObserver<MessageType>> observer) {
  auto vectorOfObserversOfASpecificMessage = observers_[message];
  std::ranges::remove_if(vectorOfObserversOfASpecificMessage, [&observer](auto&& obs) { return obs == observer; });
}

template <typename MessageType>
void IObservable<MessageType>::unSubscribeAll(std::shared_ptr<IObserver<MessageType>> observer) {
  for (auto& msg : messages_)
    unSubscribe(msg, observer);
}

template <typename MessageType>
void IObservable<MessageType>::notify(MessageType message) {
  auto vectorOfObserversOfASpecificMessage = observers_[message];
  for (auto&& obs : vectorOfObserversOfASpecificMessage)
    obs->update(message);
}

template <typename MessageType>
template <std::floating_point ScalarType>
void IObservable<MessageType>::notify(MessageType message, ScalarType val) {
  auto vectorOfObserversOfASpecificMessage = observers_[message];
  for (auto&& obs : vectorOfObserversOfASpecificMessage)
    obs->update(message, val);
}

template <typename MessageType>
void IObservable<MessageType>::notify(MessageType message, int intVal, double val1, double val2) {
  auto vectorOfObserversOfASpecificMessage = observers_[message];
  for (auto&& obs : vectorOfObserversOfASpecificMessage)
    obs->update(message, intVal, val1, val2);
}

template <typename MessageType>
template <std::floating_point ScalarType>
void IObservable<MessageType>::notify(MessageType message, Eigen::VectorX<ScalarType> vec) {
  auto vectorOfObserversOfASpecificMessage = observers_[message];
  for (auto&& obs : vectorOfObserversOfASpecificMessage)
    obs->update(message, vec);
}