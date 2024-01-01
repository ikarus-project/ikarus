// SPDX-FileCopyrightText: 2021-2024 The Ikarus Developers mueller@ibb.uni-stuttgart.de
// SPDX-License-Identifier: LGPL-3.0-or-later

#pragma once
#include <map>
#include <memory>

#include <Eigen/Core>

#include <ikarus/utils/makeenum.hh>

template <typename MessageType>
class IObserver {
public:
  virtual ~IObserver() = default;

  void update(MessageType message) {
    checkMessageType(message);
    updateImpl(message);
  }

  void update(MessageType message, double val) {
    checkMessageType(message);
    updateImpl(message, val);
  }

  void update(MessageType message, int val) {
    checkMessageType(message);
    updateImpl(message, val);
  };

  void update(MessageType message, const std::string& val) {
    checkMessageType(message);
    updateImpl(message, val);
  };

  void update(MessageType message, int val1, double val2) {
    checkMessageType(message);
    updateImpl(message, val1, val2);
  };

  void update(MessageType message, int val1, const std::string& val2) {
    checkMessageType(message);
    updateImpl(message, val1, val2);
  };

  void update(MessageType message, const Eigen::VectorXd& vec) {
    checkMessageType(message);
    updateImpl(message, vec);
  }

protected:
  virtual void updateImpl([[maybe_unused]] MessageType message){};
  virtual void updateImpl([[maybe_unused]] MessageType message, [[maybe_unused]] double val){};
  virtual void updateImpl([[maybe_unused]] MessageType message, [[maybe_unused]] int val){};
  virtual void updateImpl([[maybe_unused]] MessageType message, [[maybe_unused]] const std::string& val){};
  virtual void updateImpl([[maybe_unused]] MessageType message, [[maybe_unused]] int val1, const std::string& val2){};
  virtual void updateImpl([[maybe_unused]] MessageType message, [[maybe_unused]] int val1, double val2){};
  virtual void updateImpl([[maybe_unused]] MessageType message, [[maybe_unused]] const Eigen::VectorXd& vec) {}

private:
  void checkMessageType(MessageType message) {
    if (MessageType::END == message) DUNE_THROW(Dune::InvalidStateException, "The END enum type should not be used");
    if (MessageType::BEGIN == message)
      DUNE_THROW(Dune::InvalidStateException, "The BEGIN enum type should not be used");
  }
};

template <typename MessageType>
class IObservable {
public:
  IObservable() {
    for (MessageType msg = MessageType::BEGIN; msg != MessageType::END; Ikarus::increment(msg))
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
  void notify(MessageType message, int val);
  void notify(MessageType message, const std::string& val);
  void notify(MessageType message, int val1, double val2);
  void notify(MessageType message, int val1, const std::string& val2);
  template <std::floating_point ScalarType>
  void notify(MessageType message, Eigen::VectorX<ScalarType> vec);

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
void IObservable<MessageType>::notify(MessageType message, int val) {
  auto vectorOfObserversOfASpecificMessage = observers_[message];
  for (auto&& obs : vectorOfObserversOfASpecificMessage)
    obs->update(message, val);
}

template <typename MessageType>
void IObservable<MessageType>::notify(MessageType message, const std::string& val) {
  auto vectorOfObserversOfASpecificMessage = observers_[message];
  for (auto&& obs : vectorOfObserversOfASpecificMessage)
    obs->update(message, val);
}

template <typename MessageType>
void IObservable<MessageType>::notify(MessageType message, int val1, double val2) {
  auto vectorOfObserversOfASpecificMessage = observers_[message];
  for (auto&& obs : vectorOfObserversOfASpecificMessage)
    obs->update(message, val1, val2);
}

template <typename MessageType>
void IObservable<MessageType>::notify(MessageType message, int val1, const std::string& val2) {
  auto vectorOfObserversOfASpecificMessage = observers_[message];
  for (auto&& obs : vectorOfObserversOfASpecificMessage)
    obs->update(message, val1, val2);
}

template <typename MessageType>
template <std::floating_point ScalarType>
void IObservable<MessageType>::notify(MessageType message, Eigen::VectorX<ScalarType> vec) {
  auto vectorOfObserversOfASpecificMessage = observers_[message];
  for (auto&& obs : vectorOfObserversOfASpecificMessage)
    obs->update(message, vec);
}
