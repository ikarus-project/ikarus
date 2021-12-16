//
// Created by lex on 14/12/2021.
//

#pragma once

template <typename MessageType>
class IObserver {
public:
  virtual ~IObserver() = default;
  void update(MessageType message) {
    assert(MessageType::END != message && "The END enum type should not be used");
    updateImpl(message);
  };
  void update(MessageType message, double val) {
    assert(MessageType::END != message && "The END enum type should not be used");
    updateImpl(message, val);
  };

protected:
  virtual void updateImpl([[maybe_unused]] MessageType message){};
  virtual void updateImpl([[maybe_unused]] MessageType message, [[maybe_unused]] double val){};
};

template <typename MessageType>
MessageType& operator++(MessageType& e) {
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
    for (MessageType msg = MessageType::BEGIN; msg != MessageType::END; ++msg)
      messages_.push_back(msg);
  }
  virtual ~IObservable() = default;
  void subscribe(MessageType message, std::shared_ptr<IObserver<MessageType>> observer);
  void subscribeAll(std::shared_ptr<IObserver<MessageType>> observer);
  void unSubscribe(MessageType message, std::shared_ptr<IObserver<MessageType>> observer);
  void unSubscribeAll(std::shared_ptr<IObserver<MessageType>> observer);
  void notify(MessageType message);
  template <std::floating_point ScalarType>
  void notify(MessageType message, ScalarType val);

private:
  //  virtual subscribeImpl(MessageType message, std::shared_ptr<IObserver<MessageType>> observer);
  //  void unSubscribeImpl(MessageType message, std::shared_ptr<IObserver<MessageType>> observer);
  //  void notifyImpl(MessageType message);
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