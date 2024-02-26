// SPDX-FileCopyrightText: 2021-2024 The Ikarus Developers mueller@ibb.uni-stuttgart.de
// SPDX-License-Identifier: LGPL-3.0-or-later

/**
 * \file observer.hh
 * \brief Implementation of the observer design pattern
 */

#pragma once
#include <map>
#include <memory>

#include <Eigen/Core>

#include <ikarus/utils/makeenum.hh>
namespace Ikarus {

/**
 * \brief Generic observer interface for the Observer design pattern.
 * See \cite gamma1995design for a description of the design pattern
 * \ingroup observer
 * \tparam OBS The type of observable.
 */
template <typename OBS>
class IObserver
{
public:
  using MessageType = typename OBS::MessageType;
  using InfoType    = typename OBS::InfoType;
  /** * \brief Virtual destructor for the observer interface. */
  virtual ~IObserver() = default;

  /**
   * \brief Update method for receiving notifications with a message.
   * \param message The message to be received.
   */
  void update(MessageType message, const InfoType& info) {
    checkMessageType(message);
    updateImpl(message, info);
  }

protected:
  virtual void updateImpl(MessageType message, const InfoType&){};

private:
  void checkMessageType(MessageType message) {
    if (MessageType::END == message)
      DUNE_THROW(Dune::InvalidStateException, "The END enum type should not be used");
    if (MessageType::BEGIN == message)
      DUNE_THROW(Dune::InvalidStateException, "The BEGIN enum type should not be used");
  }
};

/**
 * \brief Generic observable interface for the Observer design pattern.
 * See \cite gamma1995design for a description of the design pattern
 * \tparam MT The type of message that the observable can handle.
 * \tparam IT The type of logging information.
 * \ingroup observer
 */
template <typename MT, typename IT>
class IObservable
{
public:
  using MessageType = MT;
  using InfoType    = IT;
  IObservable() {
    for (MessageType msg = MessageType::BEGIN; msg != MessageType::END; Ikarus::increment(msg))
      messages_.push_back(msg);
  }
  virtual ~IObservable() = default;
  /**
   * \brief Subscribe an observer to receive notifications for a specific message type.
   * \param message The message type to subscribe to.
   * \param observer The observer to be subscribed.
   */
  void subscribe(MessageType message, std::shared_ptr<IObserver<IObservable>> observer);
  /**
   * \brief Subscribe an observer to receive notifications for all message types.
   * \param observer The observer to be subscribed.
   */
  void subscribeAll(std::shared_ptr<IObserver<IObservable>> observer);
  /**
   * \brief Subscribe multiple observers to receive notifications for all message types.
   * \param observers List of observers to be subscribed.
   */
  void subscribeAll(std::initializer_list<std::shared_ptr<IObserver<IObservable>>> observers);
  /**
   * \brief Unsubscribe an observer from receiving notifications for a specific message type.
   * \param message The message type to unsubscribe from.
   * \param observer The observer to be unsubscribed.
   */
  void unSubscribe(MessageType message, std::shared_ptr<IObserver<IObservable>> observer);
  /**
   * \brief Unsubscribe an observer from receiving notifications for all message types.
   * \param observer The observer to be unsubscribed.
   */
  void unSubscribeAll(std::shared_ptr<IObserver<IObservable>> observer);
  /**
   * \brief Notify observers about a specific message type.
   * \param message The message type to notify about.
   * \param info The information to be logged.
   */
  void notify(MessageType message, const InfoType& info);

private:
  using ObserverVector = std::vector<std::shared_ptr<IObserver<IObservable>>>;
  using ObserverMap    = std::map<MessageType, ObserverVector>;
  ObserverMap observers_;
  std::vector<MessageType> messages_;
};

template <typename MT, typename IT>
void IObservable<MT, IT>::subscribe(MT message, std::shared_ptr<IObserver<IObservable>> observer) {
  observers_[message];
  auto&& vectorOfObserversOfASpecificMessage = observers_[message];
  vectorOfObserversOfASpecificMessage.push_back(observer);
}

template <typename MT, typename IT>
void IObservable<MT, IT>::subscribeAll(std::shared_ptr<IObserver<IObservable>> observer) {
  for (auto& msg : messages_)
    subscribe(msg, observer);
}

template <typename MT, typename IT>
void IObservable<MT, IT>::subscribeAll(std::initializer_list<std::shared_ptr<IObserver<IObservable>>> observers) {
  for (auto& observer : observers)
    for (auto& msg : messages_)
      subscribe(msg, observer);
}

template <typename MT, typename IT>
void IObservable<MT, IT>::unSubscribe(MT message, std::shared_ptr<IObserver<IObservable>> observer) {
  auto vectorOfObserversOfASpecificMessage = observers_[message];
  std::ranges::remove_if(vectorOfObserversOfASpecificMessage, [&observer](auto&& obs) { return obs == observer; });
}

template <typename MT, typename IT>
void IObservable<MT, IT>::unSubscribeAll(std::shared_ptr<IObserver<IObservable>> observer) {
  for (auto& msg : messages_)
    unSubscribe(msg, observer);
}

template <typename MT, typename IT>
void IObservable<MT, IT>::notify(MT message, const IT& info) {
  auto vectorOfObserversOfASpecificMessage = observers_[message];
  for (auto&& obs : vectorOfObserversOfASpecificMessage)
    obs->update(message, info);
}
} // namespace Ikarus
