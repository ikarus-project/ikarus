// SPDX-FileCopyrightText: 2021-2024 The Ikarus Developers mueller@ibb.uni-stuttgart.de
// SPDX-License-Identifier: LGPL-3.0-or-later

/**
 * \file observable.hh
 * \brief Implementation of the observable design pattern
 */

#pragma once
#include <map>
#include <memory>

#include <Eigen/Core>

#include <ikarus/controlroutines/controlstate.hh>
#include <ikarus/solver/nonlinearsolver/solverstate.hh>
#include <ikarus/utils/concepts.hh>
#include <ikarus/utils/makeenum.hh>
#include <ikarus/utils/observer/observer.hh>
#include <ikarus/utils/observer/observermessages.hh>

namespace Ikarus {

/**
 * \brief Generic observable interface for the Observer design pattern.
 * See \cite gamma1995design for a description of the design pattern
 * \tparam MT The type of message that the observable can handle.
 * \tparam ST The type of the state of the observable.
 * \ingroup observer
 */
template <typename MT, typename ST>
requires(Concepts::Observable<MT, ST>)
class IObservable
{
public:
  using MessageType = MT;
  using StateType   = ST;
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
   * \param state The state information needed for logging.
   */
  void notify(MessageType message, const StateType& state);

private:
  using ObserverVector = std::vector<std::shared_ptr<IObserver<IObservable>>>;
  using ObserverMap    = std::map<MessageType, ObserverVector>;
  ObserverMap observers_;
  std::vector<MessageType> messages_;
};

/**
 * \brief Alias for IObservable specific for control routines.
 */
using ControlObservable = IObservable<ControlMessages, ControlState>;

/**
 * \brief Alias for IObservable specific for non-linear solvers.
 */
using NonLinearSolverObservable = IObservable<NonLinearSolverMessages, NonLinearSolverState>;
} // namespace Ikarus

#include <ikarus/utils/observer/observable.inl>
