// SPDX-FileCopyrightText: 2021-2024 The Ikarus Developers mueller@ibb.uni-stuttgart.de
// SPDX-License-Identifier: LGPL-3.0-or-later

/**
 * \file observer.hh
 * \brief Implementation of the observer design pattern
 */

#pragma once

#include <dune/common/exceptions.hh>

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
  using StateType   = typename OBS::StateType;
  /** * \brief Virtual destructor for the observer interface. */
  virtual ~IObserver() = default;

  /**
   * \brief Update method for receiving notifications with a message.
   * \param message The message to be received.
   * \param state The state needed for logging.
   */
  void update(MessageType message, const StateType& state) {
    checkMessageType(message);
    updateImpl(message, state);
  }

protected:
  virtual void updateImpl(MessageType message, const StateType&) {};

private:
  void checkMessageType(MessageType message) {
    if (MessageType::END == message)
      DUNE_THROW(Dune::InvalidStateException, "The END enum type should not be used");
    if (MessageType::BEGIN == message)
      DUNE_THROW(Dune::InvalidStateException, "The BEGIN enum type should not be used");
  }
};
} // namespace Ikarus
