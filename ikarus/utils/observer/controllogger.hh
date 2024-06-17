// SPDX-FileCopyrightText: 2021-2024 The Ikarus Developers mueller@ibb.uni-stuttgart.de
// SPDX-License-Identifier: LGPL-3.0-or-later

/**
 * \file controllogger.hh
 * \brief Observer implementation for logging control routines
 */

#pragma once
#include "observable.hh"
#include "observer.hh"
#include "observermessages.hh"

#include <chrono>

#include <ikarus/controlroutines/controlstate.hh>
namespace Ikarus {
/**
 * \brief ControlLogger class for logging control messages.
 *
 * This class implements an observer for control messages and logs relevant information based on the received
 * messages.
 */
class ControlLogger : public IObserver<ControlObservable>
{
public:
  /**
   * \brief Implementation of the update method for control message logging.
   *
   * \param message The received control message.
   * \param state The state of the control routine needed for logging.
   */
  void updateImpl(MessageType message, const StateType& state) final;

private:
  using TimePoint = std::chrono::time_point<std::chrono::high_resolution_clock>;
  TimePoint start_{};
  TimePoint stop_{};
  std::chrono::milliseconds duration_{};
};
} // namespace Ikarus
