// SPDX-FileCopyrightText: 2021-2025 The Ikarus Developers ikarus@ibb.uni-stuttgart.de
// SPDX-License-Identifier: LGPL-3.0-or-later

/**
 * \file controllogger.hh
 * \brief Observer implementation for logging control routines
 */

#pragma once
#include <chrono>

#include <ikarus/utils/broadcaster/broadcastermessages.hh>
#include <ikarus/utils/listener/listener.hh>

namespace Ikarus {
/**
 * \brief ControlLogger class for logging control messages.
 *
 * This class implements an observer for control messages and logs relevant information based on the received
 * messages.
 */
class ControlLogger : public Listener
{
public:
  template <typename BC>
  ControlLogger& subscribeTo(BC& bc) {
    this->subscribe(bc, [&](ControlMessages message, const BC::State& state) { this->updateImpl(message, state); });

    return *this;
  }

  /**
   * \brief Implementation of the update method for control message logging.
   *
   * \param message The received control message.
   */
  void updateImpl(ControlMessages message);
  /**
   * \brief Implementation of the update method for logging control messages with string values.
   *
   * \param message The received control message.
   * \param val The string value associated with the message.
   */
  void updateImpl(ControlMessages message, const std::string& val);
  /**
   * \brief Implementation of the update method for logging control messages with an integer and a string value.
   *
   * \param message The received control message.
   * \param val1 The integer value associated with the message.
   * \param val2 The string value associated with the message.
   */
  void updateImpl(ControlMessages message, int val1, const std::string& val2);
  /**
   * \brief Implementation of the update method for logging control messages with an integer and a double value.
   *
   * \param message The received control message.
   * \param val1 The integer value associated with the message.
   * \param val2 The double value associated with the message.
   */
  void updateImpl(ControlMessages message, int val1, double val2);

  /**
   * \brief Implementation of the update method for logging control messages with a control routine state.
   *
   * \param message The received control message.
   * \param state The received control state.
   */
  void updateImpl(ControlMessages message, const Concepts::ControlRoutineState auto& state) {
    updateImpl(message);
    updateImpl(message, state.loadStep, state.stepSize);
    if (message == ControlMessages::CONTROL_STARTED)
      name_ = state.name;
  }

private:
  using TimePoint = std::chrono::time_point<std::chrono::high_resolution_clock>;
  TimePoint start_{};
  TimePoint stop_{};
  std::chrono::milliseconds duration_{};
  std::string name_;
};
} // namespace Ikarus
