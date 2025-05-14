// SPDX-FileCopyrightText: 2021-2025 The Ikarus Developers ikarus@ibb.uni-stuttgart.de
// SPDX-License-Identifier: LGPL-3.0-or-later

/**
 * \file controllogger.hh
 * \brief Listener implementation for logging control routines
 * \ingroup observer
 */

#pragma once
#include <chrono>

#include <ikarus/utils/broadcaster/broadcastermessages.hh>
#include <ikarus/utils/concepts.hh>
#include <ikarus/utils/listener/listener.hh>

namespace Ikarus {

/**
 * \brief Implementation of an observer for logging control routines.
 */
class ControlLogger : public Listener
{
public:
  template <typename BC>
  ControlLogger& subscribeTo(BC& bc) {
    this->subscribe(bc, [&](ControlMessages message, const BC::State& state) { this->update(message, state); });
    return *this;
  }

  /**
   * \brief Implementation of the update method for logging control messages with a control routine state.
   *
   * \param message The received control message.
   * \param state The received control state.
   */
  void update(ControlMessages message, const Concepts::ControlRoutineState auto& state) {
    switch (message) {
      case ControlMessages::STEP_ENDED:
        stepEnded();
        break;
      case ControlMessages::CONTROL_STARTED:
        controlStarted(state.information.name);
        break;
      case ControlMessages::STEP_STARTED:
        stepStarted(state.loadStep, state.stepSize);
        break;
      case ControlMessages::CONTROL_ENDED:
        controlEnded(state.information.totalIterations, state.information.name);
        break;
      default:
        break; // default: do nothing when notified
    }
  }

private:
  using TimePoint = std::chrono::time_point<std::chrono::high_resolution_clock>;
  TimePoint start_{};
  TimePoint stop_{};
  std::chrono::milliseconds duration_{};

  void stepEnded();
  void controlStarted(const std::string& name);
  void controlEnded(int totalIterations, const std::string& name);
  void stepStarted(int stepNumber, double stepSize);
};
} // namespace Ikarus
