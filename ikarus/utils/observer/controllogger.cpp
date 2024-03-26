// SPDX-FileCopyrightText: 2021-2024 The Ikarus Developers mueller@ibb.uni-stuttgart.de
// SPDX-License-Identifier: LGPL-3.0-or-later

#include "controllogger.hh"

#include "spdlog/spdlog.h"
#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wswitch-enum"

namespace Ikarus {

void ControlLogger::updateImpl(ControlMessages message, const ControlState& state) {
  switch (message) {
    case ControlMessages::CONTROL_STARTED:
      start_ = std::chrono::high_resolution_clock::now();
      spdlog::info("=====================================================================");
      spdlog::info("Started {}", state.name);
      spdlog::info("=====================================================================");
      break;
    case ControlMessages::CONTROL_ENDED:
      stop_     = std::chrono::high_resolution_clock::now();
      duration_ = duration_cast<std::chrono::milliseconds>(stop_ - start_);
      spdlog::info("End of {}", state.name);
      spdlog::info("Total number of iterations: {:3d}", state.totalIterations);
      spdlog::info("Elapsed time: {} ms", duration_.count());
      break;
    case ControlMessages::STEP_STARTED:
      spdlog::info("Load step: {:>4} {:>39} {:<.2e}", state.currentStep, "Step size = ", state.stepSize);
      spdlog::info("---------------------------------------------------------------------");
      break;
    case ControlMessages::STEP_ENDED:
      spdlog::info("Load factor after completion of load step {:>2} {:>4} {:<.2e}", state.currentStep, " is ",
                   state.lambda);
      spdlog::info("=====================================================================");
      break;
    default:
      break; // default: do nothing when notified
  }
}
} // namespace Ikarus
#pragma GCC diagnostic pop
