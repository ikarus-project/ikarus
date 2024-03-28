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
      spdlog::info("==================================================================================");
      spdlog::info("Started {}", state.name);
      break;
    case ControlMessages::CONTROL_ENDED:
      stop_     = std::chrono::high_resolution_clock::now();
      duration_ = duration_cast<std::chrono::milliseconds>(stop_ - start_);
      spdlog::info("End of {}", state.name);
      spdlog::info("Total number of iterations: {:3d}", state.totalIterations);
      spdlog::info("Elapsed time: {} ms", duration_.count());
      break;
    case ControlMessages::STEP_STARTED:
      spdlog::info("Load step: {:>4} {:>56} {:<.3e}", state.currentStep, "Step size = ", state.stepSize);
      spdlog::info("----------------------------------------------------------------------------------");
      break;
    case ControlMessages::STEP_ENDED:
      if (not state.initialConfig) {
        spdlog::info("Load factor: {:<.3e}", state.lambda);
        spdlog::info("==================================================================================");
      } else {
      }
      break;
    default:
      break; // default: do nothing when notified
  }
}
} // namespace Ikarus
#pragma GCC diagnostic pop
