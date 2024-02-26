// SPDX-FileCopyrightText: 2021-2024 The Ikarus Developers mueller@ibb.uni-stuttgart.de
// SPDX-License-Identifier: LGPL-3.0-or-later

#include "controllogger.hh"

#include "spdlog/spdlog.h"
#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wswitch-enum"

namespace Ikarus {

void ControlLogger::updateImpl(ControlMessages message, const ControlState& info) {
  switch (message) {
    case ControlMessages::CONTROL_STARTED:
      start_ = std::chrono::high_resolution_clock::now();
      spdlog::info("=====================================================================");
      spdlog::info("Started {}", info.name);
      spdlog::info("=====================================================================");
      break;
    case ControlMessages::CONTROL_ENDED:
      stop_     = std::chrono::high_resolution_clock::now();
      duration_ = duration_cast<std::chrono::milliseconds>(stop_ - start_);
      spdlog::info("End of {}", info.name);
      spdlog::info("Total number of iterations: {:3d}", info.totalIterations);
      spdlog::info("Elapsed time: {} ms", duration_.count());
      break;
    case ControlMessages::STEP_STARTED:
      spdlog::info("Load step: {:>4} {:>39} {:<.2e}", info.currentStep, "Step size = ", info.stepSize);
      spdlog::info("---------------------------------------------------------------------");
      break;
    case ControlMessages::STEP_ENDED:
      spdlog::info("Load factor after completion of load step {:>2} {:>4} {:<.2e}", info.currentStep, " is ",
                   info.lambda);
      spdlog::info("=====================================================================");
      break;
    default:
      break; // default: do nothing when notified
  }
}
} // namespace Ikarus
#pragma GCC diagnostic pop
