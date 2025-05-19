// SPDX-FileCopyrightText: 2021-2025 The Ikarus Developers ikarus@ibb.uni-stuttgart.de
// SPDX-License-Identifier: LGPL-3.0-or-later

#include "controllogger.hh"

#include <spdlog/spdlog.h>

namespace Ikarus {

void ControlLogger::stepEnded() {
  spdlog::info("===============================================================================");
}

void ControlLogger::controlStarted(const std::string& name) {
  start_ = std::chrono::high_resolution_clock::now();
  spdlog::info("===============================================================================");
  spdlog::info("Started " + name);
  spdlog::info("===============================================================================");
}

void ControlLogger::stepStarted(int stepNumber, double stepSize) {
  spdlog::info("Load step: {:>4} {:>49} {:<.2e}", stepNumber, "Step size = ", stepSize);
  spdlog::info("-------------------------------------------------------------------------------");
}

void ControlLogger::controlEnded(int totalIterations, const std::string& name) {
  stop_     = std::chrono::high_resolution_clock::now();
  duration_ = duration_cast<std::chrono::milliseconds>(stop_ - start_);
  spdlog::info("End of " + name);
  spdlog::info("Total number of iterations: {:3d}", totalIterations);
  spdlog::info("Elapsed time: {} ms", duration_.count());
}
} // namespace Ikarus
