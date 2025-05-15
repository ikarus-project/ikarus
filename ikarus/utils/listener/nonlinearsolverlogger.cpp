// SPDX-FileCopyrightText: 2021-2025 The Ikarus Developers ikarus@ibb.uni-stuttgart.de
// SPDX-License-Identifier: LGPL-3.0-or-later

#include "nonlinearsolverlogger.hh"

#include <spdlog/spdlog.h>

namespace Ikarus {

void NonLinearSolverLogger::init() {
  iters_ = 1;
  rNorm_ = 0.0;
  dNorm_ = 0.0;
  spdlog::info("Non-linear solver started:");
  spdlog::info("{:<11} {:<20} {:<20} {:<20}", "Ite", "normR", "normD", "lambda");
  spdlog::info("-------------------------------------------------------------------------------");
}

void NonLinearSolverLogger::iterationEnded() {
  if (lambda_.has_value())
    spdlog::info("{} {:<10d} {:<20.2e} {:<20.2e} {:<20.2e}", "", iters_, rNorm_, dNorm_, lambda_.value());
  else
    spdlog::info("{} {:<10d} {:<20.2e} {:<20.2e} {:^20}", "", iters_, rNorm_, dNorm_, " - ");
  ++iters_;
}

void NonLinearSolverLogger::finishedSuccessfully(int numberOfIterations) {
  spdlog::info("Number of iterations: {}", numberOfIterations);
}

} // namespace Ikarus
