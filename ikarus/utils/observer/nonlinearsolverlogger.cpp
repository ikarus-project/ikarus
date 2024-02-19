// SPDX-FileCopyrightText: 2021-2024 The Ikarus Developers mueller@ibb.uni-stuttgart.de
// SPDX-License-Identifier: LGPL-3.0-or-later

#include "nonlinearsolverlogger.hh"

#include "spdlog/spdlog.h"
#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wswitch-enum"

namespace Ikarus {

void NonLinearSolverLogger::updateImpl(NonLinearSolverMessages message, const NonLinearSolverLoggingInformation& info) {
  switch (message) {
    case NonLinearSolverMessages::INIT:
      spdlog::info("Non-linear solver started:");
      spdlog::info("{:<11} {:<20} {:<20}", "Ite", "normR", "normD");
      spdlog::info("---------------------------------------------------------------------");
      break;
    case NonLinearSolverMessages::ITERATION_STARTED:
      break;
    case NonLinearSolverMessages::ITERATION_ENDED:
      spdlog::info("{} {:<10d} {:<20.2e} {:<20.2e}", "", info.currentIter, info.residualNorm, info.correctionNorm);
      break;
    case NonLinearSolverMessages::RESIDUALNORM_UPDATED:
      break;
    case NonLinearSolverMessages::CORRECTIONNORM_UPDATED:
      break;
    case NonLinearSolverMessages::SOLUTION_CHANGED:
      break;
    case NonLinearSolverMessages::FINISHED_SUCESSFULLY:
      spdlog::info("Number of iterations by the non-linear solver: {}", info.iterations);
      break;
    default:
      break;
  }
}
} // namespace Ikarus

#pragma GCC diagnostic pop
