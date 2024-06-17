// SPDX-FileCopyrightText: 2021-2024 The Ikarus Developers mueller@ibb.uni-stuttgart.de
// SPDX-License-Identifier: LGPL-3.0-or-later

#include "nonlinearsolverlogger.hh"

#include "spdlog/spdlog.h"
#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wswitch-enum"

namespace Ikarus {
void NonLinearSolverLogger::updateImpl(MessageType message, const StateType& state) {
  switch (message) {
    case NonLinearSolverMessages::INIT:
      spdlog::info("Non-linear solver started:");
      if (state.lambda and state.subsidiaryFunction)
        spdlog::info("{:>4} {:>25} {:>19} {:>23} {:>15}", "Ite", "FirstOrderOptimality", "CorrectionNorm",
                     "SubsidiaryFunction", "LoadFactor");
      else
        spdlog::info("{:>4} {:>25} {:>19}", "Ite", "FirstOrderOptimality", "CorrectionNorm");
      spdlog::info("------------------------------------------------------------------------------------------");
      break;
    case NonLinearSolverMessages::ITERATION_STARTED:
      break;
    case NonLinearSolverMessages::ITERATION_ENDED:
      if (state.lambda and state.subsidiaryFunction)
        spdlog::info("{} {:>3d} {:>25.3e} {:>19.3e} {:>23.3e} {:>15.3e}", "", state.currentIter, state.residualNorm,
                     state.correctionNorm, state.subsidiaryFunction.value(), state.lambda.value());
      else
        spdlog::info("{} {:>3d} {:>25.3e} {:>19.3e}", "", state.currentIter, state.residualNorm, state.correctionNorm);
      break;
    case NonLinearSolverMessages::SOLVER_FINISHED:
      if (state.success)
        spdlog::info("Number of iterations by the non-linear solver: {}", state.iterations);
      else
        spdlog::warn("The non-linear solver FAILED to converge.");
      break;
    default:
      break;
  }
}
} // namespace Ikarus

#pragma GCC diagnostic pop
