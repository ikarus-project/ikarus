//
// Created by lex on 14/12/2021.
//

#pragma once
#include <string>

#include "spdlog/spdlog.h"

#include "ikarus/utils/Observer/observer.h"
#include "ikarus/utils/Observer/observerMessages.h"

class NonLinearSolverLogger : public IObserver<NonLinearSolverMessages> {
public:
  void updateImpl(NonLinearSolverMessages message) override {
    switch (message) {
      case NonLinearSolverMessages::ITERATION_STARTED:
        break;
      case NonLinearSolverMessages::INIT:
        iters = 0;
        dNorm = 0.0;
        rNorm = 0.0;
        spdlog::info("Non-linear solver started:");
        spdlog::info("  i ResidualNorm CorrectionNorm");
        break;
      case NonLinearSolverMessages::ITERATION_ENDED:
        spdlog::info("{:>3d} {:>12.2e} {:14.2e}", iters, rNorm, dNorm);
        ++iters;
        break;
      default:
        break;
    }
  }

  void updateImpl(NonLinearSolverMessages message, double val) override {
    switch (message) {
      case NonLinearSolverMessages::RESIDUALNORM_UPDATED:
        rNorm = val;
        break;
      case NonLinearSolverMessages::CORRECTIONNORM_UPDATED:
        dNorm = val;
        break;
      default:
        break;
    }
  }

  void updateImpl(NonLinearSolverMessages message, int intVal, double val1, double val2) override {
    switch (message) {
      case NonLinearSolverMessages::FINISHED_SUCESSFULLY:
        spdlog::info(
            "Nonlinear solver finished: {} iterations and final residualNorm {:03.2e} < {:01.2e} <-- desired "
            "tolerance.",
            intVal, val1, val2);
        break;
      default:
        break;
    }
  }

  void updateImpl(NonLinearSolverMessages, const Eigen::VectorXd&) override {}

private:
  int iters{0};
  double dNorm{0};
  double rNorm{0};
};
