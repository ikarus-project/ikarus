// SPDX-FileCopyrightText: 2022 The Ikarus Developers mueller@ibb.uni-stuttgart.de
// SPDX-License-Identifier: LGPL-3.0-or-later

#pragma once
#include "observer.hh"
#include "observerMessages.hh"

#include <string>

#include "spdlog/spdlog.h"

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
        sNorm = 0.0;
        spdlog::info("Non-linear solver started:");
        spdlog::info("  i ResidualNorm CorrectionNorm ScalarSubsidiaryNorm");
        break;
      case NonLinearSolverMessages::ITERATION_ENDED:
        spdlog::info("{:>3d} {:>12.2e} {:14.2e} {:20.2e}", iters, rNorm, dNorm, sNorm);
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
      case NonLinearSolverMessages::SCALARSUBSIDIARY_UPDATED:
        sNorm = val;
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
  double sNorm{0};
};
