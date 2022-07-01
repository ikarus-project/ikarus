/*
 * This file is part of the Ikarus distribution (https://github.com/IkarusRepo/Ikarus).
 * Copyright (c) 2022. The Ikarus developers.
 *
 * This library is free software; you can redistribute it and/or
 *  modify it under the terms of the GNU Lesser General Public
 *  License as published by the Free Software Foundation; either
 *  version 2.1 of the License, or (at your option) any later version.
 *
 * This library is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
 * Lesser General Public License for more details.
 *
 * You should have received a copy of the GNU Lesser General Public
 * License along with this library; if not, write to the Free Software
 * Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301 USA
 */

#pragma once
#include <string>

#include "spdlog/spdlog.h"

#include <ikarus/utils/observer/observer.hh>
#include <ikarus/utils/observer/observerMessages.hh>

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
