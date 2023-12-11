// SPDX-FileCopyrightText: 2021-2023 The Ikarus Developers mueller@ibb.uni-stuttgart.de
// SPDX-License-Identifier: LGPL-3.0-or-later

#pragma once
#include "observer.hh"
#include "observermessages.hh"

namespace Ikarus {
  class NonLinearSolverLogger : public IObserver<NonLinearSolverMessages> {
  public:
    void updateImpl(NonLinearSolverMessages message) final;
    void updateImpl(NonLinearSolverMessages message, double val) final;
    void updateImpl(NonLinearSolverMessages message, int intVal, double val1, double val2) final;
    void updateImpl(NonLinearSolverMessages, const Eigen::VectorXd&) final {}

  private:
    int iters{0};
    double dNorm{0};
    double rNorm{0};
    double sNorm{0};
  };
}  // namespace Ikarus
