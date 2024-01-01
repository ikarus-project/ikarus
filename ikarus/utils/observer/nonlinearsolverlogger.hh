// SPDX-FileCopyrightText: 2021-2024 The Ikarus Developers mueller@ibb.uni-stuttgart.de
// SPDX-License-Identifier: LGPL-3.0-or-later

#pragma once
#include "observer.hh"
#include "observermessages.hh"

namespace Ikarus {
  class NonLinearSolverLogger : public IObserver<NonLinearSolverMessages> {
  public:
    void updateImpl(NonLinearSolverMessages message) final;
    void updateImpl(NonLinearSolverMessages message, double val) final;
    void updateImpl(NonLinearSolverMessages message, int intVal) final;
    void updateImpl(NonLinearSolverMessages, const std::string&) final {}
    void updateImpl(NonLinearSolverMessages, int, const std::string&) final {}
    void updateImpl(NonLinearSolverMessages, int, double) final {}
    void updateImpl(NonLinearSolverMessages, const Eigen::VectorXd&) final {}

  private:
    int iters{0};
    double dNorm{0};
    double rNorm{0};
    double lambda{0};
  };
}  // namespace Ikarus