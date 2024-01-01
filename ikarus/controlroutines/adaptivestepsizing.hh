// SPDX-FileCopyrightText: 2021-2024 The Ikarus Developers mueller@ibb.uni-stuttgart.de
// SPDX-License-Identifier: LGPL-3.0-or-later

#pragma once

#include <ikarus/controlroutines/pathfollowingfunctions.hh>
#include <ikarus/solver/nonlinearsolver/solverinfos.hh>

namespace Ikarus::AdaptiveStepSizing {
  struct NoOp {
    template <typename NonLinearOperator>
    void operator()(const NonLinearSolverInformation&, Ikarus::SubsidiaryArgs&, const NonLinearOperator&) {}
    [[nodiscard]] int targetIterations() const { return 0; }
    void setTargetIterations(int) {}
  };

  struct IterationBased {
    template <typename NonLinearOperator>
    void operator()(const NonLinearSolverInformation& solverInfo, Ikarus::SubsidiaryArgs& subsidiaryArgs,
                    const NonLinearOperator&) {
      if (subsidiaryArgs.currentStep == 0) return;
      if (targetIterations_ == 0)
        DUNE_THROW(Dune::InvalidStateException,
                   "For IterationBased, targetIterations should not be equal to 0. Try calling "
                   "setTargetIterations(int) first.");
      const auto previousIterations = solverInfo.iterations;
      if (previousIterations > targetIterations_)
        subsidiaryArgs.stepSize
            = sqrt(static_cast<double>(targetIterations_) / previousIterations) * subsidiaryArgs.stepSize;
    }
    [[nodiscard]] int targetIterations() const { return targetIterations_; }
    void setTargetIterations(int targetIterations) { targetIterations_ = targetIterations; }

  private:
    int targetIterations_{0};
  };
}  // namespace Ikarus::AdaptiveStepSizing
