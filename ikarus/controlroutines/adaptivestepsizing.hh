// SPDX-FileCopyrightText: 2021-2023 The Ikarus Developers mueller@ibb.uni-stuttgart.de
// SPDX-License-Identifier: LGPL-3.0-or-later

#pragma once

#include <ikarus/controlroutines/pathfollowingfunctions.hh>
#include <ikarus/solver/nonlinearsolver/solverinfos.hh>

namespace Ikarus::AdaptiveStepSizing {
  struct NoOp {
  public:
    template <typename NonLinearOperator>
    void operator()(const NonLinearSolverInformation&, Ikarus::SubsidiaryArgs&, const NonLinearOperator&) {}
    [[nodiscard]] int targetIterations() const { return 0; }
    void setTargetIterations(const int) {}
  };

  struct IterationBased {
  public:
    template <typename NonLinearOperator>
    void operator()(const NonLinearSolverInformation& solverInfo, Ikarus::SubsidiaryArgs& subsidiaryArgs,
                    const NonLinearOperator&) {
      if (subsidiaryArgs.actualStep == 0) return;
      if (targetIterations_ == 0)
        DUNE_THROW(Dune::InvalidStateException,
                   "For IterationBased, targetIterations should not be equal to 0 and actualStep should be "
                   "greater than 0. Try calling setTargetIterations(int) first.");
      const auto previousIterations = solverInfo.iterations;
      if (previousIterations > targetIterations_)
        subsidiaryArgs.stepSize = sqrt(static_cast<double>(targetIterations_) / static_cast<double>(previousIterations))
                                  * subsidiaryArgs.stepSize;
    }
    [[nodiscard]] int targetIterations() const { return targetIterations_; }
    void setTargetIterations(const int targetIterations) { targetIterations_ = targetIterations; }

  private:
    int targetIterations_{0};
  };
}  // namespace Ikarus::AdaptiveStepSizing
