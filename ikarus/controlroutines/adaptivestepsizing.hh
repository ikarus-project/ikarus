// SPDX-FileCopyrightText: 2021-2023 The Ikarus Developers mueller@ibb.uni-stuttgart.de
// SPDX-License-Identifier: LGPL-3.0-or-later

#pragma once

#include <ikarus/solver/nonlinearsolver/solverinfos.hh>
#include <ikarus/utils/pathfollowingfunctions.hh>

namespace Ikarus::AdaptiveStepSizing {
  template <typename ScalarType = double>
  struct NoOp {
  public:
    template <typename NonLinearOperator>
    void operator()(const NonLinearSolverInformation&, Ikarus::SubsidiaryArgs&, const NonLinearOperator&) {}
    [[nodiscard]] int targetIterations() const { return 0; }
    void setTargetIterations(int) {}
  };

  template <typename ScalarType = double>
  struct IterationBased {
  public:
    template <typename NonLinearOperator>
    void operator()(const NonLinearSolverInformation& solverInfo, Ikarus::SubsidiaryArgs& subsidiaryArgs,
                    const NonLinearOperator&) {
      if (not((subsidiaryArgs.actualStep > 0) and (targetIterations_ != 0)))
        DUNE_THROW(Dune::InvalidStateException,
                   "For IterationBased, targetIterations should not be equal to 0 and actualStep should be "
                   "greater than 0. Try calling setTargetIterations(int) first.");
      ScalarType previousIterations = solverInfo.iterations;
      if (previousIterations > targetIterations_)
        subsidiaryArgs.stepSize
            = sqrt(static_cast<ScalarType>(targetIterations_) / previousIterations) * subsidiaryArgs.stepSize;
    }
    [[nodiscard]] int targetIterations() const { return targetIterations_; }
    void setTargetIterations(int targetIterations) { targetIterations_ = targetIterations; }

  private:
    int targetIterations_{0};
  };
}  // namespace Ikarus::AdaptiveStepSizing
