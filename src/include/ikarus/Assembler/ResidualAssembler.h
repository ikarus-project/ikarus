//
// Created by Alex on 26.06.2021.
//

#pragma once
#include <ikarus/FiniteElements/FiniteElementPolicies.h>

namespace Ikarus::Assembler {
  template <typename DofManagerType>
  class VectorAssembler {
    explicit VectorAssembler(DofManagerType& dofManager) : dofManager_{dofManager} {}

    auto getVector(Ikarus::FiniteElements::ElementVectorAffordances& vectorAffordances) {
      for (const auto&& [fe, dofs, vars] : dofManager_->elementDofsVariableTuple) {
      }
      return vec;
    }

  private:
    DofManagerType* dofManager_;
    Eigen::VectorXd vec;
  };
}  // namespace Ikarus::Assembler