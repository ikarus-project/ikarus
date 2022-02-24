//
// Created by Alex on 02.07.2021.
//

#pragma once

namespace Ikarus::Solving {
  template <typename AsssemblerType>
  concept HasAssembleMatrix = requires(AsssemblerType as, MatrixAffordances matA) {
    as.assembleMatrix(matA);
  };

  template <typename AsssemblerType>
  concept HasAssembleVector = requires(AsssemblerType as, Ikarus::FiniteElements::VectorAffordances vecA) {
    as.assembleVector(vecA);
  };

  template <typename AsssemblerType>
  concept HasAssembleScalar = requires(AsssemblerType as, ScalarAffordances scalA) {
    as.assembleScalar(scalA);
  };

}  // namespace Ikarus::Solving