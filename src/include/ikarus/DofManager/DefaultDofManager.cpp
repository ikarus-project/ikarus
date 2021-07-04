//
// Created by ac126718 on 04.07.2021.
//

#include "DefaultDofManager.h"
namespace Ikarus::DofManager {
  VariableVector& operator+=(VariableVector& varVecArg, Eigen::VectorXd& correction) {
    assert(static_cast<long long int>(varVecArg.correctionSize()) == correction.size());
    for (size_t variableIndex = 0; auto&& var : std::ranges::join_view(varVecArg.dofVecimpl))
      var += correction(varVecArg.variableIndices[variableIndex++]);
    return varVecArg;
  }
}