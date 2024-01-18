// SPDX-FileCopyrightText: 2021-2024 The Ikarus Developers mueller@ibb.uni-stuttgart.de
// SPDX-License-Identifier: LGPL-3.0-or-later

#pragma once

#include <dune/localfefunctions/manifolds/realTuple.hh>

#include <ikarus/finiteelements/ferequirements.hh>
#include <ikarus/finiteelements/physicshelper.hh>

namespace Ikarus::FEHelper {
  /**
   * @brief Gets the local solution Dune block vector
   *
   * @tparam FERequirementType Type of the FE requirements.
   * @tparam LocalView Type of the local view.
   * @tparam ScalarType Scalar type for the local solution vector.
   *
   * @param x The global solution vector.
   * @param localView Local view of the element.
   * @param dx Optional global solution vector.
   *
   * @return A Dune block vector representing the solution quantities at each node.
   * */
  template <typename FERequirementType, typename LocalView, typename ScalarType>
  auto localSolutionBlockVector(const typename FERequirementType::SolutionVectorTypeRaw& x, const LocalView& localView,
                                const std::optional<const Eigen::VectorX<ScalarType>>& dx = std::nullopt) {
    using Traits           = TraitsFromLocalView<LocalView>;
    constexpr int worldDim = Traits::worlddim;
    const auto& fe         = localView.tree().child(0).finiteElement();
    Dune::BlockVector<Dune::RealTuple<ScalarType, worldDim>> localX(fe.size());
    if (dx)
      for (auto i = 0U; i < localX.size(); ++i)
        for (auto j = 0U; j < worldDim; ++j)
          localX[i][j] = dx.value()[i * worldDim + j] + x[localView.index(localView.tree().child(j).localIndex(i))[0]];
    else
      for (auto i = 0U; i < localX.size(); ++i)
        for (auto j = 0U; j < worldDim; ++j)
          localX[i][j] = x[localView.index(localView.tree().child(j).localIndex(i))[0]];
    return localX;
  }
}  // namespace Ikarus::FEHelper