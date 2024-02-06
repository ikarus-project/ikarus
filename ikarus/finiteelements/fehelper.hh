// SPDX-FileCopyrightText: 2021-2024 The Ikarus Developers mueller@ibb.uni-stuttgart.de
// SPDX-License-Identifier: LGPL-3.0-or-later

#pragma once

#include <dune/localfefunctions/manifolds/realTuple.hh>

#include <ikarus/finiteelements/ferequirements.hh>
#include <ikarus/finiteelements/physicshelper.hh>

namespace Ikarus::FEHelper {
/**
 * \brief Gets the local solution Dune block vector
 *
 * \tparam Traits Type of the FE traits.
 * \tparam ST Scalar type for the local solution vector.
 *
 * \param x The global solution vector.
 * \param localView Local view of the element.
 * \param dx Optional global solution vector.
 *
 * \return A Dune block vector representing the solution quantities at each node.
 * */
template <typename Traits, typename ST>
auto localSolutionBlockVector(const typename Traits::FERequirementType::SolutionVectorTypeRaw& x,
                              const typename Traits::LocalView& localView,
                              const std::optional<const Eigen::VectorX<ST>>& dx = std::nullopt) {
  constexpr int worldDim = Traits::worlddim;
  const auto& fe         = localView.tree().child(0).finiteElement();
  Dune::BlockVector<Dune::RealTuple<ST, worldDim>> localX(fe.size());
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
} // namespace Ikarus::FEHelper
