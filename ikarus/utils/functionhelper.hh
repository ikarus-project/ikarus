// SPDX-FileCopyrightText: 2021-2024 The Ikarus Developers mueller@ibb.uni-stuttgart.de
// SPDX-License-Identifier: LGPL-3.0-or-later

/**
 * \file functionhelper.hh
 * \brief Helper for dune-functions
 */

#pragma once

namespace Ikarus::utils {
/**
 * \brief A function to obtain the global positions of the nodes of an element with Lagrangian basis, see Dune book
 * page 314
 * \ingroup utils
 * \tparam size Size of the nodal coordinate vector
 * \tparam LV Type of the local view
 *
 * \param localView Local view bounded to an element
 * \param lagrangeNodeCoords A vector of nodal coordinates to be updated
 */
template <int size, typename LV>
void obtainLagrangeNodePositions(const LV& localView,
                                 std::vector<Dune::FieldVector<double, size>>& lagrangeNodeCoords) {
  static_assert(Concepts::LagrangeNode<std::remove_cvref_t<decltype(localView.tree().child(0))>>,
                "This function is only supported for Lagrange basis");
  assert(localView.bound() && "The local view must be bound to an element");
  const auto& localFE = localView.tree().child(0).finiteElement();
  lagrangeNodeCoords.resize(localFE.size());
  std::vector<double> out;
  for (int i = 0; i < size; i++) {
    auto ithCoord = [&i](const Dune::FieldVector<double, size>& x) { return x[i]; };
    localFE.localInterpolation().interpolate(ithCoord, out);
    for (std::size_t j = 0; j < out.size(); j++)
      lagrangeNodeCoords[j][i] = out[j];
  }
  for (auto& nCoord : lagrangeNodeCoords)
    nCoord = localView.element().geometry().global(nCoord);
}

} // namespace Ikarus::utils
