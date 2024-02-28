// SPDX-FileCopyrightText: 2021-2024 The Ikarus Developers mueller@ibb.uni-stuttgart.de
// SPDX-License-Identifier: LGPL-3.0-or-later

/**
 * \file functionhelper.hh
 * \brief Helper for dune-functions
 */

#pragma once

#include <ranges>

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
/**
 * \brief A function to obtain the local coordinates of subentities of an FiniteElement
 * \ingroup utils
 * \tparam FE Type of the finite element
 * \param fe finite element
 * \param codim codim of requested subentity, defaults to vertices
 * \return view over the position of the subenties
 */
template <typename FE>
auto referenceElementSubEntitiyPositions(FE& fe, const int codim = FE::Traits::mydim) {
  constexpr int dim            = FE::Traits::mydim;
  const auto& element          = fe.gridElement();
  const auto& referenceElement = Dune::referenceElement<double, dim>(element.type());
  const int numberOfVertices   = referenceElement.size(codim);

  auto getPosition = [=](const int i) { return referenceElement.position(i, codim); };
  return std::views::transform(std::views::iota(0, numberOfVertices), getPosition);
}

} // namespace Ikarus::utils
