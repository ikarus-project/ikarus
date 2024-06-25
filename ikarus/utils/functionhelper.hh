// SPDX-FileCopyrightText: 2021-2024 The Ikarus Developers mueller@ibb.uni-stuttgart.de
// SPDX-License-Identifier: LGPL-3.0-or-later

/**
 * \file functionhelper.hh
 * \brief Helper for dune-functions
 */

#pragma once

#include <ranges>

#include <dune/localfefunctions/eigenDuneTransformations.hh>

#include <ikarus/utils/concepts.hh>

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
                "obtainLagrangeNodePositions is only supported for Lagrange power basis");
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
 * \brief A helper function to obtain the global index from the global positions for a Lagrange node
 * \ingroup utils
 * \tparam size Size of the nodal coordinate vector
 * \tparam FEC Type of the finite element container.
 * \tparam CI Type of the child index defining the direction to be fixed
 *
 * \param fes Finite element container.
 * \param pos Global position
 * \param childIndex Index of the child in a power basis defining the direction
 * \return Global index
 */
template <int size, typename FEC, typename CI>
auto globalIndexFromGlobalPosition(const FEC& fes, const Eigen::Vector<double, size>& pos, const CI& childIndex) {
  static_assert(Concepts::LagrangeNode<std::remove_cvref_t<decltype(fes[0].localView().tree().child(0))>>,
                "globalIndexFromGlobalPosition is only supported for Lagrange power basis");
  constexpr double tol = 1e-8;
  typename std::remove_cvref_t<decltype(fes[0])>::LocalView::MultiIndex index{0};
  bool positionFound = false;
  for (const auto& fe : fes) {
    const auto& localView = fe.localView();
    assert(localView.bound() && "The local view must be bound to an element");
    std::vector<Dune::FieldVector<double, size>> lagrangeNodeCoords;
    obtainLagrangeNodePositions(localView, lagrangeNodeCoords);
    const auto& localFE = localView.tree().child(0).finiteElement();
    for (int i = 0; i < localFE.size(); i++)
      if (Dune::toEigen(lagrangeNodeCoords[i]).isApprox(pos, tol)) {
        index = localView.index(localView.tree().child(childIndex).localIndex(i));
        positionFound = true;
      }
  }
  if (not positionFound)
    DUNE_THROW(Dune::InvalidStateException, "No Lagrange node found at the given position");
  return index;
}

/**
 * \brief A function to obtain the local coordinates of subentities of an FiniteElement
 * \ingroup utils
 * \tparam FE Type of the finite element
 * \param fe finite element
 * \param codim codim of requested subentity
 * \return view over the position of the subenties
 */
template <typename FE>
auto referenceElementSubEntityPositions(FE& fe, int codim) {
  constexpr int dim            = FE::Traits::mydim;
  const auto& element          = fe.gridElement();
  const auto& referenceElement = Dune::referenceElement<double, dim>(element.type());
  const int numberOfVertices   = referenceElement.size(codim);

  auto getPosition = [=](const int i) { return referenceElement.position(i, codim); };
  return std::views::transform(std::views::iota(0, numberOfVertices), getPosition);
}

/**
 * \brief A function to obtain the local coordinates the vertices of an FiniteElement
 * \ingroup utils
 * \tparam FE Type of the finite element
 * \param fe finite element
 * \return
 */
template <typename FE>
auto referenceElementVertexPositions(FE& fe) {
  return referenceElementSubEntityPositions(fe, FE::Traits::mydim);
}

} // namespace Ikarus::utils
