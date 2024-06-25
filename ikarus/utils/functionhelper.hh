// SPDX-FileCopyrightText: 2021-2024 The Ikarus Developers mueller@ibb.uni-stuttgart.de
// SPDX-License-Identifier: LGPL-3.0-or-later

/**
 * \file functionhelper.hh
 * \brief Helper for dune-functions
 */

#pragma once

#include <ranges>

#include <dune/grid/utility/hierarchicsearch.hh>
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
template <int size, typename LV,typename F>
void obtainLagrangeGlobalNodePositions(const LV& localView,std::vector<Dune::FieldVector<double, size>>& lagrangeNodeGlobalCoords) {
  std::vector<Dune::FieldVector<double, size>> lagrangeNodeGlobalCoords;
  auto fT= [&](const auto& localFE,int i, auto&& localCoordinate)
  {
    lagrangeNodeGlobalCoords.emplace_back(localCoordinate);
    return false;
  };
  forEachLagrangeNodePosition(localView,fT);
}




template <int size, typename LV,typename F>
void forEachLagrangeNodePosition(const LV& localView,F&& f) {
  static_assert(Concepts::LagrangeNode<std::remove_cvref_t<decltype(localView.tree().child(0))>>,
                "obtainLagrangeNodePositions is only supported for Lagrange power basis");
  assert(localView.bound() && "The local view must be bound to an element");
  const auto& localFE = localView.tree().child(0).finiteElement();
  std::vector<Dune::FieldVector<double, size>> lagrangeNodeCoords;
  lagrangeNodeCoords.resize(localFE.size());
  std::vector<double> out;
  for (int i = 0; i < size; i++) {
    std::vector<double> out;
    auto ithCoord = [&i](const Dune::FieldVector<double, size>& x) { return x[i]; };
    const auto& localFE = localView.tree().child(0).finiteElement();
    localFE.localInterpolation().interpolate(ithCoord, out);
    for (std::size_t j = 0; j < out.size(); j++)
      lagrangeNodeCoords[j][i] = out[j];
  }
  for ( auto i=0;auto& nCoord : lagrangeNodeCoords)
      if(f(localFE,i++,nCoord))
        break;
  // return  std::ranges::transform_view(std::ranges::iota_view(0,out.size()),[&](int i){  std::vector<double> out;
  // for (int i = 0; i < size; i++) {
  //   auto ithCoord = [&i](const Dune::FieldVector<double, size>& x) { return x[i]; };
  //   localFE.localInterpolation().interpolate(ithCoord, out);
  //   for (std::size_t j = 0; j < out.size(); j++)
  //     lagrangeNodeCoords[j][i] = out[j];
  // }});
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
template <int size, typename Basis, typename CI>
auto globalIndexFromGlobalPosition(const Basis& basis, const Dune::Field<double, size>& requiredPosition) {
  static_assert(Concepts::LagrangeNode<std::remove_cvref_t<decltype(fes[0].localView().tree().child(0))>>,
                "globalIndexFromGlobalPosition is only supported for Lagrange power basis");
  constexpr double tol = 1e-8;
  typename std::remove_cvref_t<decltype(fes[0])>::LocalView::MultiIndex index{0};
  bool positionFound = false;
  Dune:Dune::HierarchicSearch hSearch(basis.gridView().grid(),basis.gridView().indexSet());
  const auto e = hSearch.findEntity(requiredPosition);
  const auto& localView = e.localView();
  std::optional<std::array<  typename LV::MultiIndex , size>> globalIndices;
  const auto geo= localView.element().geometry();
  const auto& node = localView.tree();

  auto fT= [&](const auto& localFE,int i, auto&& localCoordinate)
  {
    if (Dune::FloatCmp::eq(geo.global(localCoordinate),requiredPosition,tol)) {
      globalIndices.emplace();
      for (int j = 0; j < size; j++) {
        globalIndices.value()[j]= localView.index(node.child(j).localIndex(i));
      }
      return true;
    }else
      return false;
  };
  forEachLagrangeNodePosition(localView,fT);
  return globalIndices;
}

template <int size, typename LV,typename F>
void obtainLagrangeGlobalNodePositions(const LV& localView,const Dune::FieldVector<double, size>& requiredPosition, double tol) {
  std::vector<Dune::FieldVector<double, size>> lagrangeNodeGlobalCoords;
  const auto geo= localView.element().geometry();
  const auto& node = localView.tree();


  return globalIndices;
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
