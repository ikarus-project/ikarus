// SPDX-FileCopyrightText: 2021-2026 The Ikarus Developers ikarus@ibb.uni-stuttgart.de
// SPDX-License-Identifier: LGPL-3.0-or-later

/**
 * \file functionhelper.hh
 * \brief Helper for dune-functions
 */

#pragma once

#include <ranges>

#include <dune/common/float_cmp.hh>
#include <dune/grid/utility/hierarchicsearch.hh>

#include <ikarus/utils/concepts.hh>
#include <ikarus/utils/flatprebasis.hh>
#include <ikarus/utils/traversal.hh>

namespace Ikarus::utils {
/**
 * \brief A function to obtain the global positions of the nodes of an element with Lagrangian basis
 * \ingroup utils
 * \tparam size Size of the global nodal coordinate vector
 * \tparam LV Type of the local view
 *
 * \param localView Local view bounded to an element
 * \param lagrangeNodeGlobalCoords A vector of global nodal coordinates to be updated
 */
template <int size, typename LV>
void obtainLagrangeGlobalNodePositions(const LV& localView,
                                       std::vector<Dune::FieldVector<double, size>>& lagrangeNodeGlobalCoords) {
  auto fT = [&]([[maybe_unused]] int nodeNumber, Dune::FieldVector<double, size>&& localCoordinate) {
    lagrangeNodeGlobalCoords.emplace_back(localView.element().geometry().global(localCoordinate));
    return false;
  };
  forEachLagrangeNodePosition(localView, fT);
}

/**
 * \brief A helper function to obtain the global index from the global positions for a Lagrange node
 * \ingroup utils
 * \tparam worldDim Dimension of the world space.
 * \tparam Basis Type of the basis.
 * \tparam cmpStyle Comparison style being used to check if a node exists at a certain position.
 *
 * \details As per
 * https://gitlab.mn.tu-dresden.de/spraetor/dune-common-extensions/-/blob/master/dune/common/float_cmp.hh#L58, relative
 * comparison styles do not work properly if one of the values being compared is zero. Hence, for this function,
 * absolute is being used as the default comparison style.
 *
 * \param basis The grid basis.
 * \param pos Global position
 * \return Global index
 */
template <int worldDim, typename Basis, Dune::FloatCmp::CmpStyle cmpStyle = Dune::FloatCmp::CmpStyle::absolute>
auto globalIndexFromGlobalPosition(const Basis& basis, const Dune::FieldVector<double, worldDim>& pos) {
  static_assert(Concepts::LagrangeNode<std::remove_cvref_t<decltype(basis.localView().tree().child(0))>>,
                "globalIndexFromGlobalPosition is only supported for Lagrange basis");
  constexpr double tol = 1e-8;
  using LocalView      = std::remove_cvref_t<decltype(basis.localView())>;
  using MultiIndex     = typename LocalView::MultiIndex;
  using Element        = typename LocalView::Element;
  using Tree           = LocalView::Tree;

  if constexpr (Tree::isComposite) {
    Dune::Hybrid::forEach(
        Dune::Hybrid::integralRange(Dune::index_constant<std::tuple_size_v<typename Basis::PreBasis::SubPreBases>>()),
        [&](const auto i) {
          static_assert(
              not Tree::template Child<i>::Type::isComposite,
              "globalIndexFromGlobalPosition is not implemented to handle a composite basis within a composite basis.");
        });
  }

  static constexpr std::size_t numChildren = Impl::PreBasisInfo<Tree>::size;
  constexpr int myDim                      = Element::mydimension;
  Dune::HierarchicSearch hSearch(basis.gridView().grid(), basis.gridView().indexSet());
  const auto& ele = hSearch.findEntity(pos);
  auto localView  = basis.localView();
  localView.bind(ele);
  const auto geo   = localView.element().geometry();
  const auto& node = localView.tree();
  std::optional<std::array<MultiIndex, numChildren>> globalIndices;

  auto fT = [&](int nodeNumber, Dune::FieldVector<double, myDim>&& localCoordinate) {
    if (Dune::FloatCmp::eq<Dune::FieldVector<double, worldDim>, cmpStyle>(geo.global(localCoordinate), pos, tol)) {
      globalIndices.emplace();
      for (int j = 0; j < numChildren; j++)
        globalIndices.value()[j] = localView.index(node.child(j).localIndex(nodeNumber));
      return true;
    }
    return false;
  };
  forEachLagrangeNodePosition(localView, fT);
  if (not globalIndices.has_value())
    DUNE_THROW(Dune::InvalidStateException, "No Lagrange node found at the given position in the grid.");
  return globalIndices.value();
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

/**
 * \brief if T is a pointer type, return the dereferenced value, otherwise return the value itself.
 *
 * \tparam T
 * \param t
 * \return decltype(auto)
 */
template <typename T>
decltype(auto) maybeDeref(T& t) {
  if constexpr (Concepts::PointerOrSmartPointer<T>)
    return *t;
  else
    return t;
}

/**
 * \brief A helper function to compute the forces due to inhomogeneous Dirichlet BCs.
 * \tparam A Type of the assembler.
 * \param assembler The assembler.
 */
template <typename A>
requires(std::is_pointer_v<std::remove_cvref_t<A>> || traits::isSharedPtr<std::remove_cvref_t<A>>::value)
auto obtainForcesDueToIDBC(const A& assembler) {
  const auto& K = assembler->matrix(DBCOption::Raw);
  auto& dv      = assembler->dirichletValues();
  auto newInc   = Eigen::VectorXd::Zero(dv.size()).eval();
  dv.evaluateInhomogeneousBoundaryConditionDerivative(newInc, 1.0);

  Eigen::VectorXd F_dirichlet;
  F_dirichlet = K * newInc;

  if (assembler->dBCOption() == DBCOption::Full)
    assembler->dirichletValues().setZeroAtConstrainedDofs(F_dirichlet);
  else
    F_dirichlet = assembler->createReducedVector(F_dirichlet);
  return F_dirichlet;
}

} // namespace Ikarus::utils
