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

/**
 * \brief A function which loops over all the nodes of a tree and performs different actions
 * for a power node and a leaf node depending on the corresponding functor passed.
 *
 * \details This function is inspired from the function Dune::TypeTree::Detail::forEachNode
 * available in dune/typetree/traversal.hh
 *
 * \tparam T Type of tree.
 * \tparam ChildIndex Type of the index of the child.
 * \tparam PowerFunc Type of the functor called for a power node.
 * \tparam LeafFunc Type of the functor called for a leaf node.
 * \param tree Tree whose nodes are to be looped over.
 * \param childIndex The index of the child.
 * \param powerFunc A functor to be called for a power node.
 * \param leafFunc A functor to be called for a leaf node.
 */
template <class T, class ChildIndex, class PowerFunc, class LeafFunc>
void forEachLeafOrPowerNode(T&& tree, ChildIndex&& childIndex, PowerFunc&& powerFunc, LeafFunc&& leafFunc) {
  using Tree = std::decay_t<T>;
  if constexpr (Tree::isLeaf) {
    leafFunc(tree);
  } else if constexpr (Tree::isPower) {
    if constexpr (Tree::template Child<Dune::Indices::_0>::Type::isComposite) {
      auto indices = std::make_index_sequence<Tree::degree()>{};
      Dune::Hybrid::forEach(indices, [&](auto i) { forEachLeafOrPowerNode(tree.child(i), i, powerFunc, leafFunc); });
    } else
      powerFunc(tree, childIndex);
  } else {
    auto indices = std::make_index_sequence<Tree::degree()>{};
    Dune::Hybrid::forEach(indices, [&](auto i) { forEachLeafOrPowerNode(tree.child(i), i, powerFunc, leafFunc); });
  }
}

/**
 * \brief Get the global flat indices for the provided basis.
 *
 * \details The global indices are collected in a FlatInterLeaved order. I.e. x_0, y_0, z_0, ..., x_n, y_n, z_n
 * This function can handle a scalar basis, power basis, and a composite basis.
 * The maximum number of children for the composite basis is 2.
 *
 * \tparam LocalView Type of the local view
 *
 * \param localView Local view of the element.
 * \param globalIndices Output vector to store global indices.
 */
template <typename LocalView>
void globalFlatIndices(const LocalView& localView, std::vector<typename LocalView::MultiIndex>& globalIndices) {
  globalIndices.clear();
  using namespace Dune::Indices;

  /** \brief Functor to obtain global indices of a leaf node. */
  auto leafOpFunc = [&](auto&& node) {
    const auto& fe = node.finiteElement();
    for (size_t i = 0; i < fe.size(); ++i) {
      globalIndices.push_back(localView.index(node.localIndex(i)));
    }
  };

  /** \brief Functor to obtain global indices of a power node. */
  auto powerOpFunc = [&](auto&& node, auto&& ci) {
    if constexpr (requires { node.child(ci).finiteElement(); }) {
      const auto& fe         = node.child(ci).finiteElement();
      const int childrenSize = node.degree();
      for (size_t i = 0; i < fe.size(); ++i)
        for (int j = 0; j < childrenSize; ++j)
          globalIndices.push_back(localView.index(node.child(j).localIndex(i)));
    } else {
      const int childrenSize = node.child(ci, 0).degree();
      const auto& fe         = node.child(ci, 0).finiteElement();
      for (size_t i = 0; i < fe.size(); ++i)
        for (int j = 0; j < childrenSize; ++j)
          globalIndices.push_back(localView.index(node.child(ci, j).localIndex(i)));
    }
  };
  forEachLeafOrPowerNode(localView.tree(), _0, powerOpFunc, leafOpFunc);
}
} // namespace Ikarus::FEHelper
