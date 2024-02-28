// SPDX-FileCopyrightText: 2021-2024 The Ikarus Developers mueller@ibb.uni-stuttgart.de
// SPDX-License-Identifier: LGPL-3.0-or-later

/**
 * \file traversal.hh
 * \brief Contains functions to traverse through a tree to its different nodes
 */

#pragma once

namespace Ikarus::utils {

/**
 * \brief A function which loops over all the nodes of a tree and performs different actions
 * for a power node (with leaf node as child) and a leaf node depending on the corresponding functor passed.
 *
 * \details This function is inspired from the function Dune::TypeTree::Detail::forEachNode
 * available in dune/typetree/traversal.hh
 *
 * \tparam T Type of tree.
 * \tparam TreePath Type of the HybridTreePath.
 * \tparam PowerFunc Type of the functor called for a power node.
 * \tparam LeafFunc Type of the functor called for a leaf node.
 * \param tree Tree whose nodes are to be looped over.
 * \param treePath The tree path which handles index values.
 * \param powerFunc A functor to be called for a power node.
 * \param leafFunc A functor to be called for a leaf node.
 */
template <class T, class TreePath, class PowerFunc, class LeafFunc>
void forEachLeafOrPowerLeafNode(T&& tree, TreePath&& treePath, PowerFunc&& powerFunc, LeafFunc&& leafFunc) {
  using Tree = std::decay_t<T>;
  if constexpr (Tree::isLeaf) {
    leafFunc(tree, treePath);
  } else if constexpr (Tree::isPower) {
    if constexpr (Tree::template Child<Dune::Indices::_0>::Type::isLeaf) {
      powerFunc(tree, treePath);
    } else {
      for (std::size_t i = 0; i < tree.degree(); ++i) {
        auto childTreePath = Dune::TypeTree::push_back(treePath, i);
        forEachLeafOrPowerLeafNode(tree.child(i), childTreePath, powerFunc, leafFunc);
      }
    }
  } else {
    auto indices = std::make_index_sequence<Tree::degree()>{};
    Dune::Hybrid::forEach(indices, [&](auto i) {
      auto childTreePath = Dune::TypeTree::push_back(treePath, i);
      forEachLeafOrPowerLeafNode(tree.child(i), childTreePath, powerFunc, leafFunc);
    });
  }
}
} // namespace Ikarus::utils
