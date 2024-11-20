// SPDX-FileCopyrightText: 2021-2025 The Ikarus Developers mueller@ibb.uni-stuttgart.de
// SPDX-License-Identifier: LGPL-3.0-or-later

#pragma once

#include <functional>
#include <optional>

#include <dune/localfefunctions/manifolds/realTuple.hh>

#include <ikarus/utils/traversal.hh>

namespace Ikarus::FEHelper {
/**
 * \brief Gets the local solution Dune block vector
 *
 * \tparam Traits Type of the FE traits.
 * \tparam ST Scalar type for the local solution vector.
 * \tparam Vector Global solution vector
 *
 * \param x The global solution vector.
 * \param localView Local view of the element.
 * \param dx Optional global solution vector.
 *
 * \return A Dune block vector representing the solution quantities at each node.
 * */
template <typename Traits, typename Vector, typename ST>
auto localSolutionBlockVector(
    const Vector& x, const typename Traits::LocalView& localView,
    const std::optional<std::reference_wrapper<const Eigen::VectorX<ST>>>& dx = std::nullopt) {
  constexpr int worldDim = Traits::worlddim;
  const auto& fe         = localView.tree().child(0).finiteElement();
  Dune::BlockVector<Dune::RealTuple<ST, worldDim>> localX(fe.size());
  if (dx) {
    for (auto i = 0U; i < localX.size(); ++i)
      for (auto j = 0U; j < worldDim; ++j)
        localX[i][j] =
            dx.value().get()[i * worldDim + j] + x[localView.index(localView.tree().child(j).localIndex(i))[0]];
  } else {
    for (auto i = 0U; i < localX.size(); ++i)
      for (auto j = 0U; j < worldDim; ++j)
        localX[i][j] = x[localView.index(localView.tree().child(j).localIndex(i))[0]];
  }

  return localX;
}

namespace Impl {
  /**
   * \brief A helper function to handle global indices of a scalar basis at a leaf node.
   *
   * \tparam LocalView Type of the local view
   * \tparam Node Type of the leaf node
   *
   * \param localView Local view of the element.
   * \param node Leaf node of a tree.
   * \param globalIndices Output vector to store global indices.
   */
  template <typename LocalView, typename Node>
  void leafNodeIndices(const LocalView& localView, const Node& node,
                       std::vector<typename LocalView::MultiIndex>& globalIndices) {
    const auto& fe = node.finiteElement();
    for (size_t i = 0; i < fe.size(); ++i)
      globalIndices.push_back(localView.index(node.localIndex(i)));
  }

  /**
   * \brief A helper function to handle global indices of a power basis at a power node.
   *
   * \tparam LocalView Type of the local view
   * \tparam Node Type of the power node
   *
   * \param localView Local view of the element.
   * \param node Power node of a tree.
   * \param globalIndices Output vector to store global indices.
   */
  template <typename LocalView, typename Node>
  void powerNodeIndices(const LocalView& localView, const Node& node,
                        std::vector<typename LocalView::MultiIndex>& globalIndices) {
    const auto& fe         = node.child(0).finiteElement();
    const int childrenSize = node.degree();
    for (size_t i = 0; i < fe.size(); ++i)
      for (int j = 0; j < childrenSize; ++j)
        globalIndices.push_back(localView.index(node.child(j).localIndex(i)));
  }
} // namespace Impl

/**
 * \brief Get the global indices for the provided local view of an element.
 *
 * \tparam LocalView Type of the local view
 *
 * \param localView Local view of the element.
 * \param globalIndices Output vector to store global indices.
 */
template <typename LocalView>
void globalIndicesFromLocalView(const LocalView& localView,
                                std::vector<typename LocalView::MultiIndex>& globalIndices) {
  assert(localView.bound());
  globalIndices.clear();
  using namespace Dune::Indices;
  using namespace FEHelper::Impl;

  auto leafOpFunc = [&](auto&& node, [[maybe_unused]] auto&& treePath) {
    leafNodeIndices(localView, node, globalIndices);
  };

  auto powerOpFunc = [&](auto&& node, [[maybe_unused]] auto&& treePath) {
    powerNodeIndices(localView, node, globalIndices);
  };

  utils::forEachLeafOrPowerLeafNode(localView.tree(), Dune::TypeTree::hybridTreePath(), powerOpFunc, leafOpFunc);
}

/**
 * \brief Get the global indices for the provided finite element.
 *
 * \details The global indices are collected in a FlatInterLeaved order or in BlockedInterleaved order.
 * This function can handle a scalar basis, power basis, and a composite basis.
 *
 * \tparam FiniteElement Type of the finite element.
 *
 * \param fe The finite element.
 * \param globalIndices Output vector to store global indices.
 */
template <typename FiniteElement>
void globalIndices(const FiniteElement& fe, std::vector<typename FiniteElement::LocalView::MultiIndex>& globalIndices) {
  globalIndicesFromLocalView(fe.localView(), globalIndices);
}

/**
 * \brief Maps the entries from the localMatrix of the size \f$ n_s \f$ by \f$ n_s \f$ (where \f$ n_s \f$ is the number
 of
 * coefficiets of a subspace basis) to the local index in the element matrix.
 *
 * \tparam TreePaths type of the indices used to obtain the subspacebasis which corresponds to the localMatrix, can be a
 * `index_constant`, an integer, or a `Dune::TypeTree::HybridTree` (one or more are accepted corresponding to the
 * localMatrix)
 * \param localMatrix the localMatrix containing the entries that are beeing mapped.
 * \param elementMatrix the elementMatrix where the entries are beeing mapped to.
 * \param globalBasis the global basis, from where we obtain the subspacebasis.
 * \param tps the paths to the subspacebasis'.
 */
template <typename... TreePaths>
void subspaceEntriesToElementEntries(const auto& localMatrix, auto& elementMatrix, const auto& globalBasis,
                                     TreePaths&&... tps) {
  auto treePathTuple = Dune::TupleVector(std::forward<TreePaths>(tps)...);
  for (size_t i = 0; i < localMatrix.rows(); ++i) {
    for (size_t j = 0; j < localMatrix.cols(); ++j) {
      Dune::Hybrid::forEach(treePathTuple, [&](const auto& treePath) {
        auto basis = globalBasis.child(treePath);
        elementMatrix(basis.localIndex(i), basis.localIndex(j)) += localMatrix(i, j);
      });
    }
  }
}

/**
 * \brief Forwards the indices for powerbasis to obtain the children to subspaceEntriesToElementEntries()
 *
 * \tparam nEnd index from the last child + 1
 * \tparam nStart optionnally define a starting index
 * \details see Documentation for subspaceEntriesToElementEntries()
 */
template <int nEnd, int nStart = 0>
void subspaceEntriesToElementEntries(const auto& localMatrix, auto& elementMatrix, const auto& globalBasis) {
  Dune::Hybrid::forEach(
      Dune::Hybrid::integralRange(Dune::index_constant<nStart>(), Dune::index_constant<nEnd>()),
      [&](const auto i) { subspaceEntriesToElementEntries(localMatrix, elementMatrix, globalBasis, i); });
}

} // namespace Ikarus::FEHelper
