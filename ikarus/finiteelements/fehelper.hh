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
 * \brief A helper function to handle global indices of a scalar basis.
 *
 * \tparam LocalView Type of the local view
 * \tparam Child Type of the child
 *
 * \param localView Local view of the element.
 * \param globalIndices Output vector to store global indices.
 * \param child The child of a tree.
 */
template <typename LocalView, typename Child>
void globalScalarFlatIndices(const LocalView& localView, std::vector<typename LocalView::MultiIndex>& globalIndices,
                             const Child& child) {
  const auto& fe = child.finiteElement();
  for (size_t i = 0; i < fe.size(); ++i)
    globalIndices.push_back(localView.index(child.localIndex(i)));
}

/**
 * \brief A helper function to handle global indices of a power basis.
 *
 * \tparam LocalView Type of the local view
 * \tparam Tree Type of the tree
 * \tparam ChildIndex Type of the index for the child, see dune/common/indices.hh
 *
 * \param localView Local view of the element.
 * \param globalIndices Output vector to store global indices.
 * \param tree The tree of the child.
 * \param childIndex The index of the child.
 */
template <typename LocalView, typename Tree, typename ChildIndex>
void globalPowerFlatIndices(const LocalView& localView, std::vector<typename LocalView::MultiIndex>& globalIndices,
                            const Tree& tree, const ChildIndex childIndex) {
  if constexpr (requires { tree.child(childIndex).finiteElement(); }) {
    const auto& fe             = tree.child(childIndex).finiteElement();
    constexpr int childrenSize = Tree::degree();
    for (size_t i = 0; i < fe.size(); ++i)
      for (int j = 0; j < childrenSize; ++j)
        globalIndices.push_back(localView.index(tree.child(j).localIndex(i)));
  } else {
    constexpr int childrenSize = Tree::template Child<childIndex>::Type::degree();
    const auto& fe             = tree.child(childIndex, 0).finiteElement();
    for (size_t i = 0; i < fe.size(); ++i)
      for (int j = 0; j < childrenSize; ++j)
        globalIndices.push_back(localView.index(tree.child(childIndex, j).localIndex(i)));
  }
}

/**
 * \brief A helper function to handle global indices of a composite basis.
 *
 * \tparam LocalView Type of the local view
 * \tparam Tree Type of the tree
 *
 * \param localView Local view of the element.
 * \param globalIndices Output vector to store global indices.
 * \param tree The tree of the child.
 */
template <typename Tree, typename LocalView>
void globalCompositeFlatIndices(const LocalView& localView, std::vector<typename LocalView::MultiIndex>& globalIndices,
                                const Tree& tree) {
  using namespace Dune::Indices;
  Dune::Hybrid::forEach(Dune::Hybrid::integralRange(Dune::index_constant<Tree::degree()>()), [&](const auto i) {
    using ChildType   = std::tuple_element_t<i, typename Tree::ChildTypes>;
    const auto& child = tree.child(i);
    if constexpr (ChildType::isLeaf)
      globalScalarFlatIndices(localView, globalIndices, child);
    else if constexpr (ChildType::isPower) {
      if constexpr (ChildType::template Child<_0>::Type::isComposite) {
        Dune::Hybrid::forEach(
            Dune::Hybrid::integralRange(Dune::index_constant<ChildType::degree()>()), [&](const auto j) {
              globalCompositeFlatIndices<std::remove_cvref_t<typename ChildType::template Child<j>::Type>>(
                  localView, globalIndices, child.child(j));
            });
      } else
        globalPowerFlatIndices(localView, globalIndices, tree, i);
    } else if constexpr (ChildType::isComposite)
      globalCompositeFlatIndices<std::remove_cvref_t<ChildType>>(localView, globalIndices, child);
    else
      DUNE_THROW(Dune::NotImplemented, "globalCompositeFlatIndices is not implemented for the provided basis type.");
  });
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
  using Tree = typename std::remove_cvref_t<LocalView>::Tree;

  /** \brief Bool to check if the basis is a power basis. */
  constexpr bool isPower = Tree::isPower;

  /** \brief Bool to check if the basis is a composite basis. */
  constexpr bool isComposite = Tree::isComposite;

  /** \brief Bool to check if the basis is a scalar basis. */
  constexpr bool isScalar = Tree::isLeaf;

  const auto tree = localView.tree();

  if constexpr (isPower)
    globalPowerFlatIndices(localView, globalIndices, tree, _0);
  else if constexpr (isScalar)
    globalScalarFlatIndices(localView, globalIndices, tree);
  else if constexpr (isComposite)
    globalCompositeFlatIndices<Tree>(localView, globalIndices, tree);
  else
    DUNE_THROW(Dune::NotImplemented, "globalFlatIndices is not implemented for the provided basis type.");
}
} // namespace Ikarus::FEHelper
