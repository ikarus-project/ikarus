// SPDX-FileCopyrightText: 2021-2024 The Ikarus Developers mueller@ibb.uni-stuttgart.de
// SPDX-License-Identifier: LGPL-3.0-or-later

/**
 * \file febases.hh
 * \brief Contains the FEBases class, which works with different types of bases in FlatInterLeaved elements.
 */

#pragma once

#include <ikarus/finiteelements/fetraits.hh>

namespace Ikarus {

/**
 * \brief FEBases class for working with different types of bases in FlatInterLeaved elements.
 *
 * \details A single class which can be used by the finite elements to handle scalar basis, power basis, and composite
 * basis. This class assumes the local indices in a FlatInterLeaved order.
 *
 * \tparam B The type of the basis.
 */
template <typename B>
class FEBases
{
public:
  using Traits      = FETraits<B>;                  ///< Type of the traits.
  using LocalView   = typename Traits::LocalView;   ///< Type of the local view.
  using GlobalIndex = typename Traits::GlobalIndex; ///< Type of the global index.
  using GridElement = typename Traits::Element;     ///< Type of the grid element.

  /** \brief Bool to check if the basis is a power basis. */
  static constexpr bool isPower = B::PreBasis::Node::isPower;

  /** \brief Bool to check if the basis is a composite basis. */
  static constexpr bool isComposite = B::PreBasis::Node::isComposite;

  /** \brief Bool to check if the basis is a scalar basis. */
  static constexpr bool isScalar = B::PreBasis::Node::isLeaf;

  /** \brief Number of children in the basis. */
  static constexpr int numChildren = B::PreBasis::Node::degree();

  /**
   * \brief Constructor for the FEBases class.
   *
   * \param basis The basis.
   * \param element The local element.
   */
  explicit FEBases(const B& basis, const typename LocalView::Element& element)
      : localView_{basis.flat().localView()} {
    localView_.bind(element);
  }

  /**
   * \brief Get the size of the local view.
   * \return The size of the local view.
   */
  [[nodiscard]] constexpr size_t size() const { return localView_.size(); }

  /**
   * \brief Get the global flat indices for the provided basis.
   *
   * \details The global indices are collected in a FlatInterLeaved order. I.e. x_0, y_0, z_0, ..., x_n, y_n, z_n
   * This function can handle a scalar basis, power basis, and a composite basis.
   * The maximum number of children for the composite basis is 2.
   *
   * \param globalIndices Output vector to store global indices.
   */
  void globalFlatIndices(std::vector<GlobalIndex>& globalIndices) const { globalFlatIndicesImpl(globalIndices); }

  /**
   * \brief Get the grid element associated with the local view.
   * \return The grid element.
   */
  const GridElement& gridElement() const { return localView_.element(); }

  /**
   * \brief Get the const reference to the local view.
   * \return The const reference to the local view.
   */
  const LocalView& localView() const { return localView_; }

  /**
   * \brief Get the reference to the local view.
   * \return The reference to the local view.
   */
  LocalView& localView() { return localView_; }

  /**
   * \brief Get the number of children of the basis.
   *
   * \return The number of children of the basis.
   */
  [[nodiscard]] int numberOfChildren() const { return numChildren; }

private:
  LocalView localView_;

  /**
   * \brief A helper function to handle global indices of a scalar basis.
   *
   * \tparam ChildIndex Type of the index for the child, see dune/common/indices.hh
   *
   * \param globalIndices Output vector to store global indices.
   * \param childIndex The index of the child.
   */
  template <typename ChildIndex>
  void globalScalarFlatIndices(std::vector<GlobalIndex>& globalIndices, const ChildIndex childIndex) const {
    if constexpr (requires { localView_.tree().child(childIndex).finiteElement(); }) {
      const auto& fe = localView_.tree().child(childIndex).finiteElement();
      for (size_t i = 0; i < fe.size(); ++i)
        globalIndices.push_back(localView_.index((localView_.tree().child(childIndex).localIndex(i))));
    } else {
      const auto& fe = localView_.tree().finiteElement();
      for (size_t i = 0; i < fe.size(); ++i)
        globalIndices.push_back(localView_.index(localView_.tree().localIndex(i)));
    }
  }

  /**
   * \brief A helper function to handle global indices of a power basis.
   *
   * \tparam ChildIndex Type of the index for the child, see dune/common/indices.hh
   *
   * \param globalIndices Output vector to store global indices.
   * \param childIndex The index of the child.
   */
  template <typename ChildIndex>
  void globalPowerFlatIndices(std::vector<GlobalIndex>& globalIndices, const ChildIndex childIndex) const {
    if constexpr (requires { localView_.tree().child(childIndex).finiteElement(); }) {
      const auto& fe             = localView_.tree().child(0).finiteElement();
      constexpr int childrenSize = LocalView::Tree::degree();
      for (size_t i = 0; i < fe.size(); ++i)
        for (int j = 0; j < childrenSize; ++j)
          globalIndices.push_back(localView_.index((localView_.tree().child(j).localIndex(i))));
    } else {
      const auto& fe             = localView_.tree().child(childIndex, 0).finiteElement();
      constexpr int childrenSize = LocalView::Tree::template Child<childIndex>::Type::degree();
      for (size_t i = 0; i < fe.size(); ++i)
        for (int j = 0; j < childrenSize; ++j)
          globalIndices.push_back(localView_.index((localView_.tree().child(childIndex, j).localIndex(i))));
    }
  }

  /**
   * \brief A helper function to handle global indices of a composite basis.
   *
   * \tparam ChildIndex Type of the index for the child, see dune/common/indices.hh
   *
   * \param globalIndices Output vector to store global indices.
   * \param childIndex The index of the child.
   */
  template <typename ChildIndex>
  void globalCompositeFlatIndices(std::vector<GlobalIndex>& globalIndices, const ChildIndex childIndex) const {
    constexpr int localNumChildren = LocalView::Tree::template Child<childIndex>::Type::degree();
    if constexpr (localNumChildren == 0)
      globalScalarFlatIndices(globalIndices, childIndex);
    else
      globalPowerFlatIndices(globalIndices, childIndex);
  }

  void globalFlatIndicesImpl(std::vector<GlobalIndex>& globalIndices) const {
    globalIndices.clear();
    using namespace Dune::Indices;
    if constexpr (isPower)
      globalPowerFlatIndices(globalIndices, _0);
    else if constexpr (isScalar)
      globalScalarFlatIndices(globalIndices, _0);
    else if constexpr (isComposite) {
      Dune::Hybrid::forEach(Dune::Hybrid::integralRange(Dune::index_constant<numChildren>()), [&](const auto i) {
        static_assert(not(LocalView::Tree::template Child<i>::Type::isComposite),
                      "Composite Basis cannot have a composite basis");
        globalCompositeFlatIndices(globalIndices, i);
      });
    } else
      DUNE_THROW(Dune::NotImplemented, "The provided basis type is not implemented.");
  }
};
} // namespace Ikarus
