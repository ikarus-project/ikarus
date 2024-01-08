// SPDX-FileCopyrightText: 2021-2024 The Ikarus Developers mueller@ibb.uni-stuttgart.de
// SPDX-License-Identifier: LGPL-3.0-or-later

/**
 * @file powerbasisfe.hh
 * @brief Contains the PowerBasisFE class, which works with a power basis in FlatInterLeaved elements.
 */

#pragma once

#include <ikarus/finiteelements/fetraits.hh>

namespace Ikarus {

  /**
   * @brief PowerBasisFE class for working with a power basis in FlatInterLeaved elements.
   *
   * This class assumes the local indices in a FlatInterLeaved order.
   *
   * @tparam Basis The type of the power basis.
   */
  template <typename Basis>
  class PowerBasisFE {
  public:
    using RootBasis   = Basis;                           ///< Type of the root basis.
    using LocalView   = typename Basis::LocalView;       ///< Type of the local view.
    using GlobalIndex = typename LocalView::MultiIndex;  ///< Type of the global index.
    using GridElement = typename LocalView::Element;     ///< Type of the grid element.
    using Traits      = FETraits<GridElement>;           ///< Type of the traits.

    /**
     * @brief Constructor for the PowerBasisFE class.
     *
     * @param p_basis The power basis.
     * @param element The local element.
     */
    explicit PowerBasisFE(const Basis& p_basis, const typename LocalView::Element& element)
        : localView_{p_basis.localView()} {
      static_assert(Ikarus::Concepts::PowerBasis<RootBasis>,
                    "You didn't pass a localview of a power basis to this method");
      static_assert(RootBasis::PreBasis::Node::degree() != 1,
                    "The basis has only one children. Maybe use scalarFE.hh.");

      localView_.bind(element);
    }

    /** @brief Number of children in the powerBasis. */
    static constexpr int num_children = RootBasis::PreBasis::Node::degree();

    /**
     * @brief Get the size of the local view.
     * @return The size of the local view.
     */
    [[nodiscard]] constexpr size_t size() const { return localView_.size(); }

    /**
     * @brief Get the global flat indices for the power basis.
     *
     * The global indices are collected in a FlatInterLeaved order. I.e. x_0, y_0, z_0, ..., x_n, y_n, z_n
     *
     * @param globalIndices Output vector to store global indices.
     */
    void globalFlatIndices(std::vector<GlobalIndex>& globalIndices) const {
      globalIndices.clear();

      const auto& fe = localView_.tree().child(0).finiteElement();
      for (size_t i = 0; i < fe.size(); ++i) {
        for (int j = 0; j < num_children; ++j) {
          globalIndices.push_back(localView_.index((localView_.tree().child(j).localIndex(i))));
        }
      }
    }

    /**
     * @brief Get the grid element associated with the local view.
     * @return The grid element.
     */
    const GridElement& gridElement() const { return localView_.element(); }

    /**
     * @brief Get the const reference to the local view.
     * @return The const reference to the local view.
     */
    const LocalView& localView() const { return localView_; }

    /**
     * @brief Get the reference to the local view.
     * @return The reference to the local view.
     */
    LocalView& localView() { return localView_; }

  private:
    LocalView localView_;
  };

}  // namespace Ikarus
