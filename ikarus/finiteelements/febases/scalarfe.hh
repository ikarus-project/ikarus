// SPDX-FileCopyrightText: 2021-2022 The Ikarus Developers mueller@ibb.uni-stuttgart.de
// SPDX-License-Identifier: LGPL-3.0-or-later

#pragma once

#include <ikarus/finiteelements/fetraits.hh>

namespace Ikarus {

  /**
   * scalarFE.hh includes a class which takes in a scalar basis and arranges
   * the indices one after the other
   * Single-DOF elements can be inherited from this class which uses a scalar basis
   */
  template <typename Basis>
  class ScalarFieldFE {
  public:
    using RootBasis   = Basis;
    using LocalView   = typename Basis::LocalView;
    using GlobalIndex = typename LocalView::MultiIndex;
    explicit ScalarFieldFE(const Basis& basis, const typename LocalView::Element& element)
        : localView_{basis.localView()} {
      static_assert(RootBasis::PreBasis::Node::CHILDREN == 0, "This is no scalar basis!");
      localView_.bind(element);
    }

    /** \brief Type of the Pairs of gridEntities and variable tags */
    using GridElementEntityType = typename LocalView::Element;
    using Traits                = FETraits<GridElementEntityType>;

    /** \brief Dimension of the world space */
    static constexpr int worlddim = Traits::worlddim;

    [[nodiscard]] int size() const { return localView_.size(); }

    void globalFlatIndices(std::vector<GlobalIndex>& globalIndices) const {
      const auto& fe = localView_.tree().finiteElement();
      for (size_t i = 0; i < fe.size(); ++i)
        globalIndices.push_back(localView_.index(localView_.tree().localIndex(i)));
    }

    const GridElementEntityType& getEntity() const { return localView_.element(); }
    const LocalView& localView() const { return localView_; }
    LocalView& localView() { return localView_; }

  private:
    LocalView localView_;
  };
}  // namespace Ikarus
