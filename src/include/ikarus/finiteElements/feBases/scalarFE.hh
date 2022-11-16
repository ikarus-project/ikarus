/*
 * This file is part of the Ikarus distribution (https://github.com/IkarusRepo/Ikarus).
 * Copyright (c) 2022. The Ikarus developers.
 *
 * This library is free software; you can redistribute it and/or
 *  modify it under the terms of the GNU Lesser General Public
 *  License as published by the Free Software Foundation; either
 *  version 2.1 of the License, or (at your option) any later version.
 *
 * This library is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
 * Lesser General Public License for more details.
 *
 * You should have received a copy of the GNU Lesser General Public
 * License along with this library; if not, write to the Free Software
 * Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301 USA
 */

#pragma once

#include <ikarus/finiteElements/feTraits.hh>

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
        : localView{basis.localView()} {
      static_assert(RootBasis::PreBasis::Node::CHILDREN == 0, "This is no scalar basis!");
      localView.bind(element);
    }

    /** \brief Type of the Pairs of gridEntities and variable tags */
    using GridElementEntityType = typename LocalView::Element;
    using Traits                = FETraits<GridElementEntityType>;

    /** \brief Dimension of the world space */
    static constexpr int worlddim = Traits::worlddim;

    [[nodiscard]] int size() const { return localView.size(); }

    void globalIndices(std::vector<GlobalIndex>& globalIndices) const {
      const auto& fe = localView.tree().finiteElement();
      for (size_t i = 0; i < fe.size(); ++i)
        globalIndices.push_back(localView.index(localView.tree().localIndex(i)));
    }

    const GridElementEntityType& getEntity() { return localView.element(); }

  private:
    LocalView localView;
  };
}  // namespace Ikarus
