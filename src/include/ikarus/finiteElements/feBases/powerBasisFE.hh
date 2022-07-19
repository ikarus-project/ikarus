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

#include "src/include/ikarus/finiteElements/feTraits.hh"
namespace Ikarus {

  /**
   * powerBasisFE.hh includes a class which takes in a power basis and it assumes
   * the local indices in a FlatInterLeaved
   * Multi-DOFs elements can be inherited from this class which uses a power basis
   */
  template <typename Basis>
  class PowerBasisFE {
  public:
    using RootBasis   = Basis;
    using LocalView   = typename Basis::LocalView;
    using GlobalIndex = typename LocalView::MultiIndex;
    explicit PowerBasisFE(const Basis& p_basis, const typename LocalView::Element& element)
        : localView{p_basis.localView()} {
      static_assert(Ikarus::Concepts::PowerBasis<RootBasis>,
                    "You didn't pass a localview of a power basis to this method");
      static_assert(RootBasis::PreBasis::Node::CHILDREN != 1,
                    "The basis has only one children. Maybe use scalarFE.hh. ");

      localView.bind(element);
    }

    /** \brief Type of the Pairs of gridEntities and variable tags */
    using GridElementEntityType = typename LocalView::Element;
    using Traits                = FETraits<GridElementEntityType>;

    /** \brief Number of children in the powerBasis */
    static constexpr int num_children = RootBasis::PreBasis::Node::CHILDREN;

    [[nodiscard]] constexpr int size() const { return localView.size(); }

    void globalIndices(std::vector<GlobalIndex>& globalIndices) const {
      static_assert(
          requires { localView.tree().child(0); },
          "Your basis does not provide a child accessor. Maybe use scalarFE.hh.");
      const auto& fe = localView.tree().child(0).finiteElement();
      for (size_t i = 0; i < fe.size(); ++i) {
        for (int j = 0; j < num_children; ++j) {
          globalIndices.push_back(localView.index((localView.tree().child(j).localIndex(i))));
        }
      }
    }

    const GridElementEntityType& getEntity() { return localView.element(); }

  private:
    LocalView localView;
  };
}  // namespace Ikarus