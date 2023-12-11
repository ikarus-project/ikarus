// SPDX-FileCopyrightText: 2021-2023 The Ikarus Developers mueller@ibb.uni-stuttgart.de
// SPDX-License-Identifier: LGPL-3.0-or-later

#pragma once

#include <ikarus/finiteelements/fetraits.hh>
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
        : localView_{p_basis.localView()} {
      static_assert(Ikarus::Concepts::PowerBasis<RootBasis>,
                    "You didn't pass a localview of a power basis to this method");
      static_assert(RootBasis::PreBasis::Node::degree() != 1,
                    "The basis has only one children. Maybe use scalarFE.hh. ");

      localView_.bind(element);
    }

    /** \brief Type of the Pairs of gridEntities and variable tags */
    using GridElementEntityType = typename LocalView::Element;
    using Traits                = FETraits<GridElementEntityType>;

    /** \brief Number of children in the powerBasis */
    static constexpr int num_children = RootBasis::PreBasis::Node::degree();

    [[nodiscard]] constexpr size_t size() const { return localView_.size(); }

    void globalFlatIndices(std::vector<GlobalIndex>& globalIndices) const {
      globalIndices.clear();
      static_assert(
          requires { localView_.tree().child(0); },
          "Your basis does not provide a child accessor. Maybe use scalarFE.hh.");
      const auto& fe = localView_.tree().child(0).finiteElement();
      for (size_t i = 0; i < fe.size(); ++i) {
        for (int j = 0; j < num_children; ++j) {
          globalIndices.push_back(localView_.index((localView_.tree().child(j).localIndex(i))));
        }
      }
    }

    const GridElementEntityType& gridElement() const { return localView_.element(); }
    const LocalView& localView() const { return localView_; }
    LocalView& localView() { return localView_; }

  private:
    LocalView localView_;
  };
}  // namespace Ikarus
