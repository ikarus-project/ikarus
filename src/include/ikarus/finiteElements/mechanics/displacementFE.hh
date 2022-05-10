//
// Created by Alex on 11.04.2022.
//

#pragma once

#include <ikarus/finiteElements/interface/feTraits.hh>
namespace Ikarus {
  template <typename Basis>
  class DisplacementFE {
  public:
    using RootBasis   = Basis;
    using LocalView   = typename Basis::LocalView;
    using GlobalIndex = typename LocalView::MultiIndex;
    explicit DisplacementFE(const Basis& p_basis, const typename LocalView::Element& element)
        : localView{p_basis.localView()} {
      static_assert(Ikarus::Concepts::PowerBasis<RootBasis>,
                    "You didn't pass a localview of a power basis to this method");
      static_assert(RootBasis::PreBasis::Node::CHILDREN == worlddim,
                    "The power basis children number does not coincide with the world space where the grid entity is "
                    "embedded into!");
      static_assert(
          Ikarus::Concepts::FlatIndexBasis<RootBasis>,
          "The Index merging strategy of the basis you passed has to be FlatLexicographic or FlatInterLeaved");

      localView.bind(element);
    }

    /** \brief Type of the Pairs of gridEntities and variable tags */
    using GridElementEntityType = typename LocalView::Element;
    using Traits                = FETraits<GridElementEntityType>;

    /** \brief Dimension of the world space */
    static constexpr int worlddim = Traits::worlddim;

    [[nodiscard]] constexpr int size() const { return localView.size(); }

    void globalIndices(std::vector<GlobalIndex>& globalIndices) const {
      const auto& fe = localView.tree().child(0).finiteElement();
      for (size_t i = 0; i < fe.size(); ++i) {
        for (int j = 0; j < worlddim; ++j) {
          globalIndices.push_back(localView.index((localView.tree().child(j).localIndex(i))));
        }
      }
    }

    const GridElementEntityType& getEntity() { return localView.element(); }

  private:
    LocalView localView;
  };
}  // namespace Ikarus