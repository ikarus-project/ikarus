//
// Created by Alex on 11.04.2022.
//

#pragma once

#include <ikarus/finiteElements/interface/feTraits.hh>

namespace Ikarus{
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
}
