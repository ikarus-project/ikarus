//
// Created by Alex on 28.09.2021.
//

#pragma once
#include <ikarus/Grids/GridViews/SimpleGridView.h>

template <class ContainerType>
class FEIndexSet {
public:
  FEIndexSet(const ContainerType &feContainer) {
    for (size_t index = 0, subIndex = 0; auto &fe : feContainer) {
      indexContainer_[&fe] = index++;
      auto eleDofs         = Ikarus::FiniteElements::getEntityVariableTuple(fe);
      for (auto &&[entityID, entityType, dofTypes] : eleDofs)
        if (subIndexContainer_.find(entityID) == subIndexContainer_.end()) subIndexContainer_[entityID] = subIndex++;
    }
  }

  using FEType = typename ContainerType::value_type;

  size_t index(const FEType &fe) { return indexContainer_[&fe]; }

  auto variableIndices(const FEType &fe) const {
    auto getSubIndex = [&subIndexCont = subIndexContainer_](const auto &entityIDDof) {
      return subIndexCont.at(entityIDDof.entityID);
    };
    auto entityIDDofRange = Ikarus::FiniteElements::getEntityVariableTuple(fe);
    auto res              = std::ranges::transform_view(entityIDDofRange, getSubIndex);
    auto res2             = std::vector<size_t>(res.begin(), res.end());
    return res2;
  }

  auto indexOfEntity(const size_t id) const { return subIndexContainer_.at(id); }

  int size() { return indexContainer_.size(); }

private:
  std::map<const FEType *, size_t> indexContainer_;
  std::unordered_map<size_t, size_t> subIndexContainer_;
};
