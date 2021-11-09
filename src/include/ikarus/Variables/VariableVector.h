//
// Created by Alex on 28.09.2021.
//

#pragma once

#include <cstdarg>
#include <iostream>
#include <unordered_set>

#include <ikarus/FiniteElements/FEIndexSet.h>
#include <ikarus/utils/utils/algorithms.h>

namespace Ikarus::Variable {

  template <class FEContainer>
  class VariableVector {
  public:
    explicit VariableVector(const FEContainer &feContainer) : feContainer_{&feContainer}, feIndexSet{feContainer} {
      struct EntityTypeAndVarTagSet {
        Ikarus::EntityType type;
        std::unordered_set<Ikarus::Variable::VariablesTags> unorderedSet;
      };
      using DofSet = std::unordered_map<size_t, EntityTypeAndVarTagSet>;
      DofSet dofSet;
      for (auto &&fe : feContainer)
        for (auto &&[entityID, entityType, dofTypes] :
             Ikarus::FiniteElements::getEntityVariableTuple(fe))  // Create set of Dofs for each grid entity
        {
          auto &&entitySetEntry = dofSet[entityID];
          entitySetEntry.type   = entityType.value();
          entitySetEntry.unorderedSet.insert(begin(dofTypes), end(dofTypes));
        }

      // A vector of pairs of a entity id and the corresponding unique variable tags
      using EntityVariablePairVector = std::vector<std::pair<size_t, EntityTypeAndVarTagSet>>;

      EntityVariablePairVector dofVector(begin(dofSet), end(dofSet));
      // sort vector to have increasing ids of entities
      std::ranges::sort(dofVector, [](auto &&a, auto &&b) { return a.first < b.first; });

      variablesForEachEntity.resize(dofVector.size());
      // creating vector of variables and save in variableIndexMap to which entity the belong
      for (auto &&[entityID, typeAnddofTagVector] : dofVector) {
        auto &&[entityType, dofTagVector] = typeAnddofTagVector;
        const auto &entityIndex           = feIndexSet.indexOfEntity(entityID);
        entityTypes[entityIndex]          = entityType;
        for (auto &&var : dofTagVector | std::views::transform(&Ikarus::Variable::VariableFactory::createVariable))
          variablesForEachEntity[entityIndex].emplace_back(var);
      }

      // find out how much dofs we actually have
      for (auto &&var : getValues())
        dofSizeValue += Ikarus::Variable::correctionSize(var);

      size_t indexCounter = 0;
      // add increasing degrees of freedom integer indices
      variableIndices.resize(dofVector.size());
      for (int i = 0; auto &&var : variableIndices) {
        var.resize(Ikarus::Variable::correctionSize(this->variablesForEachEntity[i++]));
        std::iota(var.begin(), var.end(), indexCounter);
        indexCounter += var.size();
      }
    }

    auto elementVariables() {
      return transform_viewOverElements([this](const auto &fe) { return this->variablesOfSingleElement(fe); });
    };

    auto elementDofVectorSize() {
      return transform_viewOverElements([this](const auto &fe) { return dofSize(fe); });
    };

    auto elementDofs() {
      return transform_viewOverElements([this](const auto &fe) { return this->dofIndicesOfSingleElement(fe); });
    };

    auto getValues() { return std::ranges::join_view(variablesForEachEntity); };

    auto variablesOfSingleElement(const typename FEContainer::value_type &fe) {
      Ikarus::FEValues elementVariables;

      auto feSubIndicesRange = feIndexSet.subIndices(fe);
      for (auto &feSubIndex : feSubIndicesRange) {
        const auto &entityType = entityTypes.at(feSubIndex);
        elementVariables.set(entityType, variablesForEachEntity[feSubIndex]);
      }
      return elementVariables;
    }

    auto dofIndicesOfSingleElement(const typename FEContainer::value_type &fe) {
      Eigen::ArrayXi indices(dofSize(fe));
      indices.setZero();
      size_t posHelper = 0;

      auto eleDofs = Ikarus::FiniteElements::getEntityVariableTuple(fe);

      for (auto &&entityID : feIndexSet.subIndices(fe)) {
        auto &currentIndices = variableIndices[feIndexSet.indexOfEntity(entityID)];
        indices.template segment(posHelper, currentIndices.size()) = currentIndices;
        posHelper += currentIndices.size();
      }
      return indices;
    }

    size_t correctionSize() { return dofSizeValue; }

    VariableVector &operator+=(const Eigen::VectorXd &correction) {
      assert(static_cast<long long int>(correctionSize()) == correction.size());
      for (size_t variableIndex = 0; auto &&var : std::ranges::join_view(this->variablesForEachEntity))
        var += correction(variableIndices[variableIndex++]);
      return *this;
    }
    template <class Functype>
    auto transform_viewOverElements(Functype fn) {
      return std::ranges::transform_view(*feContainer_, fn);
    }

  private:
    std::vector<std::vector<Ikarus::Variable::IVariable>> variablesForEachEntity;
    std::unordered_map<size_t, Ikarus::EntityType> entityTypes;
    std::vector<Eigen::ArrayXi> variableIndices;
    size_t dofSizeValue{};
    FEContainer const *feContainer_;
    FEIndexSet<FEContainer> feIndexSet;
  };
}  // namespace Ikarus::Variable