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
        std::unordered_set<Ikarus::Variable::VariableTags> unorderedSet;
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
      for (auto &&[entityID, typeAnddofTagSet] : dofVector) {
        auto &&[entityType, dofTagSet] = typeAnddofTagSet;
        const auto &entityIndex        = feIndexSet.indexOfEntity(entityID);
        entityTypes[entityIndex]       = entityType;
        for (auto &&var : dofTagSet | std::views::transform(&Ikarus::Variable::VariableFactory::createVariable))
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

    auto elementVariables() const {
      return transform_viewOverElements([this](const auto &fe) { return this->variablesOfSingleElement(fe); });
    };

    auto elementDofVectorSize() {
      return transform_viewOverElements([this](const auto &fe) { return dofSize(fe); });
    };

    auto elementDofs() {
      return transform_viewOverElements([this](const auto &fe) { return this->dofIndicesOfElement(fe); });
    };

    auto getValues() { return std::ranges::join_view(variablesForEachEntity); };
    auto getValues() const { return std::ranges::join_view(variablesForEachEntity); };
    auto size() const {
      auto v   = std::ranges::join_view(variablesForEachEntity);
      int size = std::accumulate(v.begin(), v.end(), 0, [](auto &&acc, auto &&var) { return acc + valueSize(var); });
      return size;
    };

    [[nodiscard]] Eigen::VectorXd getUnderlyingVector() const {
      Eigen::VectorXd vec(size());
      for (int offset = 0; auto &&v : this->getValues()) {
        const auto curvar                  = getValue(v);
        vec.segment(offset, curvar.size()) = curvar;
        offset += curvar.size();
      }
      return vec;
    }

    auto variablesOfSingleElement(const typename FEContainer::value_type &fe) {
      FiniteElements::FEValues elementVariables;

      auto feSubIndicesRange = feIndexSet.variableIDs(fe);
      for (auto &feSubIndex : feSubIndicesRange) {
        const auto &entityType = entityTypes.at(feSubIndex);
        elementVariables.add(entityType, variablesForEachEntity[feSubIndex]);
      }
      return elementVariables;
    }

    auto variablesOfSingleElement(const typename FEContainer::value_type &fe) const {
      FiniteElements::FEValues elementVariables;

      auto feSubIndicesRange = feIndexSet.variableIDs(fe);
      for (auto &feSubIndex : feSubIndicesRange) {
        const auto &entityType = entityTypes.at(feSubIndex);
        elementVariables.add(entityType, variablesForEachEntity[feSubIndex]);
      }
      return elementVariables;
    }

    auto dofIndicesOfElement(const typename FEContainer::value_type &fe) const {
      Eigen::ArrayX<size_t> indices(dofSize(fe));
      indices.setZero();
      size_t posHelper = 0;

      auto eleDofs = Ikarus::FiniteElements::getEntityVariableTuple(fe);

      for (auto &&entityID : feIndexSet.variableIDs(fe)) {
        auto &currentIndices                                       = variableIndices[entityID];
        indices.template segment(posHelper, currentIndices.size()) = currentIndices;
        posHelper += currentIndices.size();
      }
      return indices;
    }

    const auto &dofIndicesOfEntity(const size_t &gridIndexOfEntity) const { return variableIndices[gridIndexOfEntity]; }

    [[nodiscard]] size_t correctionSize() const { return dofSizeValue; }
    [[nodiscard]] size_t elementSize() const { return feContainer_->size(); }
    const auto &getFeContainer() const { return *feContainer_; }

    VariableVector &operator+=(const Eigen::VectorXd &correction) {
      if (static_cast<decltype(correction.size())>(correctionSize()) == correction.size())
      for (size_t variableIndex = 0; auto &&var : std::ranges::join_view(this->variablesForEachEntity))
        var += correction(variableIndices[variableIndex++]);
      return *this;
    }


    VariableVector &operator-=(const Eigen::VectorXd &correction) {
      assert(static_cast<long long int>(correctionSize()) == correction.size());
      for (size_t variableIndex = 0; auto &&var : std::ranges::join_view(this->variablesForEachEntity))
        var -= correction(variableIndices[variableIndex++]);
      return *this;
    }

    template<typename Range> requires (!std::is_same_v<Range,Eigen::VectorXd>)
    VariableVector &operator+=(Range&& r) {
      assert(static_cast<decltype(r.size())>(correctionSize()) == r.size());
      for (size_t variableIndex = 0; auto &&var : std::ranges::join_view(this->variablesForEachEntity)) {
        const auto& variableIndices_ = variableIndices[variableIndex++];
        for (int i = 0; i < variableIndices_.size(); ++i) {
          var[i]+= r[variableIndices_[i]];
        }
      }
      return *this;
    }

    template<typename Range> requires (!std::is_same_v<Range,Eigen::VectorXd>)
    VariableVector &operator-=(Range&& r) {
      assert(static_cast<decltype(r.size())>(correctionSize()) == r.size());
      for (size_t variableIndex = 0; auto &&var : std::ranges::join_view(this->variablesForEachEntity)) {
        const auto& variableIndices_ = variableIndices[variableIndex++];
        for (int i = 0; i < variableIndices_.size(); ++i) {
          var[i]-= r[variableIndices_[i]];
        }
      }
      return *this;
    }

    template <class Functype>
    auto transform_viewOverElements(Functype fn) {
      return std::ranges::transform_view(*feContainer_, fn);
    }

    template <class Functype>
    auto transform_viewOverElements(Functype fn) const {
      return std::ranges::transform_view(*feContainer_, fn);
    }

  private:
    template <class FEContainer1>
    friend VariableVector<FEContainer1> operator+(const VariableVector<FEContainer1>& varVec,const Eigen::VectorXd &correction) ;
        mutable std::vector<std::vector<Ikarus::Variable::IVariable>> variablesForEachEntity;
        template<typename FEContainer1,typename Range> requires (!std::is_same_v<Range,Eigen::VectorXd>)
            friend VariableVector<FEContainer1> operator+(const VariableVector<FEContainer1>& varVec,Range&& r);
    std::unordered_map<size_t, Ikarus::EntityType> entityTypes;
    std::vector<Eigen::ArrayX<size_t>> variableIndices;
    size_t dofSizeValue{};
    FEContainer const *feContainer_;
    FEIndexSet<FEContainer> feIndexSet;
  };

  template <class FEContainer>
  VariableVector<FEContainer> operator+(const VariableVector<FEContainer>& varVec,const Eigen::VectorXd &correction) {
    VariableVector res=varVec;
    assert (static_cast<decltype(correction.size())>(res.correctionSize()) == correction.size());
      for (size_t variableIndex = 0; auto &&var : std::ranges::join_view(res.variablesForEachEntity))
        var += correction(res.variableIndices[variableIndex++]);
    return res;
  }

  template<typename FEContainer,typename Range> requires (!std::is_same_v<Range,Eigen::VectorXd>)
      VariableVector<FEContainer> operator+(const VariableVector<FEContainer>& varVec,Range&& r) {
    VariableVector res=varVec;
    assert(static_cast<decltype(r.size())>(varVec.correctionSize()) == r.size());
    for (size_t variableIndex = 0; auto &&var : std::ranges::join_view(res.variablesForEachEntity)) {
      const auto& variableIndices_ = res.variableIndices[variableIndex++];
      for (int i = 0; i < variableIndices_.size(); ++i) {
        var[i]+= r[variableIndices_[i]];
      }
    }
    return res;
  }

  template <class FEContainer>
  inline std::ostream &operator<<(std::ostream &o, const VariableVector<FEContainer> &var) {
    for (auto &&v : var.getValues()) {
      o << getValue(v) << std::endl;
    }
    return o;
  }

}  // namespace Ikarus::Variable