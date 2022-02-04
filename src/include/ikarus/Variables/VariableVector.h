//
// Created by Alex on 28.09.2021.
//

#pragma once

#include <cstdarg>
#include <iostream>
#include <unordered_set>

#include <ikarus/FiniteElements/FEIndexSet.h>
#include <ikarus/FiniteElements/FiniteElementFunctionConcepts.h>
#include <ikarus/utils/utils/algorithms.h>

namespace Ikarus::Variable {

  template <class FEContainer>
  class VariableVector {
  public:
    explicit VariableVector(const FEContainer &feContainer) : feContainer_{&feContainer} {
      std::cout << "VariableVector1" << std::endl;
      using DofSet = std::map<Ikarus::FiniteElements::VariableIndicesPair::Indices,
                              std::unordered_set<Ikarus::Variable::VariableTags>>;
      DofSet dofSet;
      for (auto &&fe : feContainer) {
        for (auto &&[indices, dofTypes] : Ikarus::FiniteElements::getEntityVariableTuple(fe)) {
          auto &&entitySetEntry = dofSet[indices];
          entitySetEntry.unorderedSet.insert(begin(dofTypes), end(dofTypes));
        }
      }

      for (int pos = 0; auto &[indices, varTagsSet] : dofSet) {
        auto viewOverTags = varTagsSet | std::views::transform([](auto &&tag) {
                              return Ikarus::Variable::VariableFactory::createVariable(tag);
                            });
        for (auto var : viewOverTags) {
          variablesForEachNode.push_back(std::move(var));
          auto &currentVar = variablesForEachNode.back();
          variableID.insert({&currentVar, variablesForEachNode.size() - 1});
          auto corsize = correctionSize(currentVar);
          Eigen::Matrix<size_t, Eigen::Dynamic, 1, 8, 1> indicesOfvariable;
          for (int i = 0; i < corsize; ++i)
            indicesOfvariable(i) = indices[i + pos];
          pos += corsize;
          variableIndices.insert({&currentVar, indicesOfvariable});
          variableOfIndices.insert({indicesOfvariable,&currentVar });
        }
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

    auto getValues() { return variablesForEachNode; };
    auto getValues() const { return variablesForEachNode; };
    auto size() const {
      int size = std::accumulate(variablesForEachNode.begin(), variablesForEachNode.end(), 0,
                                 [](auto &&acc, auto &&var) { return acc + valueSize(var); });
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
      std::vector<Ikarus::Variable::IVariable*> elementVariables;

      for (auto &&[indices, dofTypes] : Ikarus::FiniteElements::getEntityVariableTuple(fe)) {

      }
      return elementVariables;
    }

    auto variablesOfSingleElement(const typename FEContainer::value_type &fe) const {
      FiniteElements::FEValues elementVariables;

      //      auto feSubIndicesRange = feIndexSet.variableIDs(fe);
      //      for (auto &feSubIndex : feSubIndicesRange) {
      //        const auto &entityType = entityTypes.at(feSubIndex);
      //        elementVariables.add(entityType, variablesForEachNode[variablesForEachEntityfeSubIndex]);
      //      }
      return elementVariables;
    }

    auto correctionSize() const {
      int size = std::accumulate(variablesForEachNode.begin(), variablesForEachNode.end(), 0,
                                 [](auto &&acc, auto &&var) { return acc + Variable::correctionSize(var); });
      return size;
    };
    //    [[nodiscard]] size_t correctionSize() const { return dofSizeValue; }
    [[nodiscard]] size_t elementSize() const { return feContainer_->size(); }
    const auto &getFeContainer() const { return *feContainer_; }

    VariableVector &operator+=(const Eigen::VectorXd &correction) {
      if (static_cast<decltype(correction.size())>(correctionSize()) == correction.size())
        for (size_t variableIndex = 0; auto &&var : std::ranges::join_view(this->variablesForEachNode))
          var += correction(variableIndices[variableIndex++]);
      return *this;
    }

    VariableVector &operator-=(const Eigen::VectorXd &correction) {
      assert(static_cast<long long int>(correctionSize()) == correction.size());
      for (size_t variableIndex = 0; auto &&var : std::ranges::join_view(this->variablesForEachNode))
        var -= correction(variableIndices[variableIndex++]);
      return *this;
    }

    template <typename Range>
    requires(!std::is_same_v<Range, Eigen::VectorXd>) VariableVector &operator+=(Range &&r) {
      assert(static_cast<decltype(r.size())>(correctionSize()) == r.size());
      for (size_t variableIndex = 0; auto &&var : std::ranges::join_view(this->variablesForEachNode)) {
        const auto &variableIndices_ = variableIndices[variableIndex++];
        for (int i = 0; i < variableIndices_.size(); ++i)
          var[i] += r[variableIndices_[i]];
      }
      return *this;
    }

    template <typename Range>
    requires(!std::is_same_v<Range, Eigen::VectorXd>) VariableVector &operator-=(Range &&r) {
      assert(static_cast<decltype(r.size())>(correctionSize()) == r.size());
      for (size_t variableIndex = 0; auto &&var : std::ranges::join_view(this->variablesForEachNode)) {
        const auto &variableIndices_ = variableIndices[variableIndex++];
        for (int i = 0; i < variableIndices_.size(); ++i)
          var[i] -= r[variableIndices_[i]];
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
    friend VariableVector<FEContainer1> operator+(const VariableVector<FEContainer1> &varVec,
                                                  const Eigen::VectorXd &correction);

    mutable std::vector<Ikarus::Variable::IVariable> variablesForEachNode;
    template <typename FEContainer1, typename Range>
    requires(!std::is_same_v<Range, Eigen::VectorXd>) friend VariableVector<FEContainer1> operator+(
        const VariableVector<FEContainer1> &varVec, Range &&r);
    std::map<Variable::IVariable *, Eigen::Matrix<size_t, Eigen::Dynamic, 1, 8, 1>> variableIndices;
    std::map<Eigen::Matrix<size_t, Eigen::Dynamic, 1, 8, 1>,Variable::IVariable * > variableOfIndices;
    std::map<Variable::IVariable *, size_t> variableID;
    FEContainer const *feContainer_;
  };

  template <class FEContainer>
  VariableVector<FEContainer> operator+(const VariableVector<FEContainer> &varVec, const Eigen::VectorXd &correction) {
    VariableVector res = varVec;
    assert(static_cast<decltype(correction.size())>(res.correctionSize()) == correction.size());
    for (size_t variableIndex = 0; auto &&var : std::ranges::join_view(res.variablesForEachNode))
      var += correction(res.variableIndices[variableIndex++]);
    return res;
  }

  template <typename FEContainer, typename Range>
  requires(!std::is_same_v<Range, Eigen::VectorXd>) VariableVector<FEContainer>
  operator+(const VariableVector<FEContainer> &varVec, Range &&r) {
    VariableVector res = varVec;
    assert(static_cast<decltype(r.size())>(varVec.correctionSize()) == r.size());
    for (size_t variableIndex = 0; auto &&var : std::ranges::join_view(res.variablesForEachNode)) {
      const auto &variableIndices_ = res.variableIndices[variableIndex++];
      for (int i = 0; i < variableIndices_.size(); ++i)
        var[i] += r[variableIndices_[i]];
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