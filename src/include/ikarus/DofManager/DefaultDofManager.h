//
// Created by Alex on 26.05.2021.
//

#pragma once

#include <ranges>
#include <unordered_set>

#include <dune/common/hybridutilities.hh>
#include <dune/geometry/dimension.hh>

#include <ikarus/FiniteElements/FiniteElementPolicies.h>
#include <ikarus/FiniteElements/InterfaceFiniteElement.h>
#include <ikarus/Variables/InterfaceVariable.h>
#include <ikarus/Variables/VariableDefinitions.h>
#include <ikarus/utils/std/algorithms.h>
#include <ikarus/utils/std/hashs.h>

namespace Ikarus::DofManager {

  class VariableVector {
  public:
    size_t correctionSize() {
      if (!dofSizeCalculated) {
        dofSizeValue      = std::accumulate(begin(variableIndices), end(variableIndices), 0,
                                       [](size_t acc, auto&& varVecLoc) { return acc + varVecLoc.size(); });
        dofSizeCalculated = true;
      }
      return dofSizeValue;
    }

    auto getValues() { return std::ranges::join_view(dofVecimpl); }

  private:
    std::vector<std::vector<Ikarus::Variable::IVariable>> dofVecimpl;
    std::vector<Eigen::ArrayXi> variableIndices;
    size_t dofSizeValue;
    bool dofSizeCalculated{false};

    template <typename FEContainer, typename GridViewType>
    requires Concepts::HasgetEntityVariablePairs<typename std::remove_cvref_t<FEContainer>::value_type> || Concepts::
        HasFreegetEntityVariablePairs<typename std::remove_cvref_t<FEContainer>::value_type>
    friend class DefaultDofManager;

    friend VariableVector& operator+=(VariableVector& varVecArg, Eigen::VectorXd& correction);
  };

  VariableVector& operator+=(VariableVector& varVecArg, Eigen::VectorXd& correction) ;


  template <typename FEContainer, typename GridViewType>
  requires Concepts::HasgetEntityVariablePairs<typename std::remove_cvref_t<FEContainer>::value_type> || Concepts::
      HasFreegetEntityVariablePairs<typename std::remove_cvref_t<FEContainer>::value_type>
  class DefaultDofManager {
  public:
    using DofSet        = std::unordered_map<size_t, std::unordered_set<Ikarus::Variable::VariablesTags>>;
    using IndexMap      = std::unordered_map<size_t, size_t>;
    using DofVector     = std::vector<std::vector<Ikarus::Variable::IVariable>>;
    using VariableIndex = std::vector<Eigen::ArrayXi>;

    DefaultDofManager(FEContainer& feContainer, GridViewType& gv) : fe_container_{&feContainer}, gridView_{&gv} {}

    void createElementDofRelationship() {
      size_t indexCounter = 0;
      variableIndexMap.clear();
      indexSet.clear();
      varVec.dofVecimpl.clear();
      varVec.dofSizeCalculated = false;
      DofSet dofSet;
      for (auto&& fe : (*fe_container_)) {
        gridEntityFEmap.insert({getEntityID(fe), &fe});
        auto eleDofs = getEntityVariablePairs(fe);  // get list of entity id and list of variable tags
        for (auto&& [entityID, dofTagVector] : eleDofs)
          // Create set of Dofs for each grid entity
          dofSet[entityID].insert(begin(dofTagVector), end(dofTagVector));
      }
      std::vector<std::pair<size_t, std::unordered_set<Ikarus::Variable::VariablesTags>>> dofVector2(begin(dofSet),
                                                                                                     end(dofSet));
      std::ranges::sort(dofVector2, [](auto&& a, auto&& b) { return a.first < b.first; });

      for (auto&& [entityID, dofTagVector] : dofVector2) {
        variableIndexMap[entityID] = indexCounter++;

        varVec.dofVecimpl.emplace_back();
        for (auto&& var : dofTagVector | std::views::transform(&Ikarus::Variable::VariableFactory::createVariable))
          varVec.dofVecimpl.back().push_back(var);
      }

      size_t indexCounter2 = 0;
      varVec.variableIndices.resize(varVec.dofVecimpl.size());
      for (int i = 0; auto&& var : varVec.variableIndices) {
        var.resize(Ikarus::Variable::correctionSize(varVec.dofVecimpl[i++]));
        std::iota(var.begin(), var.end(), indexCounter2);
        indexCounter2 += var.size();
      }
      isElementDofRelationshipCreated = true;
    }

    size_t correctionSize() {
      if (isElementDofRelationshipCreated == true)
        return varVec.correctionSize();
      else
        throw std::logic_error("You first need to call createElementDofRelationship before this");
    }

  public:
    template <typename GridEntityType>
    size_t elementDofVectorSize(GridEntityType& ele) {
      size_t dofsize = 0;

      auto dofsizeIncrementFromEntity = [this](auto&& ent) {
        auto entID = ent->getID();
        if (hasIndex(entID))
          return Ikarus::Variable::correctionSize(varVec.dofVecimpl.at(variableIndexMap.at(entID)));
        else
          return size_t{0};
      };
      constexpr size_t dim = GridEntityType::dimension;
      // Loop over all entities of the grid and collect the sizes of their degrees of freedom
      Dune::Hybrid::forEach(Dune::Hybrid::integralRange(Dune::index_constant<dim>()), [&](auto i) {
        for (auto&& ent : entities(ele, i))
          dofsize += dofsizeIncrementFromEntity(ent);
        dofsize += dofsizeIncrementFromEntity(&ele);
      });
      return dofsize;
    }

    template <typename GridEntityType>
    auto elementDofIndices(GridEntityType& ele) {
      Eigen::ArrayXi indices(elementDofVectorSize(ele));
      size_t posHelper = 0;

      auto dofIndicesFromEntity = [this, &posHelper, &indices](auto&& ent) {
        auto entID = ent->getID();
        if (hasIndex(entID)) {
          auto sizeOfVar                                 = varVec.variableIndices[variableIndexMap.at(entID)].size();
          indices.template segment(posHelper, sizeOfVar) = varVec.variableIndices[variableIndexMap.at(entID)];
          posHelper += sizeOfVar;
        }
      };

      constexpr size_t dim = GridEntityType::dimension;
      // Loop over all entities of the grid and collect the indices of their degrees of freedom
      Dune::Hybrid::forEach(Dune::Hybrid::integralRange(Dune::index_constant<dim>()), [&](auto i) {
        for (auto&& ent : entities(ele, i))
          dofIndicesFromEntity(ent);
        dofIndicesFromEntity(&ele);
      });

      return indices;
    }

    template <typename GridEntityType>
    auto elementVariables(GridEntityType& ele) {
      std::vector<Ikarus::Variable::IVariable*> elementVariables;

      auto variablesOfEntity = [this, &elementVariables](auto&& ent) {
        auto entID = ent->getID();
        if (hasIndex(entID)) {
          auto& entityDofVector = varVec.dofVecimpl.at(variableIndexMap.at(entID));
          auto pointerRange     = Ikarus::stl::transformValueRangeToPointerRange(entityDofVector);
          elementVariables.insert(elementVariables.end(), pointerRange.begin(), pointerRange.end());
        }
      };

      constexpr size_t dim = GridEntityType::dimension;
      Dune::Hybrid::forEach(Dune::Hybrid::integralRange(Dune::index_constant<dim>()), [&](auto i) {
        for (auto&& ent : entities(ele, i))
          variablesOfEntity(ent);
        variablesOfEntity(&ele);
      });

      return elementVariables;
    }

    auto elementDofsVariableTuple() {
      constexpr int gridDim = GridViewType::dimension;
      return std::ranges::transform_view(entities(gridView_, Dune::index_constant<gridDim>()), [&](auto&& ge) {
        auto elementsDegreesOfFreedom = elementDofIndices(ge);
        auto elementsVariables        = elementVariables(ge);
        return std::make_tuple(gridEntityFEmap[ge.getID()], elementsDegreesOfFreedom, elementsVariables);
      });
    }

    auto& getVariables() { return varVec; }

  private:
    bool hasIndex(size_t iD) {
      try {
        variableIndexMap.at(iD);
        return true;
      } catch (const std::out_of_range&) {
        return false;
      }
    }
    FEContainer* fe_container_;
    GridViewType* gridView_;
    bool isElementDofRelationshipCreated{false};
    IndexMap variableIndexMap;
    std::unordered_map<size_t, Ikarus::FiniteElements::IFiniteElement*> gridEntityFEmap;
    std::vector<std::pair<size_t, std::vector<size_t>>> indexSet;
    VariableVector varVec;
  };

}  // namespace Ikarus::DofManager