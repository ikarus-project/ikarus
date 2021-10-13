
#pragma once

#include <dune/common/hybridutilities.hh>

#include <ikarus/Grids/EntityHelperFunctions.h>
#include <ikarus/Variables/VariableDefinitions.h>
#include <ikarus/utils/utils/algorithms.h>
#include <ikarus/utils/utils/traits.h>

template <typename... Args>
struct At {
  std::tuple<Args...> args;
};

template <typename... Args>
auto at(Args &&...args) {
  return At<Args &&...>{std::forward_as_tuple(std::forward<Args>(args)...)};
}

template <typename... Args>
struct Data {
  std::tuple<Args...> args;
};

template <typename... Args>
auto data(Args &&...args) {
  return Data<Args &&...>{std::forward_as_tuple(std::forward<Args>(args)...)};
}

template <typename... Vars>
auto datatuple(const Vars &...vars) requires
    std::conjunction<std::is_convertible<Vars, Ikarus::Variable::VariablesTags>...>::value {
  return std::array<Ikarus::Variable::VariablesTags, sizeof...(Vars)>({vars...});
}

template <typename IndexSetType>
class GridData {
public:
  static constexpr int gridDim      = IndexSetType::dimension;
  static constexpr int gridWorldDim = IndexSetType::worldDimension;

  explicit GridData(const IndexSetType &indexSet) : indexSet_{&indexSet} {}

  void resizeIfEmpty(int entityTypeId) {
    if (variablesForEachEntity[entityTypeId].empty())
      variablesForEachEntity[entityTypeId].resize(indexSet_->size(gridDim - entityTypeId));
  }

  template <typename... DataArgs, typename... Args>
  void add(const Data<DataArgs...> &data, const At<Args...> &at) {
    static_assert(sizeof...(DataArgs) >= 1);

    constexpr auto n = sizeof...(Args);
    static_assert(n >= 1);
    static_assert(n == sizeof...(DataArgs), "The entity types and datatype size should be the same!");

    Dune::Hybrid::forEach(Dune::Hybrid::integralRange(Dune::index_constant<n>()), [&](const auto i) {
      static constexpr int entityTypeId = Ikarus::resolveEntityType<decltype(std::get<i>(at.args)), gridDim>();
      resizeIfEmpty(entityTypeId);

      for (auto &vectorOfEntityVariables : variablesForEachEntity[entityTypeId]) {
        const auto entityData = std::get<i>(data.args);

        if constexpr (Ikarus::utils::is_std_array<std::decay_t<decltype(entityData)>>::value)
          for (auto singleVar : entityData | std::views::transform(&Ikarus::Variable::VariableFactory::createVariable))
            Ikarus::utils::appendUnique(vectorOfEntityVariables, singleVar);
        else
          Ikarus::utils::appendUnique(vectorOfEntityVariables,
                                      Ikarus::Variable::VariableFactory::createVariable(entityData));
      }
    });
  }

  template <typename... DataArgs, typename... Args>
  void remove(const Data<DataArgs...> &data, const At<Args...> &at) {
    static_assert(sizeof...(DataArgs) >= 1);

    constexpr auto n = sizeof...(Args);
    static_assert(n >= 1);
    static_assert(n == sizeof...(DataArgs), "The entity types and datatype size should be the same!");

    Dune::Hybrid::forEach(Dune::Hybrid::integralRange(Dune::index_constant<n>()), [&](const auto i) {
      static constexpr int entityTypeId = Ikarus::resolveEntityType<decltype(std::get<i>(at.args)), gridDim>();
      resizeIfEmpty(entityTypeId);

      for (auto &vectorOfEntityVariables : variablesForEachEntity[entityTypeId]) {
        const auto entityData = std::get<i>(data.args);

        if constexpr (Ikarus::utils::is_std_array<std::decay_t<decltype(entityData)>>::value)
          for (auto singleVar : entityData | std::views::transform(&Ikarus::Variable::VariableFactory::createVariable))
            vectorOfEntityVariables.erase(
                std::remove(vectorOfEntityVariables.begin(), vectorOfEntityVariables.end(), singleVar),
                vectorOfEntityVariables.end());
        else
          vectorOfEntityVariables.erase(std::remove(vectorOfEntityVariables.begin(), vectorOfEntityVariables.end(),
                                                    Ikarus::Variable::VariableFactory::createVariable(entityData)),
                                        vectorOfEntityVariables.end());
      }
    });
  }

  template <typename EntityType>
  auto &getData(const EntityType &entity) {
    static constexpr int entityDim = EntityType::mydimension;
    return variablesForEachEntity[entityDim][indexSet_->index(entity)];
  }

  template <typename EntityType>
  auto getAllSubEntityData(const EntityType &entity) {
    std::vector<Ikarus::Variable::IVariable *> varVector;

    for (int dimSubEntity = 0; dimSubEntity < EntityType::mydimension; ++dimSubEntity)
      for (size_t i = 0; i < entity.subEntities(gridWorldDim - dimSubEntity); ++i) {
        auto subEntityIndex = indexSet_->subIndex(entity, i, gridWorldDim - dimSubEntity);
        for (auto &var : variablesForEachEntity[dimSubEntity][subEntityIndex])
          varVector.push_back(&var);
      }

    for (auto &var : getData(entity))  // add vars of the entity itself
      varVector.push_back(&var);
    return varVector;
  }

  template <typename FEType>
  auto getAllSubEntityDataOfFE(const FEType &fe) {
    std::vector<Ikarus::Variable::IVariable *> varVector;

    for (unsigned int dimSubEntity = 0; dimSubEntity < dimension(fe); ++dimSubEntity)
      for (size_t i = 0; i < subEntities(fe, gridWorldDim - dimSubEntity); ++i) {
        auto subEntityIndex = subIndex(fe, i, gridWorldDim - dimSubEntity);
        for (auto &var : variablesForEachEntity[dimSubEntity][subEntityIndex])
          varVector.push_back(&var);
      }

    for (auto &var : variablesForEachEntity[dimension(fe)][subIndex(
             fe, 0, gridWorldDim - dimension(fe))])  // add vars of the entity itself
      varVector.push_back(&var);
    return varVector;
  }

private:
  std::array<std::vector<std::vector<Ikarus::Variable::IVariable>>, gridDim + 1> variablesForEachEntity;
  IndexSetType const *const indexSet_;
};