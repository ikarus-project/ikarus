
#pragma once

#include <dune/common/hybridutilities.hh>

#include <ikarus/Grids/EntityHelperFunctions.h>
#include <ikarus/Variables/VariableDefinitions.h>
#include <ikarus/utils/utils/algorithms.h>
#include <ikarus/utils/utils/traits.h>

template <typename... Args>
struct Data {
  std::tuple<Args...> args;
};

template <typename... Args>
auto data(Args &&...args) {
  return Data<Args &&...>{std::forward_as_tuple(std::forward<Args>(args)...)};
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

  template <typename... DataArgs, int codim>
  requires(sizeof...(DataArgs) >= 1) void add(const Data<DataArgs...> &data, Ikarus::EntitiesWithCoDim<codim> &&) {
    add(data, Ikarus::getEntityTypeFromCodim(gridDim, codim));
  }

  template <typename... DataArgs>
  requires(sizeof...(DataArgs) >= 1) void add(const Data<DataArgs...> &data, const Ikarus::EntityType &at) {
    constexpr auto n = sizeof...(DataArgs);
    Dune::Hybrid::forEach(Dune::Hybrid::integralRange(Dune::index_constant<n>()), [&](const auto i) {
      const int entityTypeId = Ikarus::determineEntityDimension(at);
      resizeIfEmpty(entityTypeId);

      for (auto &vectorOfEntityVariables : variablesForEachEntity[entityTypeId]) {
        const auto entityData = std::get<i>(data.args);
        Ikarus::utils::appendUnique(vectorOfEntityVariables,
                                    Ikarus::Variable::VariableFactory::createVariable(entityData));
      }
    });
  }

  template <typename... DataArgs>
  requires(sizeof...(DataArgs) >= 1) void remove(const Data<DataArgs...> &data, const Ikarus::EntityType &at) {
    constexpr auto n = sizeof...(DataArgs);
    Dune::Hybrid::forEach(Dune::Hybrid::integralRange(Dune::index_constant<n>()), [&](const auto i) {
      const int entityTypeId = Ikarus::determineEntityDimension(at);
      resizeIfEmpty(entityTypeId);

      for (auto &vectorOfEntityVariables : variablesForEachEntity[entityTypeId]) {
        const auto entityData = std::get<i>(data.args);
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
    Ikarus::FiniteElements::FEValues entityValues;

    for (int dimSubEntity = 0; dimSubEntity < EntityType::mydimension; ++dimSubEntity)
      for (size_t i = 0; i < entity.subEntities(gridWorldDim - dimSubEntity); ++i) {
        auto subEntityIndex = indexSet_->subIndex(entity, i, gridWorldDim - dimSubEntity);
        for (auto &var : variablesForEachEntity[dimSubEntity][subEntityIndex])
          entityValues.add(Ikarus::getEntityTypeFromCodim(gridDim, gridDim - dimSubEntity), var);
      }

    for (auto &var : getData(entity))  // add vars of the entity itself
      entityValues.add(Ikarus::getEntityTypeFromCodim(gridDim, 0), var);
    return entityValues;
  }

  template <typename FEType>
  auto getAllSubEntityDataOfFE(const FEType &fe) {
    Ikarus::FiniteElements::FEValues entityValues;

    for (unsigned int dimSubEntity = 0; dimSubEntity < dimension(fe); ++dimSubEntity)
      if (!variablesForEachEntity[dimSubEntity].empty())
        for (size_t i = 0; i < subEntities(fe, gridWorldDim - dimSubEntity); ++i) {
          auto subEntityIndex = subIndex(fe, i, gridWorldDim - dimSubEntity);
          for (auto &var : variablesForEachEntity[dimSubEntity][subEntityIndex])
            entityValues.add(Ikarus::getEntityTypeFromCodim(gridDim, gridDim - dimSubEntity), var);
        }

    if (!variablesForEachEntity[dimension(fe)].empty())
      for (auto &var : variablesForEachEntity[dimension(fe)][subIndex(
               fe, 0, gridWorldDim - dimension(fe))])  // add vars of the entity itself
        entityValues.add(Ikarus::getEntityTypeFromCodim(gridDim, 0), var);
    return entityValues;
  }

private:
  std::array<std::vector<std::vector<Ikarus::Variable::IVariable>>, gridDim + 1> variablesForEachEntity;
  IndexSetType const *const indexSet_;
};